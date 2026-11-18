"""Retrieve or check the initial hourly point requests in the local weather plan."""
import argparse
import csv
import hashlib
import json
import os
import tempfile
import zipfile
from datetime import datetime, timedelta
from pathlib import Path

import netCDF4
import numpy as np


DATASETS = {"era5-land": "reanalysis-era5-land-timeseries",
            "era5": "reanalysis-era5-single-levels-timeseries"}
VARIABLES = ["2m_temperature", "2m_dewpoint_temperature"]


def validate_response(payload, row):
    def check_dataset(nc, found, expected):
        time = next((nc.variables[k] for k in ("valid_time", "time") if k in nc.variables), None)
        if time is None or not hasattr(time, "units"):
            raise ValueError("Missing hourly time axis")
        start, hours = expected
        origin, next_hour = netCDF4.date2num([start, start + timedelta(hours=1)], time.units,
                                             getattr(time, "calendar", "standard"))
        wanted = origin + np.arange(hours) * (next_hour-origin)
        if time.size != hours or not np.allclose(time[:], wanted, rtol=0, atol=abs(next_hour-origin)/3600):
            raise ValueError("Returned dates or hourly axis differ from request")
        for key in ("latitude", "longitude"):
            v = np.asarray(nc[key][:]).ravel()
            if len(v) != 1:
                raise ValueError("Expected a single grid point")
            delta = float(v[0])-float(row[key])
            if key == "longitude":
                delta = (delta+180) % 360-180
            if not np.isfinite(delta) or abs(delta) > 0.001:
                raise ValueError("Returned location differs from request")
        for short, long in zip(("t2m", "d2m"), VARIABLES):
            name = next((k for k in (short, long) if k in nc.variables), None)
            if name:
                v = nc.variables[name]
                if short in found or v.size != hours or getattr(v, "units", "").lower() not in ("k", "kelvin"):
                    raise ValueError("Unexpected hourly variable, length or unit")
                vals = np.ma.asarray(v[:]).compressed()
                vals = vals[np.isfinite(vals)]
                if np.any((vals < 180) | (vals > 340)):
                    raise ValueError("Temperature outside the validation range")
                found.add(short)

    first = datetime.fromisoformat(row["start_date"])
    hours = 24 * ((datetime.fromisoformat(row["end_date"])-first).days + 1)
    found = set()
    if zipfile.is_zipfile(payload):
        with zipfile.ZipFile(payload) as z:
            for info in z.infolist():
                if info.filename.lower().endswith(".nc"):
                    with netCDF4.Dataset("member.nc", memory=z.read(info)) as nc:
                        check_dataset(nc, found, (first, hours))
    else:
        with netCDF4.Dataset(payload) as nc:
            check_dataset(nc, found, (first, hours))
    if found != {"t2m", "d2m"}:
        raise ValueError("The weather response does not contain both variables")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dataset", choices=DATASETS)
    parser.add_argument("--download", action="store_true", help="submit missing requests to Copernicus")
    parser.add_argument("--limit", type=int, help="process only the first N windows")
    args = parser.parse_args()
    if args.limit is not None and args.limit < 1:
        parser.error("--limit must be positive")
    local = Path("data/processed/local")
    plan = local / "initial_weather_requests.csv"
    if not plan.exists():
        parser.error("Run 03_weather_requests.R first.")
    with plan.open(encoding="utf-8", newline="") as f:
        rows = [r for r in csv.DictReader(f) if r["dataset"] == DATASETS[args.dataset]]
    if args.limit:
        rows = rows[:args.limit]
    cache = local / "weather" / args.dataset
    checked = missing = 0
    client = None
    for row in rows:
        request = {"variable": VARIABLES, "date": [row["start_date"] + "/" + row["end_date"]],
                   "location": {"longitude": float(row["longitude"]), "latitude": float(row["latitude"])},
                   "data_format": "netcdf"}
        digest = hashlib.sha256(json.dumps({"dataset": row["dataset"], "request": request},
                                            sort_keys=True).encode()).hexdigest()
        final = cache / (digest + ".response")
        if final.exists():
            validate_response(final, row)
            checked += 1
            continue
        if not args.download:
            missing += 1
            continue
        if client is None:
            import cdsapi
            client = cdsapi.Client(quiet=True)
        cache.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(dir=cache, delete=False) as tmp:
            tmp_path = Path(tmp.name)
        try:
            client.retrieve(row["dataset"], request, str(tmp_path))
            validate_response(tmp_path, row)
            os.replace(tmp_path, final)
            checked += 1
        finally:
            tmp_path.unlink(missing_ok=True)
        if checked % 100 == 0:
            print(f"Validated {checked} windows", flush=True)
    print(f"{len(rows)} initial windows selected; {checked} weather responses validated; {missing} absent.")


if __name__ == "__main__":
    main()
