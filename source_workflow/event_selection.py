"""Reconstruct the pre-weather cohort from the original sensor records."""
import csv
from collections import defaultdict
from datetime import datetime, timezone
from hashlib import sha256
from itertools import combinations
from pathlib import Path
from tempfile import NamedTemporaryFile

import numpy as np
from pyproj import Geod


FIELDS = ("timestamp", "heart_rate", "latitude", "longitude", "altitude", "speed")
GEO = Geod(ellps="WGS84")
FIRST_UTC = datetime(1990, 1, 1, tzinfo=timezone.utc).timestamp()
LAST_UTC = datetime(2027, 1, 1, tzinfo=timezone.utc).timestamp()


def sensor_features(row, simultaneous, metadata_span):
    # Preserve original encounter order and use every original numeric channel
    # to decide whether two traces or two rows of readings are identical.
    arrays = []
    parse_ok = True
    for field in FIELDS:
        try:
            value = row.get(field)
            col = np.asarray([] if value is None else value, dtype="<f8")
            if col.ndim != 1:
                raise ValueError()
        except (ValueError, TypeError):
            col = np.array([], dtype="<f8")
            parse_ok = False
        col[col == 0] = 0.0
        col[np.isnan(col)] = np.nan
        arrays.append(col)
    t, hr, lat, lon = arrays[:4]
    n = len(t)
    aligned = n and all(len(col) == n for col in arrays[1:4])
    signature = None
    coverage = 0.0
    coherent = False
    if simultaneous:
        if parse_ok:
            digest = sha256()
            for col in arrays:
                digest.update(len(col).to_bytes(8, "little"))
                digest.update(col.tobytes())
            signature = digest.digest()
        if aligned:
            good_hr = np.isfinite(hr) & (hr >= 40) & (hr <= 220)
            good_gps = (np.isfinite(lat) & np.isfinite(lon) & (np.abs(lat) <= 90)
                        & (np.abs(lon) <= 180) & ((lat != 0) | (lon != 0)))
            coverage = min(good_hr.mean(), good_gps.mean())
        coherent = (n >= 2 and np.isfinite(t).all() and np.all(np.diff(t) >= 0)
                    and np.isfinite(metadata_span) and abs(metadata_span - (t[-1] - t[0])) <= 1)
    if aligned and n >= 2 and np.isfinite(t).all() and len(np.unique(np.floor(t))) < n:
        padded = [col if len(col) == n else np.full(n, np.nan) for col in arrays]
        matrix = np.ascontiguousarray(np.column_stack(padded))
        # Drop duplicate numeric rows before combining runs of equal seconds.
        keys = matrix.view(np.dtype((np.void, matrix.itemsize * matrix.shape[1]))).ravel()
        unique = matrix[np.sort(np.unique(keys, return_index=True)[1])]
        edges = np.r_[0, np.flatnonzero(np.diff(np.floor(unique[:, 0])) != 0) + 1, len(unique)]
        gps = unique[edges[:-1], 2:4].copy()
        for i in np.flatnonzero(np.diff(edges) > 1):
            group = unique[edges[i]:edges[i+1], 2:4]
            for j in range(2):
                values = group[np.isfinite(group[:, j]), j]
                gps[i, j] = np.median(values) if len(values) else np.nan
        lat, lon = gps.T
    elif len(lat) != n or len(lon) != n:
        lat = lon = np.empty(0)
    ok = np.isfinite(lat) & np.isfinite(lon) & (np.abs(lat) <= 90) & (np.abs(lon) <= 180) & ((lat != 0) | (lon != 0))
    lat, lon = lat[ok], lon[ok]
    point = bounds = None
    cells = set()
    if len(lat):
        point = (float(np.median(lat)), float(np.median(lon)))
        if simultaneous:
            bounds = (lon[0], lat[0], lon[-1], lat[-1])
            phi, lam = np.radians(lat), np.radians(lon)
            xyz = np.column_stack((np.cos(phi) * np.cos(lam), np.cos(phi) * np.sin(lam), np.sin(phi)))
            cells = set(map(tuple, np.floor(xyz * (6371008.8 / 100)).astype("int32")))
    return {"signature": signature, "point": point, "ends": bounds, "cells": cells,
            "coherent": coherent, "coverage": coverage}


def keep_event(members):
    if len(members) == 1:
        return members[0]
    probable = distinguishable = True
    for a, b in combinations(members, 2):
        if a["signature"] is not None and a["signature"] == b["signature"]:
            distinguishable = False
            continue
        if a["ends"] is None or b["ends"] is None:
            return None
        distances = [GEO.inv(*a["ends"][i:i+2], *b["ends"][i:i+2])[2] for i in (0, 2)]
        cells = a["cells"] | b["cells"]
        jaccard = len(a["cells"] & b["cells"]) / len(cells) if cells else 0.
        spans = [a["meta_span"], b["meta_span"]]
        relative = abs(spans[0] - spans[1]) / max(spans) if min(spans) > 0 else np.inf
        probable &= max(distances) <= 50 and jaccard >= .9 and relative <= .02
        distinguishable &= min(distances) >= 1000
    if not probable and not distinguishable:
        return None
    return min(members, key=lambda x: (not x["eligible"], not x["coherent"], -x["coverage"], x["order"]))


def preweather_cohort(raw_dir, metadata, read_rows, output=None):
    groups = defaultdict(list)
    for session, info in metadata.items():
        runner, ordinal, start, end, eligible, run_span = info
        if np.isfinite(start):
            groups[(runner, start)].append(session)
    simultaneous = {s for members in groups.values() if len(members) > 1 for s in members}
    for row in read_rows(raw_dir / "endomondoHR.json.gz"):
        session = row["id"]
        if session not in metadata:
            continue
        runner, order, start, end, eligible, span = metadata[session]
        features = sensor_features(row, session in simultaneous, end - start)
        features.update(runner=runner, order=order, start=start, end=start + span,
                        meta_span=end - start, eligible=eligible)
        metadata[session] = features
    cohort = []
    for sessions in groups.values():
        event = keep_event([metadata[s] for s in sessions])
        if (event is not None and event["eligible"] and event["point"] is not None and
                event["start"] >= FIRST_UTC and event["end"] <= LAST_UTC and
                0 < event["end"] - event["start"] <= 86400):
            cohort.append(event)
    runners = sorted({e["runner"] for e in cohort})
    if (len(cohort), len(runners)) != (115036, 862):
        raise ValueError(f"Table S1 cohort differs: {len(cohort)} runs, {len(runners)} runners")
    if output:
        # Write to this clone's ignored data directory only. The weather input
        # contains positions and event times and must remain on this machine.
        output = Path(output)
        output.parent.mkdir(parents=True, exist_ok=True)
        runner_code = {runner: i for i, runner in enumerate(runners, 1)}
        with NamedTemporaryFile(dir=output.parent, delete=False) as pending:
            pending_path = Path(pending.name)
        try:
            with pending_path.open("w", encoding="utf-8", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(("session_index", "runner_index", "start_utc_seconds", "end_utc_seconds",
                                 "median_latitude", "median_longitude"))
                for index, event in enumerate(sorted(cohort, key=lambda e: e["order"]), 1):
                    writer.writerow((index, runner_code[event["runner"]], event["start"], event["end"], *event["point"]))
            pending_path.replace(output)
        finally:
            pending_path.unlink(missing_ok=True)
    return len(cohort), len(runners)
