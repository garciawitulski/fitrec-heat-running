"""Count and link the two original FitRec files and independently resolve run clocks.

The provider distributes one Python dictionary literal per gzip-compressed line.
The default output contains aggregate counts. An explicit weather-preparation
option also saves indexed event times and centroids under the ignored data tree.
"""
import argparse
import ast
import csv
import functools
import gzip
import json
import re
from array import array
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


tokens = re.compile(r"""'(?:[^'\\]|\\.)*'|"(?:[^"\\]|\\.)*"|\b(?:None|True|False)\b""")
literal_words = re.compile(r"\b(?:None|True|False)\b")
json_words = {"None": "null", "True": "true", "False": "false"}


@functools.lru_cache(maxsize=4096)
def json_token(match_text):
    if match_text in json_words:
        return json_words[match_text]
    if match_text.startswith("'") and "\\" not in match_text and '"' not in match_text:
        return '"' + match_text[1:-1] + '"'
    return json.dumps(ast.literal_eval(match_text), ensure_ascii=True)


def source_rows(path):
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        for line_no, line in enumerate(stream, 1):
            try:
                if '"' not in line and "\\" not in line:
                    # With a single unescaped quote delimiter, counting quotes
                    # distinguishes a literal None/boolean from ordinary text.
                    if any(word in line for word in json_words):
                        parts = []
                        start = 0
                        quoted = False
                        for token in literal_words.finditer(line):
                            quoted ^= bool(line.count("'", start, token.start()) % 2)
                            parts.extend([line[start:token.start()], token.group() if quoted else json_words[token.group()]])
                            start = token.end()
                        parts.append(line[start:])
                        line = "".join(parts)
                    yield json.loads(line.replace("'", '"'))
                    continue
                try:
                    parsed = json.loads(line)
                except ValueError:
                    parsed = json.loads(tokens.sub(lambda m: json_token(m.group()), line))
                yield parsed
            except (SyntaxError, ValueError, TypeError):
                try:
                    yield ast.literal_eval(line)
                except (SyntaxError, ValueError) as exc:
                    raise ValueError(f"Invalid record {line_no} in {path.name}") from exc


def clock_agrees(numeric_times, iso_times):
    if len(numeric_times) < 2 or not isinstance(iso_times, list) or len(numeric_times) != len(iso_times):
        return False
    observed = np.asarray(numeric_times, dtype=float)
    if not np.isfinite(observed).all() or observed[-1] <= observed[0]:
        return False
    try:
        if all(isinstance(t, str) and t.endswith("Z") for t in iso_times):
            dates = np.asarray([t[:-1] for t in iso_times], dtype="datetime64[us]")
            if np.isnat(dates).any():
                return False
            anchor = dates.astype("int64") / 1e6
        else:
            dates = [datetime.fromisoformat(t.replace("Z", "+00:00")) for t in iso_times]
            if any(t.tzinfo is None for t in dates):
                return False
            anchor = np.asarray([t.timestamp() for t in dates])
    except (AttributeError, TypeError, ValueError, OverflowError, OSError):
        return False
    if (anchor[0] < -631152000 or anchor[-1] > datetime.now(timezone.utc).timestamp()
            or anchor[-1] <= anchor[0] or np.any(np.diff(anchor) < 0)):
        return False
    # Anchor the recorded elapsed time to the metadata UTC start. Allow at most
    # one second of disagreement at every reading, including both endpoints.
    return np.max(np.abs(anchor[0] + (observed-observed[0]) - anchor)) <= 1


def utc_seconds(text):
    try:
        stamp = datetime.fromisoformat(text.replace("Z", "+00:00"))
        return stamp.timestamp() if stamp.tzinfo is not None else np.nan
    except (ValueError, AttributeError, OSError, OverflowError):
        return np.nan


def rebuild(raw_dir, cohort=False, weather_file=None):
    hr_users = set()
    hr_sessions = {}
    running_clocks = {}
    for row in source_rows(raw_dir / "endomondoHR.json.gz"):
        session, runner = row["id"], row["userId"]
        if type(session) is not int or type(runner) is not int or session in hr_sessions:
            raise ValueError("Invalid or repeated measurement-file identifier")
        hr_sessions[session] = runner
        hr_users.add(runner)
        if row["sport"] == "run":
            times = row.get("timestamp", [])
            running_clocks[session] = array("d", times) if isinstance(times, list) else array("d")
    print(f"Scanned {len(hr_sessions):,} measurement records", flush=True)

    meta_users = set()
    meta_sessions = set()
    matched = meta_runs = matched_runs = eligible = 0
    eligible_runners = set()
    meta_clocks = {}
    for ordinal, row in enumerate(source_rows(raw_dir / "endomondoMeta.json.gz"), 1):
        session, runner = row["id"], row["userId"]
        if type(session) is not int or type(runner) is not int or session in meta_sessions:
            raise ValueError("Invalid or repeated metadata identifier")
        meta_sessions.add(session)
        meta_users.add(runner)
        is_run = row["sport"] == "run"
        meta_runs += is_run
        if session not in hr_sessions:
            continue
        if hr_sessions[session] != runner or ((session in running_clocks) != is_run):
            raise ValueError("Runner or running classification disagrees across sources")
        matched += 1
        if is_run:
            matched_runs += 1
            times = row.get("timestamp")
            numeric = running_clocks.pop(session)
            confirmed = clock_agrees(numeric, times)
            if confirmed:
                eligible += 1
                eligible_runners.add(runner)
            if cohort:
                start = utc_seconds(times[0]) if times else np.nan
                end = utc_seconds(times[-1]) if times else np.nan
                span = numeric[-1]-numeric[0] if numeric else np.nan
                meta_clocks[session] = (runner, ordinal, start, end, confirmed, span)
    if running_clocks:
        raise ValueError("Some measurement runs have no metadata")
    counts = {
        "measurement_workouts": len(hr_sessions), "measurement_users": len(hr_users),
        "measurement_running_workouts": matched_runs,
        "metadata_workouts": len(meta_sessions), "metadata_users": len(meta_users),
        "metadata_running_workouts": meta_runs, "linked_measurement_workouts": matched,
        "metadata_only_workouts": len(meta_sessions)-matched,
        "metadata_only_running_workouts": meta_runs-matched_runs,
        "running_workouts_with_confirmed_clock": eligible,
        "runners_with_confirmed_clock": len(eligible_runners),
    }
    if cohort:
        from event_selection import preweather_cohort
        print(f"Linked {matched:,} workouts; reconstructing Table S1", flush=True)
        records, users = preweather_cohort(raw_dir, meta_clocks, source_rows, weather_file)
        counts.update(table_s1_preweather_runs=records, table_s1_preweather_runners=users)
    return counts


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-dir", type=Path, default=Path("data/raw"))
    parser.add_argument("--output", type=Path, default=Path("outputs/validation/source_counts.csv"))
    parser.add_argument("--cohort", action="store_true", help="apply the Table S1 duplicate and geography rules")
    parser.add_argument("--prepare-weather-input", action="store_true",
                        help="also write locally indexed event times and centroids under the ignored data directory")
    args = parser.parse_args()
    cohort = args.cohort or args.prepare_weather_input
    weather_file = Path("data/processed/local/preweather_sessions.csv") if args.prepare_weather_input else None
    results = rebuild(args.raw_dir, cohort, weather_file)
    expected = (253020, 1104, 117902, 962190, 1364, 347558, 253020, 709170, 229656, 117781, 863)
    if cohort:
        expected += (115036, 862)
    checks = list(zip(results.items(), expected, strict=True))
    if any(actual != ref for (_, actual), ref in checks):
        raise RuntimeError("Source cohort mismatch: " + str([(k, v, ref) for (k, v), ref in checks if v != ref]))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(("quantity", "calculated", "reference"))
        writer.writerows((quantity, actual, ref) for (quantity, actual), ref in checks)
    print(f"All {len(checks)} counts match the source audit.")
