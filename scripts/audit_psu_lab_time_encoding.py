#!/usr/bin/env python3
"""Aggregate-only diagnostic for PSU LAB_RESULT_CM time encodings.

Tests whether lab time fields are standard clock strings, numeric minutes since midnight,
numeric seconds since midnight, or HHMM-style integers. Exports aggregate summaries only.
No identifiers, patient rows, free text values, propensity scores, or treatment effects are exported.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import duckdb
import pandas as pd

MIN_CELL = 11
CORE_LOINCS = ["2524-7", "32693-4", "19239-3", "2160-0", "38483-4", "6690-2", "777-3", "1975-2"]


def q(v: str) -> str:
    return "'" + str(v).replace("'", "''") + "'"


def qi(v: str) -> str:
    return '"' + v.replace('"', '""') + '"'


def safe(v):
    if v is None:
        return None
    try:
        n = int(v)
    except Exception:
        return None
    return n if n >= MIN_CELL else None


def find_parquet(root: Path, stems: list[str]) -> Path:
    for stem in stems:
        cands = sorted(root.glob(f"PCORnet/parquet/**/{stem}*.parquet"))
        if not cands:
            cands = sorted(root.glob(f"PCORnet/parquet/**/{stem.replace('_','')}*.parquet"))
        if cands:
            exact = [p for p in cands if p.stem.lower() == stem.lower()]
            return exact[0] if exact else max(cands, key=lambda p: p.stat().st_size)
    raise FileNotFoundError(stems)


def first(cols: set[str], names: list[str]) -> str | None:
    lut = {c.upper(): c for c in cols}
    for n in names:
        if n.upper() in lut:
            return lut[n.upper()]
    return None


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("data_root", type=Path)
    ap.add_argument("--output-dir", type=Path, required=True)
    args = ap.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    lab_path = find_parquet(args.data_root, ["lab_reduced", "lab_result_cm"])
    con = duckdb.connect()
    con.execute("PRAGMA threads=4")
    con.execute(f"CREATE VIEW lab AS SELECT * FROM read_parquet({q(str(lab_path))})")
    cols = set(con.execute("DESCRIBE lab").fetchdf()["column_name"].astype(str))

    loinc = first(cols, ["LAB_LOINC"])
    result_num = first(cols, ["RESULT_NUM", "LAB_RESULT_NUM", "RESULT_NUMERIC"])
    fields = {
        "order": (first(cols, ["LAB_ORDER_DATE"]), first(cols, ["LAB_ORDER_TIME"])),
        "specimen": (first(cols, ["SPECIMEN_DATE"]), first(cols, ["SPECIMEN_TIME"])),
        "result": (first(cols, ["RESULT_DATE"]), first(cols, ["RESULT_TIME"])),
    }
    if not loinc:
        raise RuntimeError("LAB_LOINC not found")

    code_sql = ",".join(q(x) for x in CORE_LOINCS)
    base_where = f"cast({qi(loinc)} as varchar) in ({code_sql})"
    if result_num:
        base_where += f" and try_cast({qi(result_num)} as double) is not null"

    rows = []
    interpretation_rows = []
    for clock, (date_col, time_col) in fields.items():
        if not date_col or not time_col:
            rows.append({"clock": clock, "date_field_present": bool(date_col), "time_field_present": bool(time_col)})
            continue
        t = f"trim(cast({qi(time_col)} as varchar))"
        n = con.execute(f"SELECT count(*) FROM lab WHERE {base_where}").fetchone()[0]
        nonnull = con.execute(f"SELECT count(*) FROM lab WHERE {base_where} AND {qi(time_col)} IS NOT NULL AND {t}<>''").fetchone()[0]
        colon = con.execute(f"SELECT count(*) FROM lab WHERE {base_where} AND {qi(time_col)} IS NOT NULL AND strpos({t}, ':')>0").fetchone()[0]
        numeric = con.execute(f"SELECT count(*) FROM lab WHERE {base_where} AND try_cast({t} AS DOUBLE) IS NOT NULL").fetchone()[0]
        integerish = con.execute(f"SELECT count(*) FROM lab WHERE {base_where} AND try_cast({t} AS DOUBLE) IS NOT NULL AND abs(try_cast({t} AS DOUBLE)-round(try_cast({t} AS DOUBLE)))<1e-9").fetchone()[0]
        stats = con.execute(f"""
            SELECT min(try_cast({t} AS DOUBLE)),
                   quantile_cont(try_cast({t} AS DOUBLE),0.05),
                   quantile_cont(try_cast({t} AS DOUBLE),0.50),
                   quantile_cont(try_cast({t} AS DOUBLE),0.95),
                   max(try_cast({t} AS DOUBLE)),
                   sum(case when try_cast({t} AS DOUBLE) between 0 and 23 then 1 else 0 end),
                   sum(case when try_cast({t} AS DOUBLE) between 0 and 1439 then 1 else 0 end),
                   sum(case when try_cast({t} AS DOUBLE) between 0 and 2359 then 1 else 0 end),
                   sum(case when try_cast({t} AS DOUBLE) between 0 and 86399 then 1 else 0 end)
            FROM lab WHERE {base_where} AND try_cast({t} AS DOUBLE) IS NOT NULL
        """).fetchone()
        rows.append({
            "clock": clock,
            "date_field_present": True,
            "time_field_present": True,
            "core_lab_rows": safe(n),
            "nonnull_time_rows": safe(nonnull),
            "colon_string_rows": safe(colon),
            "numeric_rows": safe(numeric),
            "integer_like_rows": safe(integerish),
            "numeric_min": stats[0], "numeric_p05": stats[1], "numeric_median": stats[2], "numeric_p95": stats[3], "numeric_max": stats[4],
            "numeric_0_23": safe(stats[5]), "numeric_0_1439": safe(stats[6]), "numeric_0_2359": safe(stats[7]), "numeric_0_86399": safe(stats[8]),
        })

        d = f"try_cast({qi(date_col)} as date)"
        txt_ts = f"try_cast(cast({qi(date_col)} as varchar)||' '||{t} as timestamp)"
        num = f"try_cast({t} as double)"
        min_ts = f"cast({d} as timestamp) + ({num} * interval 1 minute)"
        sec_ts = f"cast({d} as timestamp) + ({num} * interval 1 second)"
        hh = f"floor({num}/100)"
        mm = f"mod({num},100)"
        hhmm_valid = f"({num} between 0 and 2359 and {mm} between 0 and 59 and {hh} between 0 and 23)"
        hhmm_ts = f"cast({d} as timestamp) + ({hh} * interval 1 hour) + ({mm} * interval 1 minute)"
        for label, expr, valid in [
            ("standard_text", txt_ts, "true"),
            ("minutes_since_midnight", min_ts, f"{num} between 0 and 1439"),
            ("seconds_since_midnight", sec_ts, f"{num} between 0 and 86399"),
            ("hhmm_numeric", hhmm_ts, hhmm_valid),
        ]:
            count = con.execute(f"SELECT count(*) FROM lab WHERE {base_where} AND {d} IS NOT NULL AND {valid} AND {expr} IS NOT NULL").fetchone()[0]
            interpretation_rows.append({"clock": clock, "interpretation": label, "parseable_core_lab_rows": safe(count)})

    pd.DataFrame(rows).to_csv(args.output_dir / "time_field_distribution.csv", index=False)
    pd.DataFrame(interpretation_rows).to_csv(args.output_dir / "candidate_interpretations.csv", index=False)
    summary = {
        "privacy_mode": "aggregate_only_no_ids_no_patient_rows_no_free_text_export",
        "core_loincs": CORE_LOINCS,
        "purpose": "test whether PSU lab TIME fields are clock strings, minutes since midnight, seconds since midnight, or HHMM-style numerics",
        "guardrail": "Diagnostic only. Do not fit propensity scores or treatment effects from this task."
    }
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
