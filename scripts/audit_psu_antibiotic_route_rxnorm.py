#!/usr/bin/env python3
"""Aggregate-only audit of PSU antibiotic route coding and legacy RxNorm coverage."""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import duckdb
import pandas as pd

MIN_CELL = 11
BROAD_PATTERN = (
    "vancomycin|piperacillin|tazobactam|zosyn|cefepime|ceftazidime|meropenem|"
    "imipenem|aztreonam|linezolid|daptomycin|ceftolozane|avibactam"
)


def q(value: str) -> str:
    return "'" + value.replace("'", "''") + "'"


def suppress(df: pd.DataFrame, col: str = "count") -> pd.DataFrame:
    out = df.copy()
    if col in out.columns:
        mask = out[col].fillna(0) < MIN_CELL
        out.loc[mask, col] = pd.NA
        out["suppressed"] = mask
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("data_root", type=Path)
    ap.add_argument("--output-dir", type=Path, required=True)
    args = ap.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    root = args.data_root / "PCORnet"
    prescribing = root / "parquet" / "prescribing.parquet"
    med_admin = root / "parquet" / "med_admin.parquet"
    old_code = root / "code" / "config" / "codes_antibiotics.py"
    if not prescribing.exists() or not med_admin.exists():
        raise FileNotFoundError("Required prescribing or med_admin parquet missing")

    old_text = old_code.read_text(encoding="utf-8", errors="ignore") if old_code.exists() else ""
    rxnorm_literals = sorted(set(re.findall(r"(?<![A-Za-z0-9_])(\d{4,8})(?![A-Za-z0-9_])", old_text)))
    rxnorm_sql = ",".join(q(x) for x in rxnorm_literals) or "''"

    con = duckdb.connect()
    con.execute("PRAGMA threads=4")
    con.execute(f"CREATE VIEW p AS SELECT * FROM read_parquet({q(str(prescribing))})")
    con.execute(f"CREATE VIEW m AS SELECT * FROM read_parquet({q(str(med_admin))})")

    route_queries = [
        ("PRESCRIBING", "p", "raw_rx_med_name", "rx_route"),
        ("MED_ADMIN", "m", "raw_medadmin_med_name", "medadmin_route"),
    ]
    route_frames = []
    for source, view, text_col, route_col in route_queries:
        df = con.execute(f"""
            SELECT {q(source)} AS source,
                   coalesce(cast({route_col} AS VARCHAR), '<MISSING>') AS route,
                   count(*)::BIGINT AS count
            FROM {view}
            WHERE regexp_matches(lower(coalesce(cast({text_col} AS VARCHAR),'')), {q(BROAD_PATTERN)})
            GROUP BY 1,2
            ORDER BY count DESC
            LIMIT 100
        """).fetchdf()
        route_frames.append(suppress(df))
    pd.concat(route_frames, ignore_index=True).to_csv(args.output_dir / "broad_route_value_counts.csv", index=False)

    mapping_rows = []
    specs = [
        ("PRESCRIBING", "p", "raw_rx_med_name", "rxnorm_cui", "rx_route", "encounterid"),
        ("MED_ADMIN", "m", "raw_medadmin_med_name", "medadmin_code", "medadmin_route", "encounterid"),
    ]
    for source, view, text_col, code_col, route_col, enc_col in specs:
        row = con.execute(f"""
            SELECT
              {q(source)} AS source,
              count(*)::BIGINT AS total_rows,
              sum(regexp_matches(lower(coalesce(cast({text_col} AS VARCHAR),'')), {q(BROAD_PATTERN)}))::BIGINT AS broad_text_rows,
              sum(cast({code_col} AS VARCHAR) IN ({rxnorm_sql}))::BIGINT AS legacy_rxnorm_rows,
              sum(regexp_matches(lower(coalesce(cast({text_col} AS VARCHAR),'')), {q(BROAD_PATTERN)})
                  AND cast({code_col} AS VARCHAR) IN ({rxnorm_sql}))::BIGINT AS broad_text_and_legacy_rxnorm_rows,
              count(DISTINCT CASE WHEN regexp_matches(lower(coalesce(cast({text_col} AS VARCHAR),'')), {q(BROAD_PATTERN)}) THEN {enc_col} END)::BIGINT AS broad_text_encounters,
              count(DISTINCT CASE WHEN cast({code_col} AS VARCHAR) IN ({rxnorm_sql}) THEN {enc_col} END)::BIGINT AS legacy_rxnorm_encounters
            FROM {view}
        """).fetchdf().iloc[0].to_dict()
        mapping_rows.append(row)
    pd.DataFrame(mapping_rows).to_csv(args.output_dir / "rxnorm_mapping_summary.csv", index=False)

    code_counts = []
    for source, view, code_col in [("PRESCRIBING", "p", "rxnorm_cui"), ("MED_ADMIN", "m", "medadmin_code")]:
        df = con.execute(f"""
            SELECT {q(source)} AS source, cast({code_col} AS VARCHAR) AS code, count(*)::BIGINT AS count
            FROM {view}
            WHERE cast({code_col} AS VARCHAR) IN ({rxnorm_sql})
            GROUP BY 1,2
            ORDER BY count DESC
        """).fetchdf()
        code_counts.append(suppress(df))
    pd.concat(code_counts, ignore_index=True).to_csv(args.output_dir / "legacy_rxnorm_code_counts.csv", index=False)

    summary = {
        "privacy_mode": "aggregate_only_no_ids_no_patient_rows_no_free_text_export",
        "minimum_reported_cell": MIN_CELL,
        "frozen_mimic_broad_pattern": BROAD_PATTERN,
        "legacy_rxnorm_literal_count": len(rxnorm_literals),
        "legacy_code_present": old_code.exists(),
        "purpose": "Determine whether PSU IV route fields can be used directly and whether the prior local RxNorm mapping provides a stronger route-independent broad-spectrum phenotype than literal IV-route coding.",
        "guardrail": "Diagnostic only; no antibiotic source or phenotype is frozen from this audit alone.",
    }
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
