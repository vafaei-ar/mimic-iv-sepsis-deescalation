#!/usr/bin/env python3
"""Aggregate-only audit to freeze PSU microbiology semantics for external replication.

Reads the local PSU PCORnet lab_reduced parquet and exports only aggregate test-name/
LOINC/specimen counts, aggregate result-semantic counts, and timing coverage. It never
exports patient identifiers, patient rows, or raw result text.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import duckdb
import pandas as pd

MIN_CELL = 11


def suppress(df: pd.DataFrame, count_col: str = "count") -> pd.DataFrame:
    out = df.copy()
    if count_col in out:
        mask = out[count_col].fillna(0) < MIN_CELL
        out.loc[mask, count_col] = pd.NA
        out["suppressed"] = mask
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("data_root", type=Path)
    ap.add_argument("--output-dir", type=Path, required=True)
    args = ap.parse_args()
    outdir = args.output_dir
    outdir.mkdir(parents=True, exist_ok=True)

    lab = args.data_root / "PCORnet" / "parquet" / "lab" / "lab_reduced.parquet"
    if not lab.exists():
        raise FileNotFoundError(lab)

    con = duckdb.connect()
    con.execute("PRAGMA threads=4")
    # DuckDB does not allow prepared parameters inside CREATE VIEW table functions.
    # Quote the already validated local path as a SQL string literal instead.
    lab_sql = str(lab).replace("'", "''")
    con.execute(f"CREATE VIEW lab AS SELECT * FROM read_parquet('{lab_sql}')")

    # Candidate culture-like tests. Keep standardized/name fields only, never result free text.
    candidate_where = r"""
      lower(coalesce(RAW_LAB_NAME,'')) LIKE '%culture%'
      OR lower(coalesce(RAW_LAB_NAME,'')) LIKE '%organism%'
      OR lower(coalesce(RAW_LAB_NAME,'')) LIKE '%susceptib%'
      OR lower(coalesce(RAW_LAB_NAME,'')) LIKE '%gram stain%'
      OR lower(coalesce(RAW_LAB_NAME,'')) LIKE '%afb%'
      OR lower(coalesce(RAW_LAB_NAME,'')) LIKE '%fungal%'
    """

    tests = con.execute(f"""
      SELECT
        coalesce(LAB_LOINC,'<MISSING>') AS loinc,
        coalesce(RAW_LAB_NAME,'<MISSING>') AS raw_lab_name,
        coalesce(SPECIMEN_SOURCE,'<MISSING>') AS specimen_source,
        count(*)::BIGINT AS count
      FROM lab
      WHERE {candidate_where}
      GROUP BY 1,2,3
      HAVING count(*) >= {MIN_CELL}
      ORDER BY count DESC
      LIMIT 500
    """).fetchdf()
    tests["suppressed"] = False
    tests.to_csv(outdir / "candidate_culture_tests.csv", index=False)

    # Result semantics are classified internally; raw result strings are never exported.
    semantics = con.execute(f"""
      WITH c AS (
        SELECT
          upper(trim(coalesce(RESULT_QUAL,''))) AS rq,
          lower(coalesce(RAW_RESULT,'')) AS rr
        FROM lab
        WHERE {candidate_where}
      ), cls AS (
        SELECT CASE
          WHEN rq IN ('NEGATIVE','NOT DETECTED','NO GROWTH')
               OR regexp_matches(rr, '(^|\\b)(no growth|negative|not detected)(\\b|$)') THEN 'negative_or_no_growth'
          WHEN rq IN ('POSITIVE','PRESUMPTIVE POSITIVE','DETECTED','PRESENT')
               OR regexp_matches(rr, '(^|\\b)(positive|detected|growth of|isolated)(\\b|$)') THEN 'positive_or_growth_signal'
          WHEN rq = '' AND rr = '' THEN 'missing_result_semantics'
          ELSE 'other_or_unresolved'
        END AS semantic_class
        FROM c
      )
      SELECT semantic_class, count(*)::BIGINT AS count
      FROM cls GROUP BY 1 ORDER BY count DESC
    """).fetchdf()
    semantics = suppress(semantics)
    semantics.to_csv(outdir / "result_semantics.csv", index=False)

    timing = con.execute(f"""
      SELECT
        count(*)::BIGINT AS candidate_rows,
        sum(CASE WHEN SPECIMEN_DATE IS NOT NULL THEN 1 ELSE 0 END)::BIGINT AS specimen_date_present,
        sum(CASE WHEN SPECIMEN_TIME IS NOT NULL THEN 1 ELSE 0 END)::BIGINT AS specimen_time_present,
        sum(CASE WHEN RESULT_DATE IS NOT NULL THEN 1 ELSE 0 END)::BIGINT AS result_date_present,
        sum(CASE WHEN RESULT_TIME IS NOT NULL THEN 1 ELSE 0 END)::BIGINT AS result_time_present,
        sum(CASE WHEN SPECIMEN_DATE IS NOT NULL AND RESULT_DATE IS NOT NULL THEN 1 ELSE 0 END)::BIGINT AS both_dates_present
      FROM lab WHERE {candidate_where}
    """).fetchdf()
    timing.to_csv(outdir / "timing_coverage.csv", index=False)

    specimen = con.execute(f"""
      SELECT coalesce(SPECIMEN_SOURCE,'<MISSING>') AS specimen_source,
             count(*)::BIGINT AS count
      FROM lab WHERE {candidate_where}
      GROUP BY 1 ORDER BY count DESC LIMIT 100
    """).fetchdf()
    specimen = suppress(specimen)
    specimen.to_csv(outdir / "specimen_source_counts.csv", index=False)

    summary = {
        "privacy_mode": "aggregate_only_no_ids_no_rows_no_raw_result_text_export",
        "minimum_reported_cell": MIN_CELL,
        "candidate_definition": "culture/organism/susceptibility/gram stain/AFB/fungal terms in RAW_LAB_NAME",
        "candidate_test_groups_exported": int(len(tests)),
        "candidate_rows": int(timing.loc[0, "candidate_rows"]),
        "timing": {k: int(timing.loc[0, k]) for k in timing.columns},
        "decision_rule_status": "not_yet_frozen",
        "next_decision": "Use aggregate candidate test groups and result-semantic classes to define a validated clinical-culture subset and a positive-result-by-day-3 rule before cohort construction.",
    }
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
