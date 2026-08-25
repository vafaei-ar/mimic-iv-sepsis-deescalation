#!/usr/bin/env python3
"""Aggregate-only audit of how PSU culture results encode positivity/negativity.

Raw result text and identifiers are inspected only locally. Exports contain only
predefined semantic classes and aggregate counts with small-cell suppression.
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
    if count_col in out.columns:
        mask = out[count_col].fillna(0) < MIN_CELL
        out.loc[mask, count_col] = pd.NA
        out["suppressed"] = mask
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("data_root", type=Path)
    ap.add_argument("--output-dir", type=Path, required=True)
    args = ap.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    lab = args.data_root / "PCORnet" / "parquet" / "lab" / "lab_reduced.parquet"
    if not lab.exists():
        raise FileNotFoundError(lab)

    con = duckdb.connect()
    con.execute("PRAGMA threads=4")
    p = str(lab).replace("'", "''")
    con.execute(f"CREATE VIEW lab AS SELECT * FROM read_parquet('{p}')")

    # Restrict to actual culture tests; exclude urinalysis reflex-panel components.
    culture_where = r"""
      (
        lower(coalesce(RAW_LAB_NAME,'')) LIKE '%culture%'
        OR lower(coalesce(RAW_LAB_NAME,'')) LIKE '%cx%'
      )
      AND lower(coalesce(RAW_LAB_NAME,'')) NOT LIKE 'ua w culture if ind.%'
      AND lower(coalesce(RAW_LAB_NAME,'')) NOT LIKE '%culture if ind.%|%'
    """

    group_expr = r"""CASE
      WHEN lower(coalesce(RAW_LAB_NAME,'')) LIKE '%blood%' THEN 'blood'
      WHEN lower(coalesce(RAW_LAB_NAME,'')) LIKE '%urine%' THEN 'urine'
      WHEN lower(coalesce(RAW_LAB_NAME,'')) LIKE '%resp%' OR lower(coalesce(RAW_LAB_NAME,'')) LIKE '%sputum%' THEN 'respiratory'
      WHEN lower(coalesce(RAW_LAB_NAME,'')) LIKE '%wound%' THEN 'wound'
      WHEN lower(coalesce(RAW_LAB_NAME,'')) LIKE '%tissue%' THEN 'tissue'
      WHEN lower(coalesce(RAW_LAB_NAME,'')) LIKE '%csf%' THEN 'csf'
      WHEN lower(coalesce(RAW_LAB_NAME,'')) LIKE '%fung%' THEN 'fungal'
      WHEN lower(coalesce(RAW_LAB_NAME,'')) LIKE '%afb%' OR lower(coalesce(RAW_LAB_NAME,'')) LIKE '%mycobact%' THEN 'afb_mycobacterial'
      WHEN lower(coalesce(RAW_LAB_NAME,'')) LIKE '%stool%' THEN 'stool'
      WHEN lower(coalesce(RAW_LAB_NAME,'')) LIKE '%fluid%' THEN 'other_fluid'
      ELSE 'other_culture' END"""

    # Classify raw result locally. Never export raw result strings.
    semantic_expr = r"""CASE
      WHEN RAW_RESULT IS NULL OR trim(cast(RAW_RESULT AS VARCHAR)) = '' THEN 'missing_or_blank'
      WHEN regexp_matches(lower(cast(RAW_RESULT AS VARCHAR)), 'no growth|no organisms? isolated|negative|not detected|none isolated') THEN 'explicit_negative_no_growth'
      WHEN regexp_matches(lower(cast(RAW_RESULT AS VARCHAR)), 'mixed (urogenital|skin|respiratory|oral|flora)|normal (respiratory|oral|skin|urogenital) flora|commensal flora') THEN 'mixed_or_normal_flora'
      WHEN regexp_matches(lower(cast(RAW_RESULT AS VARCHAR)), 'contaminant|contamination') THEN 'contamination_language'
      WHEN regexp_matches(lower(cast(RAW_RESULT AS VARCHAR)), 'positive|growth|isolated|isolate|organism') THEN 'explicit_growth_or_organism_language'
      WHEN regexp_matches(lower(cast(RAW_RESULT AS VARCHAR)), 'susceptib|resistan|\bmic\b') THEN 'susceptibility_language'
      WHEN regexp_matches(lower(cast(RAW_RESULT AS VARCHAR)), '^[0-9.+<> =-]+$') THEN 'numeric_or_symbolic'
      ELSE 'other_text_pattern' END"""

    by_group = con.execute(f"""
      SELECT {group_expr} AS culture_group,
             {semantic_expr} AS result_pattern,
             count(*)::BIGINT AS count
      FROM lab
      WHERE {culture_where}
      GROUP BY 1,2
      ORDER BY 1, count DESC
    """).fetchdf()
    suppress(by_group).to_csv(args.output_dir / "result_pattern_by_culture_group.csv", index=False)

    coverage = con.execute(f"""
      SELECT {group_expr} AS culture_group,
             count(*)::BIGINT AS culture_rows,
             sum(CASE WHEN RAW_RESULT IS NOT NULL AND trim(cast(RAW_RESULT AS VARCHAR)) <> '' THEN 1 ELSE 0 END)::BIGINT AS raw_result_present,
             sum(CASE WHEN RESULT_QUAL IS NOT NULL THEN 1 ELSE 0 END)::BIGINT AS result_qual_present,
             sum(CASE WHEN RESULT_DATE IS NOT NULL AND RESULT_TIME IS NOT NULL THEN 1 ELSE 0 END)::BIGINT AS result_datetime_present
      FROM lab
      WHERE {culture_where}
      GROUP BY 1 ORDER BY culture_rows DESC
    """).fetchdf()
    coverage.to_csv(args.output_dir / "culture_result_coverage.csv", index=False)

    # Aggregate standardized metadata that may distinguish result components.
    meta = con.execute(f"""
      SELECT {group_expr} AS culture_group,
             coalesce(LAB_LOINC,'<MISSING>') AS loinc,
             coalesce(LAB_RESULT_SOURCE,'<MISSING>') AS lab_result_source,
             coalesce(LAB_LOINC_SOURCE,'<MISSING>') AS loinc_source,
             coalesce(RESULT_LOC,'<MISSING>') AS result_loc,
             count(*)::BIGINT AS count
      FROM lab
      WHERE {culture_where}
      GROUP BY 1,2,3,4,5
      HAVING count(*) >= {MIN_CELL}
      ORDER BY count DESC
      LIMIT 500
    """).fetchdf()
    meta["suppressed"] = False
    meta.to_csv(args.output_dir / "culture_metadata_patterns.csv", index=False)

    # Safely quantify whether plausible positive text patterns agree with follow-on microbiology rows
    # within the same encounter and a 7-day window. No identifiers are exported.
    con.execute(f"""
      CREATE TEMP VIEW cultures AS
      SELECT PATID, ENCOUNTERID, SPECIMEN_DATE, RESULT_DATE,
             {group_expr} AS culture_group,
             {semantic_expr} AS result_pattern
      FROM lab WHERE {culture_where}
    """)
    follow_where = r"""
      lower(coalesce(RAW_LAB_NAME,'')) LIKE '%organism%'
      OR lower(coalesce(RAW_LAB_NAME,'')) LIKE '%identification%'
      OR lower(coalesce(RAW_LAB_NAME,'')) LIKE '%susceptib%'
      OR lower(coalesce(RAW_LAB_NAME,'')) LIKE '%minimum inhibitory%'
      OR lower(coalesce(RAW_LAB_NAME,'')) LIKE '% mic%'
    """
    con.execute(f"""
      CREATE TEMP VIEW followup AS
      SELECT PATID, ENCOUNTERID, RESULT_DATE,
             CASE WHEN lower(coalesce(RAW_LAB_NAME,'')) LIKE '%susceptib%'
                       OR lower(coalesce(RAW_LAB_NAME,'')) LIKE '%minimum inhibitory%'
                       OR lower(coalesce(RAW_LAB_NAME,'')) LIKE '% mic%'
                  THEN 'susceptibility_or_mic' ELSE 'organism_or_identification' END AS followup_type
      FROM lab WHERE {follow_where}
    """)

    concord = con.execute("""
      SELECT c.culture_group, c.result_pattern,
             count(*)::BIGINT AS culture_rows,
             sum(CASE WHEN EXISTS (
               SELECT 1 FROM followup f
               WHERE f.PATID=c.PATID
                 AND (f.ENCOUNTERID=c.ENCOUNTERID OR c.ENCOUNTERID IS NULL OR f.ENCOUNTERID IS NULL)
                 AND f.RESULT_DATE BETWEEN c.SPECIMEN_DATE AND c.SPECIMEN_DATE + INTERVAL 7 DAY
             ) THEN 1 ELSE 0 END)::BIGINT AS rows_with_followup_signal
      FROM cultures c
      GROUP BY 1,2
      ORDER BY 1, culture_rows DESC
    """).fetchdf()
    concord.to_csv(args.output_dir / "result_pattern_followup_concordance.csv", index=False)

    total = int(coverage["culture_rows"].sum())
    present = int(coverage["raw_result_present"].sum())
    summary = {
        "privacy_mode": "aggregate_only_internal_raw_result_and_identifier_inspection_no_raw_text_no_ids_no_rows_export",
        "minimum_reported_cell": MIN_CELL,
        "culture_rows": total,
        "raw_result_present": present,
        "raw_result_present_fraction": (present / total) if total else None,
        "purpose": "Determine whether direct RAW_RESULT semantic patterns can recover culture positivity/negativity beyond organism/susceptibility linkage.",
        "decision_rule_status": "diagnostic_not_yet_frozen",
        "next_decision": "Freeze PSU culture positivity only if direct result-pattern and follow-up linkage evidence supports a clinically plausible rule; otherwise treat microbiology positivity as unresolved and modify external eligibility accordingly.",
    }
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
