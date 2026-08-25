#!/usr/bin/env python3
"""Aggregate-only audit of how positive clinical cultures are represented in PSU PCORnet.

Patient/encounter identifiers may be used internally only to link a culture test to
organism/susceptibility follow-on rows. No identifiers, patient rows, or raw result
text are exported. All reported cells are aggregate and small cells are suppressed.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import duckdb
import pandas as pd

MIN_CELL = 11


def qident(name: str) -> str:
    return '"' + name.replace('"', '""') + '"'


def qstr(value: str) -> str:
    return "'" + value.replace("'", "''") + "'"


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
    con.execute(f"CREATE VIEW lab AS SELECT * FROM read_parquet({qstr(str(lab))})")

    cols = [r[0] for r in con.execute("DESCRIBE lab").fetchall()]
    cmap = {c.lower(): c for c in cols}
    required = [
        "patid", "encounterid", "raw_lab_name", "lab_loinc", "specimen_source",
        "specimen_date", "specimen_time", "result_date", "result_time",
        "result_qual", "raw_result",
    ]
    missing = [c for c in required if c not in cmap]
    if missing:
        raise RuntimeError(f"Missing required lab columns: {missing}")

    C = {k: qident(cmap[k]) for k in required}

    # Restrict to actual culture tests. Explicitly exclude urinalysis components whose
    # panel label merely contains 'culture if indicated'.
    culture_where = f"""
      lower(coalesce({C['raw_lab_name']},'')) LIKE '%culture%'
      AND lower(coalesce({C['raw_lab_name']},'')) NOT LIKE 'ua w culture if ind.%'
      AND lower(coalesce({C['raw_lab_name']},'')) NOT LIKE '%culture if indicated%'
    """

    follow_where = f"""
      lower(coalesce({C['raw_lab_name']},'')) LIKE '%organism%'
      OR lower(coalesce({C['raw_lab_name']},'')) LIKE '%susceptib%'
      OR lower(coalesce({C['raw_lab_name']},'')) LIKE '%sensitivity%'
      OR lower(coalesce({C['raw_lab_name']},'')) LIKE '%isolate%'
      OR lower(coalesce({C['raw_lab_name']},'')) LIKE '%identification%'
      OR lower(coalesce({C['raw_lab_name']},'')) LIKE '% mic %'
    """

    con.execute(f"""
      CREATE TEMP TABLE cultures AS
      SELECT
        {C['patid']} AS patid,
        {C['encounterid']} AS encounterid,
        coalesce({C['lab_loinc']}, '<MISSING>') AS loinc,
        coalesce({C['specimen_source']}, '<MISSING>') AS specimen_source,
        {C['specimen_date']} AS specimen_date,
        {C['specimen_time']} AS specimen_time,
        {C['result_date']} AS result_date,
        {C['result_time']} AS result_time,
        upper(trim(coalesce({C['result_qual']},''))) AS rq,
        lower(trim(coalesce({C['raw_result']},''))) AS rr,
        CASE
          WHEN lower(coalesce({C['raw_lab_name']},'')) LIKE '%blood culture%' THEN 'blood'
          WHEN lower(coalesce({C['raw_lab_name']},'')) LIKE '%urine culture%' THEN 'urine'
          WHEN lower(coalesce({C['raw_lab_name']},'')) LIKE '%respir%' THEN 'respiratory'
          WHEN lower(coalesce({C['raw_lab_name']},'')) LIKE '%wound%' THEN 'wound'
          WHEN lower(coalesce({C['raw_lab_name']},'')) LIKE '%tissue%' THEN 'tissue'
          WHEN lower(coalesce({C['raw_lab_name']},'')) LIKE '%csf%' THEN 'csf'
          WHEN lower(coalesce({C['raw_lab_name']},'')) LIKE '%fluid%' THEN 'other_fluid'
          WHEN lower(coalesce({C['raw_lab_name']},'')) LIKE '%stool%' THEN 'stool'
          WHEN lower(coalesce({C['raw_lab_name']},'')) LIKE '%fung%' THEN 'fungal'
          WHEN lower(coalesce({C['raw_lab_name']},'')) LIKE '%afb%' OR lower(coalesce({C['raw_lab_name']},'')) LIKE '%mycobacter%' THEN 'afb_mycobacterial'
          ELSE 'other_culture'
        END AS culture_group
      FROM lab
      WHERE {culture_where}
    """)

    con.execute(f"""
      CREATE TEMP TABLE followup AS
      SELECT
        {C['patid']} AS patid,
        {C['encounterid']} AS encounterid,
        {C['specimen_date']} AS specimen_date,
        {C['result_date']} AS result_date,
        upper(trim(coalesce({C['result_qual']},''))) AS rq,
        lower(trim(coalesce({C['raw_result']},''))) AS rr,
        CASE
          WHEN lower(coalesce({C['raw_lab_name']},'')) LIKE '%susceptib%' OR lower(coalesce({C['raw_lab_name']},'')) LIKE '%sensitivity%' OR lower(coalesce({C['raw_lab_name']},'')) LIKE '% mic %' THEN 'susceptibility_or_mic'
          WHEN lower(coalesce({C['raw_lab_name']},'')) LIKE '%organism%' OR lower(coalesce({C['raw_lab_name']},'')) LIKE '%isolate%' OR lower(coalesce({C['raw_lab_name']},'')) LIKE '%identification%' THEN 'organism_or_identification'
          ELSE 'other_followup'
        END AS followup_type
      FROM lab
      WHERE {follow_where}
    """)

    # Direct result semantics on the culture row itself.
    con.execute("""
      CREATE TEMP TABLE culture_classified AS
      SELECT *,
        CASE
          WHEN rq IN ('NEGATIVE','NOT DETECTED','NO GROWTH')
            OR regexp_matches(rr, '(^|\\b)(no growth|negative|not detected)(\\b|$)')
            THEN 'direct_negative'
          WHEN rq IN ('POSITIVE','PRESUMPTIVE POSITIVE','DETECTED','PRESENT')
            OR regexp_matches(rr, '(^|\\b)(positive|detected|growth of|isolated)(\\b|$)')
            THEN 'direct_positive'
          ELSE 'direct_unresolved'
        END AS direct_class
      FROM cultures
    """)

    # A follow-on organism identification or susceptibility/MIC row in the same
    # encounter within 0-7 days of culture specimen is treated only as a positivity
    # *signal* for this audit, not yet as the frozen phenotype.
    con.execute("""
      CREATE TEMP TABLE linked AS
      SELECT c.*,
        EXISTS (
          SELECT 1 FROM followup f
          WHERE f.patid = c.patid
            AND f.encounterid = c.encounterid
            AND f.result_date >= c.specimen_date
            AND f.result_date <= c.specimen_date + INTERVAL 7 DAY
            AND (f.rr <> '' OR f.rq <> '')
            AND NOT regexp_matches(f.rr, '(^|\\b)(no growth|negative|not detected)(\\b|$)')
        ) AS has_followup_signal,
        EXISTS (
          SELECT 1 FROM followup f
          WHERE f.patid = c.patid
            AND f.encounterid = c.encounterid
            AND f.result_date >= c.specimen_date
            AND f.result_date <= c.specimen_date + INTERVAL 7 DAY
            AND f.followup_type = 'susceptibility_or_mic'
        ) AS has_susceptibility_signal,
        EXISTS (
          SELECT 1 FROM followup f
          WHERE f.patid = c.patid
            AND f.encounterid = c.encounterid
            AND f.result_date >= c.specimen_date
            AND f.result_date <= c.specimen_date + INTERVAL 7 DAY
            AND f.followup_type = 'organism_or_identification'
        ) AS has_organism_signal
      FROM culture_classified c
    """)

    mechanism = con.execute("""
      SELECT
        culture_group,
        direct_class,
        has_followup_signal,
        has_organism_signal,
        has_susceptibility_signal,
        count(*)::BIGINT AS count
      FROM linked
      GROUP BY 1,2,3,4,5
      ORDER BY count DESC
    """).fetchdf()
    mechanism = suppress(mechanism)
    mechanism.to_csv(args.output_dir / "positivity_mechanism_counts.csv", index=False)

    rq = con.execute("""
      SELECT CASE WHEN rq = '' THEN '<EMPTY>' ELSE rq END AS result_qual,
             count(*)::BIGINT AS count
      FROM cultures
      GROUP BY 1 ORDER BY count DESC
      LIMIT 100
    """).fetchdf()
    rq = suppress(rq)
    rq.to_csv(args.output_dir / "culture_result_qual_counts.csv", index=False)

    follow = con.execute("""
      SELECT followup_type, count(*)::BIGINT AS count
      FROM followup GROUP BY 1 ORDER BY count DESC
    """).fetchdf()
    follow = suppress(follow)
    follow.to_csv(args.output_dir / "followup_row_type_counts.csv", index=False)

    by_loinc = con.execute("""
      SELECT loinc, culture_group,
             count(*)::BIGINT AS culture_rows,
             sum(CASE WHEN direct_class='direct_negative' THEN 1 ELSE 0 END)::BIGINT AS direct_negative,
             sum(CASE WHEN direct_class='direct_positive' THEN 1 ELSE 0 END)::BIGINT AS direct_positive,
             sum(CASE WHEN has_followup_signal THEN 1 ELSE 0 END)::BIGINT AS linked_followup_signal
      FROM linked
      GROUP BY 1,2
      HAVING count(*) >= 11
      ORDER BY culture_rows DESC
      LIMIT 200
    """).fetchdf()
    by_loinc.to_csv(args.output_dir / "culture_loinc_linkage_summary.csv", index=False)

    totals = con.execute("""
      SELECT
        count(*)::BIGINT AS culture_rows,
        sum(CASE WHEN direct_class='direct_negative' THEN 1 ELSE 0 END)::BIGINT AS direct_negative,
        sum(CASE WHEN direct_class='direct_positive' THEN 1 ELSE 0 END)::BIGINT AS direct_positive,
        sum(CASE WHEN direct_class='direct_unresolved' THEN 1 ELSE 0 END)::BIGINT AS direct_unresolved,
        sum(CASE WHEN has_followup_signal THEN 1 ELSE 0 END)::BIGINT AS any_followup_signal,
        sum(CASE WHEN has_organism_signal THEN 1 ELSE 0 END)::BIGINT AS organism_signal,
        sum(CASE WHEN has_susceptibility_signal THEN 1 ELSE 0 END)::BIGINT AS susceptibility_signal,
        sum(CASE WHEN direct_class='direct_unresolved' AND has_followup_signal THEN 1 ELSE 0 END)::BIGINT AS unresolved_rescued_by_followup
      FROM linked
    """).fetchdf().iloc[0].to_dict()
    totals = {k: int(v) for k, v in totals.items()}

    summary = {
        "privacy_mode": "aggregate_only_internal_identifier_linkage_no_ids_no_rows_no_raw_result_text_export",
        "minimum_reported_cell": MIN_CELL,
        "culture_definition": "actual RAW_LAB_NAME culture tests excluding UA reflex-panel components",
        "linkage_window": "same PATID+ENCOUNTERID, follow-on organism/susceptibility result date 0-7 days after culture specimen date",
        "totals": totals,
        "interpretation_guardrail": "Follow-on organism/susceptibility rows are positivity signals for semantic validation, not automatically the final positive-culture phenotype.",
        "next_decision": "Assess whether direct-positive plus validated organism/susceptibility linkage can define the PSU positive clinical culture available-by-day-3 rule; freeze only after aggregate plausibility review.",
    }
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
