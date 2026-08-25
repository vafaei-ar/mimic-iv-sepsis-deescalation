#!/usr/bin/env python3
"""Aggregate-only audit of PSU culture-to-followup linkage strictness.

Uses PATID/ENCOUNTERID internally only. Exports no identifiers, patient rows, or raw
result text. Compares increasingly strict linkage rules to explain prior count
inconsistency and identify a defensible positivity signal.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import duckdb
import pandas as pd

MIN_CELL = 11


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

    culture_where = r"""
      lower(coalesce(RAW_LAB_NAME,'')) LIKE '%culture%'
      AND lower(coalesce(RAW_LAB_NAME,'')) NOT LIKE 'ua w culture if ind.%'
      AND lower(coalesce(RAW_LAB_NAME,'')) NOT LIKE '%culture if indicated%'
    """

    group_expr = r"""CASE
      WHEN lower(coalesce(RAW_LAB_NAME,'')) LIKE '%blood culture%' THEN 'blood'
      WHEN lower(coalesce(RAW_LAB_NAME,'')) LIKE '%urine culture%' THEN 'urine'
      WHEN lower(coalesce(RAW_LAB_NAME,'')) LIKE '%respir%' OR lower(coalesce(RAW_LAB_NAME,'')) LIKE '%sputum%' THEN 'respiratory'
      WHEN lower(coalesce(RAW_LAB_NAME,'')) LIKE '%wound%' THEN 'wound'
      WHEN lower(coalesce(RAW_LAB_NAME,'')) LIKE '%tissue%' THEN 'tissue'
      WHEN lower(coalesce(RAW_LAB_NAME,'')) LIKE '%csf%' THEN 'csf'
      WHEN lower(coalesce(RAW_LAB_NAME,'')) LIKE '%fung%' THEN 'fungal'
      WHEN lower(coalesce(RAW_LAB_NAME,'')) LIKE '%afb%' OR lower(coalesce(RAW_LAB_NAME,'')) LIKE '%mycobact%' THEN 'afb_mycobacterial'
      WHEN lower(coalesce(RAW_LAB_NAME,'')) LIKE '%stool%' THEN 'stool'
      WHEN lower(coalesce(RAW_LAB_NAME,'')) LIKE '%fluid%' THEN 'other_fluid'
      ELSE 'other_culture' END"""

    follow_where = r"""
      lower(coalesce(RAW_LAB_NAME,'')) LIKE '%organism%'
      OR lower(coalesce(RAW_LAB_NAME,'')) LIKE '%isolate%'
      OR lower(coalesce(RAW_LAB_NAME,'')) LIKE '%identification%'
      OR lower(coalesce(RAW_LAB_NAME,'')) LIKE '%susceptib%'
      OR lower(coalesce(RAW_LAB_NAME,'')) LIKE '%sensitivity%'
      OR lower(coalesce(RAW_LAB_NAME,'')) LIKE '%minimum inhibitory%'
      OR lower(coalesce(RAW_LAB_NAME,'')) LIKE '% mic %'
    """

    con.execute(f"""
      CREATE TEMP TABLE cultures AS
      SELECT
        PATID AS patid,
        ENCOUNTERID AS encounterid,
        SPECIMEN_DATE AS specimen_date,
        coalesce(SPECIMEN_SOURCE,'<MISSING>') AS specimen_source,
        {group_expr} AS culture_group
      FROM lab
      WHERE {culture_where}
        AND PATID IS NOT NULL
        AND SPECIMEN_DATE IS NOT NULL
    """)

    con.execute(f"""
      CREATE TEMP TABLE followup AS
      SELECT
        PATID AS patid,
        ENCOUNTERID AS encounterid,
        SPECIMEN_DATE AS specimen_date,
        RESULT_DATE AS result_date,
        coalesce(SPECIMEN_SOURCE,'<MISSING>') AS specimen_source,
        CASE
          WHEN lower(coalesce(RAW_LAB_NAME,'')) LIKE '%susceptib%'
            OR lower(coalesce(RAW_LAB_NAME,'')) LIKE '%sensitivity%'
            OR lower(coalesce(RAW_LAB_NAME,'')) LIKE '%minimum inhibitory%'
            OR lower(coalesce(RAW_LAB_NAME,'')) LIKE '% mic %'
          THEN 'susceptibility_or_mic'
          ELSE 'organism_or_identification'
        END AS followup_type,
        CASE WHEN trim(coalesce(cast(RAW_RESULT AS VARCHAR),'')) <> ''
                  OR trim(coalesce(cast(RESULT_QUAL AS VARCHAR),'')) <> ''
             THEN 1 ELSE 0 END AS has_nonblank_result
      FROM lab
      WHERE {follow_where}
        AND PATID IS NOT NULL
        AND RESULT_DATE IS NOT NULL
    """)

    rules = [
        ("patient_only_0_7d", "f.patid=c.patid AND f.result_date BETWEEN c.specimen_date AND c.specimen_date + INTERVAL 7 DAY"),
        ("patient_encounter_exact_0_7d", "f.patid=c.patid AND f.encounterid=c.encounterid AND f.result_date BETWEEN c.specimen_date AND c.specimen_date + INTERVAL 7 DAY"),
        ("patient_encounter_exact_0_7d_nonblank", "f.patid=c.patid AND f.encounterid=c.encounterid AND f.result_date BETWEEN c.specimen_date AND c.specimen_date + INTERVAL 7 DAY AND f.has_nonblank_result=1"),
        ("patient_encounter_same_day", "f.patid=c.patid AND f.encounterid=c.encounterid AND f.result_date=c.specimen_date"),
        ("patient_encounter_0_1d", "f.patid=c.patid AND f.encounterid=c.encounterid AND f.result_date BETWEEN c.specimen_date AND c.specimen_date + INTERVAL 1 DAY"),
        ("patient_encounter_0_2d", "f.patid=c.patid AND f.encounterid=c.encounterid AND f.result_date BETWEEN c.specimen_date AND c.specimen_date + INTERVAL 2 DAY"),
        ("patient_encounter_0_7d_specimen_match", "f.patid=c.patid AND f.encounterid=c.encounterid AND f.result_date BETWEEN c.specimen_date AND c.specimen_date + INTERVAL 7 DAY AND f.specimen_source=c.specimen_source AND c.specimen_source <> '<MISSING>'"),
        ("patient_encounter_0_2d_specimen_match", "f.patid=c.patid AND f.encounterid=c.encounterid AND f.result_date BETWEEN c.specimen_date AND c.specimen_date + INTERVAL 2 DAY AND f.specimen_source=c.specimen_source AND c.specimen_source <> '<MISSING>'"),
    ]

    rows = []
    for rule_name, predicate in rules:
        q = f"""
          SELECT c.culture_group,
                 f.followup_type,
                 count(DISTINCT (c.patid, c.encounterid, c.specimen_date, c.culture_group, c.specimen_source))::BIGINT AS linked_culture_units
          FROM cultures c
          JOIN followup f ON {predicate}
          GROUP BY 1,2
        """
        d = con.execute(q).fetchdf()
        d.insert(0, "linkage_rule", rule_name)
        rows.append(d)
    detail = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    detail = detail.rename(columns={"linked_culture_units": "count"})
    detail = suppress(detail)
    detail.to_csv(args.output_dir / "linkage_rule_counts.csv", index=False)

    total_rows = []
    for rule_name, predicate in rules:
        q = f"""
          SELECT
            count(DISTINCT (c.patid, c.encounterid, c.specimen_date, c.culture_group, c.specimen_source))::BIGINT AS linked_units,
            count(DISTINCT CASE WHEN f.followup_type='organism_or_identification'
                     THEN (c.patid, c.encounterid, c.specimen_date, c.culture_group, c.specimen_source) END)::BIGINT AS organism_units,
            count(DISTINCT CASE WHEN f.followup_type='susceptibility_or_mic'
                     THEN (c.patid, c.encounterid, c.specimen_date, c.culture_group, c.specimen_source) END)::BIGINT AS susceptibility_units
          FROM cultures c
          JOIN followup f ON {predicate}
        """
        r = con.execute(q).fetchdf().iloc[0].to_dict()
        total_rows.append({"linkage_rule": rule_name, **{k: int(v) for k, v in r.items()}})
    totals = pd.DataFrame(total_rows)
    totals.to_csv(args.output_dir / "linkage_rule_totals.csv", index=False)

    # Quantify missing ENCOUNTERID and specimen-source availability because permissive
    # null encounter matching was the likely source of inflation in the prior audit.
    completeness = con.execute("""
      SELECT 'culture' AS row_type,
             count(*)::BIGINT AS rows,
             sum(CASE WHEN encounterid IS NULL THEN 1 ELSE 0 END)::BIGINT AS encounterid_missing,
             sum(CASE WHEN specimen_source='<MISSING>' THEN 1 ELSE 0 END)::BIGINT AS specimen_source_missing
      FROM cultures
      UNION ALL
      SELECT 'followup', count(*)::BIGINT,
             sum(CASE WHEN encounterid IS NULL THEN 1 ELSE 0 END)::BIGINT,
             sum(CASE WHEN specimen_source='<MISSING>' THEN 1 ELSE 0 END)::BIGINT
      FROM followup
    """).fetchdf()
    completeness.to_csv(args.output_dir / "linkage_field_completeness.csv", index=False)

    summary = {
        "privacy_mode": "aggregate_only_internal_identifier_linkage_no_ids_no_rows_no_raw_result_text_export",
        "minimum_reported_cell": MIN_CELL,
        "purpose": "Explain discrepant prior follow-up linkage counts and identify a defensible culture-to-organism/susceptibility linkage rule.",
        "culture_unit": "distinct PATID+ENCOUNTERID+SPECIMEN_DATE+culture_group+specimen_source, used internally only",
        "rules_compared": [name for name, _ in rules],
        "guardrail": "This is a linkage-diagnostic audit only. No rule is frozen as the positive-culture phenotype until aggregate plausibility is reviewed.",
        "expected_explanation": "The earlier inflated counts may reflect permissive matching when ENCOUNTERID was null and/or broad 7-day same-patient windows; exact-encounter and specimen-matched rules test this directly.",
    }
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
