#!/usr/bin/env python3
"""Aggregate-only audit of exact specimen-timestamp linkage for PSU cultures.

Tests whether culture parent rows and nearby organism/susceptibility/Gram-stain
components can be linked using exact specimen date/time and/or lab order date.
Identifiers and raw results remain local. Exports contain aggregate counts only.
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

    culture_group = r"""CASE
      WHEN lower(coalesce(RAW_LAB_NAME,'')) LIKE '%blood culture%' THEN 'blood'
      WHEN lower(coalesce(RAW_LAB_NAME,'')) LIKE '%urine culture%' THEN 'urine'
      WHEN lower(coalesce(RAW_LAB_NAME,'')) LIKE '%respir%' OR lower(coalesce(RAW_LAB_NAME,'')) LIKE '%sputum%' THEN 'respiratory'
      WHEN lower(coalesce(RAW_LAB_NAME,'')) LIKE '%wound%' THEN 'wound'
      WHEN lower(coalesce(RAW_LAB_NAME,'')) LIKE '%tissue%' THEN 'tissue'
      WHEN lower(coalesce(RAW_LAB_NAME,'')) LIKE '%csf%' THEN 'csf'
      WHEN lower(coalesce(RAW_LAB_NAME,'')) LIKE '%fung%' THEN 'fungal'
      WHEN lower(coalesce(RAW_LAB_NAME,'')) LIKE '%afb%' OR lower(coalesce(RAW_LAB_NAME,'')) LIKE '%mycobacter%' THEN 'afb_mycobacterial'
      WHEN lower(coalesce(RAW_LAB_NAME,'')) LIKE '%stool%' THEN 'stool'
      WHEN lower(coalesce(RAW_LAB_NAME,'')) LIKE '%fluid%' THEN 'other_fluid'
      ELSE 'other_culture' END"""

    class_expr = r"""CASE
      WHEN regexp_matches(lower(coalesce(l.RAW_LAB_NAME,'')), '(^|[^a-z])(organism|isolate)([^a-z]|$)') THEN 'organism_or_isolate'
      WHEN lower(coalesce(l.RAW_LAB_NAME,'')) LIKE '%susceptib%'
        OR lower(coalesce(l.RAW_LAB_NAME,'')) LIKE '%minimum inhibitory%'
        OR regexp_matches(lower(coalesce(l.RAW_LAB_NAME,'')), '(^|[^a-z])mic([^a-z]|$)') THEN 'susceptibility_or_mic'
      WHEN lower(coalesce(l.RAW_LAB_NAME,'')) LIKE '%gram stain%'
        OR lower(coalesce(l.RAW_LAB_NAME,'')) LIKE '%gramstain%' THEN 'gram_stain'
      ELSE 'other'
    END"""

    con.execute(f"""
      CREATE TEMP TABLE cultures AS
      SELECT DISTINCT PATID, ENCOUNTERID, SPECIMEN_DATE, SPECIMEN_TIME, LAB_ORDER_DATE,
             {culture_group} AS culture_group
      FROM lab
      WHERE {culture_where}
        AND PATID IS NOT NULL AND ENCOUNTERID IS NOT NULL AND SPECIMEN_DATE IS NOT NULL
    """)

    con.execute(f"""
      CREATE TEMP TABLE candidates AS
      SELECT DISTINCT
        c.PATID, c.ENCOUNTERID, c.SPECIMEN_DATE AS culture_specimen_date,
        c.SPECIMEN_TIME AS culture_specimen_time, c.LAB_ORDER_DATE AS culture_order_date,
        c.culture_group,
        l.SPECIMEN_DATE AS sibling_specimen_date, l.SPECIMEN_TIME AS sibling_specimen_time,
        l.LAB_ORDER_DATE AS sibling_order_date, l.RESULT_DATE AS sibling_result_date,
        {class_expr} AS sibling_class
      FROM cultures c
      JOIN lab l
        ON l.PATID=c.PATID AND l.ENCOUNTERID=c.ENCOUNTERID
       AND l.RESULT_DATE BETWEEN c.SPECIMEN_DATE AND c.SPECIMEN_DATE + INTERVAL 7 DAY
      WHERE NOT ({culture_where})
        AND ({class_expr}) <> 'other'
    """)

    rules = {
        "same_specimen_date_time": "sibling_specimen_date=culture_specimen_date AND sibling_specimen_time=culture_specimen_time AND culture_specimen_time IS NOT NULL AND sibling_specimen_time IS NOT NULL",
        "same_specimen_date": "sibling_specimen_date=culture_specimen_date",
        "same_order_date": "sibling_order_date=culture_order_date AND culture_order_date IS NOT NULL",
        "same_specimen_datetime_or_order_date": "((sibling_specimen_date=culture_specimen_date AND sibling_specimen_time=culture_specimen_time AND culture_specimen_time IS NOT NULL AND sibling_specimen_time IS NOT NULL) OR (sibling_order_date=culture_order_date AND culture_order_date IS NOT NULL))",
    }

    rows = []
    totals = []
    for rule, cond in rules.items():
        df = con.execute(f"""
          SELECT culture_group, sibling_class,
                 count(DISTINCT (PATID, ENCOUNTERID, culture_specimen_date, culture_specimen_time, culture_group))::BIGINT AS count
          FROM candidates
          WHERE {cond}
          GROUP BY 1,2
          ORDER BY 1,2
        """).fetchdf()
        df.insert(0, "linkage_rule", rule)
        rows.append(df)
        tot = con.execute(f"""
          SELECT count(DISTINCT (PATID, ENCOUNTERID, culture_specimen_date, culture_specimen_time, culture_group))::BIGINT AS linked_units
          FROM candidates WHERE {cond}
        """).fetchone()[0]
        totals.append({"linkage_rule": rule, "linked_units": int(tot)})

    out = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    suppress(out).to_csv(args.output_dir / "timestamp_linkage_counts.csv", index=False)
    pd.DataFrame(totals).to_csv(args.output_dir / "timestamp_linkage_totals.csv", index=False)

    completeness = con.execute("""
      SELECT
        count(*)::BIGINT AS culture_units,
        sum(CASE WHEN SPECIMEN_TIME IS NULL THEN 1 ELSE 0 END)::BIGINT AS specimen_time_missing,
        sum(CASE WHEN LAB_ORDER_DATE IS NULL THEN 1 ELSE 0 END)::BIGINT AS lab_order_date_missing
      FROM cultures
    """).fetchdf()
    completeness.to_csv(args.output_dir / "culture_timestamp_completeness.csv", index=False)

    summary = {
        "privacy_mode": "aggregate_only_internal_identifier_linkage_no_ids_no_rows_no_raw_result_text_export",
        "minimum_reported_cell": MIN_CELL,
        "purpose": "Test exact specimen timestamp and order-date linkage as structural keys for PSU culture positivity components.",
        "guardrail": "Diagnostic only; no positivity phenotype is frozen from this audit alone.",
    }
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
