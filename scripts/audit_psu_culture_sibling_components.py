#!/usr/bin/env python3
"""Aggregate-only audit of lab components surrounding PSU clinical-culture tests.

This diagnostic searches for non-culture lab rows in the same patient+encounter and
near the culture specimen date to identify where microbiology positivity may be
encoded. Identifiers and RAW_RESULT are used only locally. Exports contain only
aggregate metadata/counts with small-cell suppression and no raw result text.
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

    con.execute(f"""
      CREATE TEMP TABLE culture_units AS
      SELECT DISTINCT PATID, ENCOUNTERID, SPECIMEN_DATE,
             {culture_group} AS culture_group
      FROM lab
      WHERE {culture_where}
        AND PATID IS NOT NULL AND ENCOUNTERID IS NOT NULL AND SPECIMEN_DATE IS NOT NULL
    """)

    # Candidate sibling rows: same encounter, result date from specimen day through +2d,
    # excluding rows that are themselves parent culture tests. This deliberately starts
    # broad so we can see what metadata components actually accompany cultures.
    con.execute(f"""
      CREATE TEMP TABLE siblings AS
      SELECT DISTINCT
        c.PATID, c.ENCOUNTERID, c.SPECIMEN_DATE AS culture_specimen_date,
        c.culture_group,
        coalesce(l.LAB_LOINC,'<MISSING>') AS loinc,
        coalesce(l.RAW_LAB_NAME,'<MISSING>') AS raw_lab_name,
        coalesce(l.SPECIMEN_SOURCE,'<MISSING>') AS specimen_source,
        coalesce(l.RESULT_QUAL,'<MISSING>') AS result_qual,
        CASE WHEN l.RAW_RESULT IS NOT NULL AND trim(cast(l.RAW_RESULT AS VARCHAR)) <> '' THEN 1 ELSE 0 END AS raw_result_present,
        CASE
          WHEN lower(coalesce(l.RAW_LAB_NAME,'')) LIKE '%susceptib%'
            OR lower(coalesce(l.RAW_LAB_NAME,'')) LIKE '%sensitivity%'
            OR lower(coalesce(l.RAW_LAB_NAME,'')) LIKE '%minimum inhibitory%'
            OR lower(coalesce(l.RAW_LAB_NAME,'')) LIKE '% mic %' THEN 'susceptibility_or_mic'
          WHEN lower(coalesce(l.RAW_LAB_NAME,'')) LIKE '%organism%'
            OR lower(coalesce(l.RAW_LAB_NAME,'')) LIKE '%isolate%'
            OR lower(coalesce(l.RAW_LAB_NAME,'')) LIKE '%identification%' THEN 'organism_or_identification'
          WHEN lower(coalesce(l.RAW_LAB_NAME,'')) LIKE '%gram stain%'
            OR lower(coalesce(l.RAW_LAB_NAME,'')) LIKE '%gramstain%'
            OR lower(coalesce(l.RAW_LAB_NAME,'')) LIKE '%gram %' THEN 'gram_stain'
          WHEN lower(coalesce(l.RAW_LAB_NAME,'')) LIKE '%culture%' THEN 'other_culture_component'
          WHEN lower(coalesce(l.RAW_LAB_NAME,'')) LIKE '%micro%'
            OR lower(coalesce(l.RAW_LAB_NAME,'')) LIKE '%bacter%'
            OR lower(coalesce(l.RAW_LAB_NAME,'')) LIKE '%yeast%'
            OR lower(coalesce(l.RAW_LAB_NAME,'')) LIKE '%fung%'
            OR lower(coalesce(l.RAW_LAB_NAME,'')) LIKE '%afb%'
            OR lower(coalesce(l.RAW_LAB_NAME,'')) LIKE '%mycobacter%' THEN 'other_microbiology_named'
          ELSE 'other_lab'
        END AS sibling_class,
        datediff('day', c.SPECIMEN_DATE, l.RESULT_DATE) AS day_offset
      FROM culture_units c
      JOIN lab l
        ON l.PATID=c.PATID
       AND l.ENCOUNTERID=c.ENCOUNTERID
       AND l.RESULT_DATE BETWEEN c.SPECIMEN_DATE AND c.SPECIMEN_DATE + INTERVAL 2 DAY
      WHERE NOT ({culture_where})
    """)

    class_counts = con.execute("""
      SELECT culture_group, sibling_class, day_offset,
             count(*)::BIGINT AS count,
             count(DISTINCT (PATID, ENCOUNTERID, culture_specimen_date, culture_group))::BIGINT AS culture_units_with_component,
             sum(raw_result_present)::BIGINT AS rows_with_raw_result
      FROM siblings
      GROUP BY 1,2,3
      ORDER BY culture_group, sibling_class, day_offset
    """).fetchdf()
    class_counts = suppress(class_counts)
    class_counts.to_csv(args.output_dir / "sibling_class_counts.csv", index=False)

    top_components = con.execute("""
      SELECT culture_group, sibling_class, loinc, raw_lab_name,
             count(*)::BIGINT AS count,
             count(DISTINCT (PATID, ENCOUNTERID, culture_specimen_date, culture_group))::BIGINT AS culture_units_with_component,
             sum(raw_result_present)::BIGINT AS rows_with_raw_result
      FROM siblings
      WHERE sibling_class <> 'other_lab'
      GROUP BY 1,2,3,4
      HAVING count(*) >= 11
      ORDER BY count DESC
      LIMIT 500
    """).fetchdf()
    top_components["suppressed"] = False
    top_components.to_csv(args.output_dir / "top_microbiology_sibling_components.csv", index=False)

    result_qual = con.execute("""
      SELECT sibling_class,
             CASE WHEN result_qual='' THEN '<EMPTY>' ELSE result_qual END AS result_qual,
             count(*)::BIGINT AS count
      FROM siblings
      WHERE sibling_class <> 'other_lab'
      GROUP BY 1,2
      ORDER BY 1, count DESC
    """).fetchdf()
    result_qual = suppress(result_qual)
    result_qual.to_csv(args.output_dir / "sibling_result_qual_counts.csv", index=False)

    totals = con.execute("""
      SELECT
        (SELECT count(*) FROM culture_units)::BIGINT AS culture_units,
        count(*)::BIGINT AS sibling_rows,
        count(DISTINCT (PATID, ENCOUNTERID, culture_specimen_date, culture_group))::BIGINT AS culture_units_with_any_sibling,
        count(DISTINCT CASE WHEN sibling_class <> 'other_lab' THEN (PATID, ENCOUNTERID, culture_specimen_date, culture_group) END)::BIGINT AS culture_units_with_micro_named_sibling,
        count(DISTINCT CASE WHEN sibling_class IN ('organism_or_identification','susceptibility_or_mic','gram_stain') THEN (PATID, ENCOUNTERID, culture_specimen_date, culture_group) END)::BIGINT AS culture_units_with_specific_micro_signal
      FROM siblings
    """).fetchdf().iloc[0].to_dict()
    totals = {k: int(v) for k, v in totals.items()}

    summary = {
        "privacy_mode": "aggregate_only_internal_identifier_and_raw_result_presence_inspection_no_ids_no_rows_no_raw_result_text_export",
        "minimum_reported_cell": MIN_CELL,
        "window": "same PATID+ENCOUNTERID; sibling RESULT_DATE from culture SPECIMEN_DATE through +2 days",
        "culture_unit": "distinct PATID+ENCOUNTERID+SPECIMEN_DATE+culture_group",
        "totals": totals,
        "purpose": "Identify the actual non-parent lab components that carry PSU microbiology result information before freezing culture positivity semantics.",
        "guardrail": "This is diagnostic only. No sibling component is automatically treated as a positive culture without clinical plausibility review.",
    }
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
