#!/usr/bin/env python3
"""Refined aggregate-only audit of PSU culture-adjacent microbiology components.

This version excludes common false positives from broad substring matching (for
example urinalysis "microscopic" components and "high sensitivity" CRP) and
focuses on clinically plausible microbiology result components. Identifiers and
RAW_RESULT are used only locally. Exports contain aggregate metadata/counts only.
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

    # Strict component classes. Do not use generic '%micro%' or '%sensitivity%'.
    class_expr = r"""CASE
      WHEN regexp_matches(lower(coalesce(l.RAW_LAB_NAME,'')), '(^|[^a-z])(organism|isolate|identification)([^a-z]|$)')
        THEN 'organism_or_identification'
      WHEN lower(coalesce(l.RAW_LAB_NAME,'')) LIKE '%susceptib%'
        OR lower(coalesce(l.RAW_LAB_NAME,'')) LIKE '%minimum inhibitory%'
        OR regexp_matches(lower(coalesce(l.RAW_LAB_NAME,'')), '(^|[^a-z])mic([^a-z]|$)')
        THEN 'susceptibility_or_mic'
      WHEN lower(coalesce(l.RAW_LAB_NAME,'')) LIKE '%gram stain%'
        OR lower(coalesce(l.RAW_LAB_NAME,'')) LIKE '%gramstain%'
        THEN 'gram_stain'
      WHEN lower(coalesce(l.RAW_LAB_NAME,'')) LIKE '%afb%'
        OR lower(coalesce(l.RAW_LAB_NAME,'')) LIKE '%mycobacter%'
        THEN 'afb_related'
      WHEN lower(coalesce(l.RAW_LAB_NAME,'')) LIKE '%fungal%'
        OR lower(coalesce(l.RAW_LAB_NAME,'')) LIKE '%fungus%'
        OR lower(coalesce(l.RAW_LAB_NAME,'')) LIKE '%yeast%'
        THEN 'fungal_related'
      WHEN lower(coalesce(l.RAW_LAB_NAME,'')) LIKE '%bacterial culture%'
        OR lower(coalesce(l.RAW_LAB_NAME,'')) LIKE '%bacteria culture%'
        THEN 'bacterial_culture_component'
      ELSE 'not_specific_microbiology'
    END"""

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
        {class_expr} AS sibling_class,
        datediff('day', c.SPECIMEN_DATE, l.RESULT_DATE) AS day_offset
      FROM culture_units c
      JOIN lab l
        ON l.PATID=c.PATID
       AND l.ENCOUNTERID=c.ENCOUNTERID
       AND l.RESULT_DATE BETWEEN c.SPECIMEN_DATE AND c.SPECIMEN_DATE + INTERVAL 2 DAY
      WHERE NOT ({culture_where})
    """)

    specific = "sibling_class <> 'not_specific_microbiology'"

    counts = con.execute(f"""
      SELECT culture_group, sibling_class, day_offset,
             count(*)::BIGINT AS count,
             count(DISTINCT (PATID, ENCOUNTERID, culture_specimen_date, culture_group))::BIGINT AS culture_units_with_component,
             sum(raw_result_present)::BIGINT AS rows_with_raw_result
      FROM siblings
      WHERE {specific}
      GROUP BY 1,2,3
      ORDER BY culture_group, sibling_class, day_offset
    """).fetchdf()
    suppress(counts).to_csv(args.output_dir / "refined_sibling_class_counts.csv", index=False)

    comps = con.execute(f"""
      SELECT culture_group, sibling_class, loinc, raw_lab_name,
             count(*)::BIGINT AS count,
             count(DISTINCT (PATID, ENCOUNTERID, culture_specimen_date, culture_group))::BIGINT AS culture_units_with_component,
             sum(raw_result_present)::BIGINT AS rows_with_raw_result
      FROM siblings
      WHERE {specific}
      GROUP BY 1,2,3,4
      HAVING count(*) >= {MIN_CELL}
      ORDER BY count DESC
      LIMIT 500
    """).fetchdf()
    comps["suppressed"] = False
    comps.to_csv(args.output_dir / "refined_top_components.csv", index=False)

    qual = con.execute(f"""
      SELECT sibling_class,
             CASE WHEN result_qual='' THEN '<EMPTY>' ELSE result_qual END AS result_qual,
             count(*)::BIGINT AS count
      FROM siblings
      WHERE {specific}
      GROUP BY 1,2
      ORDER BY 1, count DESC
    """).fetchdf()
    suppress(qual).to_csv(args.output_dir / "refined_result_qual_counts.csv", index=False)

    totals = con.execute(f"""
      SELECT
        (SELECT count(*) FROM culture_units)::BIGINT AS culture_units,
        count(*)::BIGINT AS specific_sibling_rows,
        count(DISTINCT (PATID, ENCOUNTERID, culture_specimen_date, culture_group))::BIGINT AS culture_units_with_specific_sibling,
        count(DISTINCT CASE WHEN sibling_class='organism_or_identification' THEN (PATID, ENCOUNTERID, culture_specimen_date, culture_group) END)::BIGINT AS units_with_organism,
        count(DISTINCT CASE WHEN sibling_class='susceptibility_or_mic' THEN (PATID, ENCOUNTERID, culture_specimen_date, culture_group) END)::BIGINT AS units_with_susceptibility,
        count(DISTINCT CASE WHEN sibling_class='gram_stain' THEN (PATID, ENCOUNTERID, culture_specimen_date, culture_group) END)::BIGINT AS units_with_gram_stain
      FROM siblings
      WHERE {specific}
    """).fetchdf().iloc[0].to_dict()
    totals = {k: int(v) for k, v in totals.items()}

    summary = {
        "privacy_mode": "aggregate_only_internal_identifier_and_raw_result_presence_inspection_no_ids_no_rows_no_raw_result_text_export",
        "minimum_reported_cell": MIN_CELL,
        "window": "same PATID+ENCOUNTERID; sibling RESULT_DATE from culture SPECIMEN_DATE through +2 days",
        "culture_unit": "distinct PATID+ENCOUNTERID+SPECIMEN_DATE+culture_group",
        "refinement": "Removed generic micro/sensitivity substring rules that misclassified urinalysis microscopic components and high-sensitivity CRP.",
        "totals": totals,
        "guardrail": "Diagnostic only; specific sibling components are not yet the frozen positive-culture phenotype.",
    }
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
