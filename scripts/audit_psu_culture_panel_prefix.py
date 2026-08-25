#!/usr/bin/env python3
"""Aggregate-only audit of PSU culture panel-prefix linkage.

Uses the RAW_LAB_NAME prefix before the first pipe as a local assay/panel family key.
Identifiers and raw result text remain local. Exports only aggregate counts and
component metadata with small-cell suppression.
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

    prefix_expr = r"""lower(trim(CASE
      WHEN strpos(coalesce(RAW_LAB_NAME,''), '|') > 0
        THEN split_part(coalesce(RAW_LAB_NAME,''), '|', 1)
      ELSE coalesce(RAW_LAB_NAME,'') END))"""

    con.execute(f"""
      CREATE TEMP TABLE cultures AS
      SELECT DISTINCT PATID, ENCOUNTERID, SPECIMEN_DATE,
             {culture_group} AS culture_group,
             {prefix_expr} AS panel_prefix
      FROM lab
      WHERE {culture_where}
        AND PATID IS NOT NULL AND ENCOUNTERID IS NOT NULL AND SPECIMEN_DATE IS NOT NULL
    """)

    con.execute(f"""
      CREATE TEMP TABLE linked AS
      SELECT DISTINCT
        c.PATID, c.ENCOUNTERID, c.SPECIMEN_DATE AS culture_specimen_date,
        c.culture_group, c.panel_prefix,
        coalesce(l.LAB_LOINC,'<MISSING>') AS loinc,
        coalesce(l.RAW_LAB_NAME,'<MISSING>') AS raw_lab_name,
        {prefix_expr.replace('RAW_LAB_NAME', 'l.RAW_LAB_NAME')} AS sibling_prefix,
        coalesce(l.RESULT_QUAL,'<MISSING>') AS result_qual,
        CASE WHEN l.RAW_RESULT IS NOT NULL AND trim(cast(l.RAW_RESULT AS VARCHAR)) <> '' THEN 1 ELSE 0 END AS raw_result_present,
        datediff('day', c.SPECIMEN_DATE, l.RESULT_DATE) AS day_offset,
        CASE
          WHEN regexp_matches(lower(coalesce(l.RAW_LAB_NAME,'')), '(^|[^a-z])(organism|isolate)([^a-z]|$)') THEN 'organism_or_isolate'
          WHEN lower(coalesce(l.RAW_LAB_NAME,'')) LIKE '%susceptib%'
            OR lower(coalesce(l.RAW_LAB_NAME,'')) LIKE '%minimum inhibitory%'
            OR regexp_matches(lower(coalesce(l.RAW_LAB_NAME,'')), '(^|[^a-z])mic([^a-z]|$)') THEN 'susceptibility_or_mic'
          WHEN lower(coalesce(l.RAW_LAB_NAME,'')) LIKE '%gram stain%' OR lower(coalesce(l.RAW_LAB_NAME,'')) LIKE '%gramstain%' THEN 'gram_stain'
          ELSE 'other_same_panel_component'
        END AS component_class
      FROM cultures c
      JOIN lab l
        ON l.PATID=c.PATID
       AND l.ENCOUNTERID=c.ENCOUNTERID
       AND l.RESULT_DATE BETWEEN c.SPECIMEN_DATE AND c.SPECIMEN_DATE + INTERVAL 2 DAY
       AND {prefix_expr.replace('RAW_LAB_NAME', 'l.RAW_LAB_NAME')} = c.panel_prefix
      WHERE c.panel_prefix <> ''
        AND NOT ({culture_where.replace('RAW_LAB_NAME', 'l.RAW_LAB_NAME')})
    """)

    counts = con.execute("""
      SELECT culture_group, component_class, day_offset,
             count(*)::BIGINT AS count,
             count(DISTINCT (PATID, ENCOUNTERID, culture_specimen_date, culture_group, panel_prefix))::BIGINT AS culture_units_with_component,
             sum(raw_result_present)::BIGINT AS rows_with_raw_result
      FROM linked
      GROUP BY 1,2,3
      ORDER BY 1,2,3
    """).fetchdf()
    suppress(counts).to_csv(args.output_dir / "panel_prefix_class_counts.csv", index=False)

    comps = con.execute(f"""
      SELECT culture_group, panel_prefix, component_class, loinc, raw_lab_name,
             count(*)::BIGINT AS count,
             count(DISTINCT (PATID, ENCOUNTERID, culture_specimen_date, culture_group, panel_prefix))::BIGINT AS culture_units_with_component,
             sum(raw_result_present)::BIGINT AS rows_with_raw_result
      FROM linked
      GROUP BY 1,2,3,4,5
      HAVING count(*) >= {MIN_CELL}
      ORDER BY count DESC
      LIMIT 500
    """).fetchdf()
    comps["suppressed"] = False
    comps.to_csv(args.output_dir / "panel_prefix_top_components.csv", index=False)

    prefixes = con.execute(f"""
      SELECT culture_group, panel_prefix, count(*)::BIGINT AS culture_units
      FROM cultures
      WHERE panel_prefix <> ''
      GROUP BY 1,2
      HAVING count(*) >= {MIN_CELL}
      ORDER BY culture_units DESC
      LIMIT 300
    """).fetchdf()
    prefixes["suppressed"] = False
    prefixes.to_csv(args.output_dir / "culture_panel_prefix_counts.csv", index=False)

    totals = con.execute("""
      SELECT
        (SELECT count(*) FROM cultures)::BIGINT AS culture_units,
        (SELECT count(*) FROM cultures WHERE panel_prefix <> '')::BIGINT AS culture_units_with_prefix,
        count(*)::BIGINT AS linked_rows,
        count(DISTINCT (PATID, ENCOUNTERID, culture_specimen_date, culture_group, panel_prefix))::BIGINT AS culture_units_with_same_prefix_component,
        count(DISTINCT CASE WHEN component_class='organism_or_isolate' THEN (PATID, ENCOUNTERID, culture_specimen_date, culture_group, panel_prefix) END)::BIGINT AS units_with_organism_or_isolate,
        count(DISTINCT CASE WHEN component_class='susceptibility_or_mic' THEN (PATID, ENCOUNTERID, culture_specimen_date, culture_group, panel_prefix) END)::BIGINT AS units_with_susceptibility_or_mic,
        count(DISTINCT CASE WHEN component_class='gram_stain' THEN (PATID, ENCOUNTERID, culture_specimen_date, culture_group, panel_prefix) END)::BIGINT AS units_with_gram_stain
      FROM linked
    """).fetchdf().iloc[0].to_dict()
    totals = {k: int(v) for k, v in totals.items()}

    summary = {
        "privacy_mode": "aggregate_only_internal_identifier_linkage_no_ids_no_rows_no_raw_result_text_export",
        "minimum_reported_cell": MIN_CELL,
        "window": "same PATID+ENCOUNTERID, sibling RESULT_DATE 0-2 days after culture SPECIMEN_DATE",
        "linkage_key": "normalized RAW_LAB_NAME prefix before first pipe, exact match",
        "purpose": "Test whether assay/panel-prefix linkage isolates true culture sibling components more cleanly than broad same-encounter matching.",
        "totals": totals,
        "guardrail": "Diagnostic only; no positivity phenotype is frozen from this audit alone.",
    }
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
