#!/usr/bin/env python3
"""Aggregate-only audit of route ambiguity for vancomycin, linezolid, and aztreonam.

Uses conservative drug-name matching plus known legacy RxNorm include/exclude codes.
Exports aggregate counts only; no patient identifiers, rows, or free-text values.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import duckdb
import pandas as pd

MIN_CELL = 11
AGENTS = {
    "vancomycin": {
        "name_terms": ["vancomycin", "vancocin"],
        "include_codes": ["11124", "202368", "239209"],
        "exclude_codes": ["313570", "313571", "2000134"],
    },
    "linezolid": {
        "name_terms": ["linezolid", "zyvox"],
        "include_codes": ["190376", "261710"],
        "exclude_codes": [],
    },
    "aztreonam": {
        "name_terms": ["aztreonam", "azactam", "cayston"],
        "include_codes": ["1272", "202561"],
        "exclude_codes": [],
    },
}


def q(value: str) -> str:
    return "'" + value.replace("'", "''") + "'"


def suppress(df: pd.DataFrame, col: str = "count") -> pd.DataFrame:
    out = df.copy()
    if col in out.columns:
        mask = out[col].fillna(0) < MIN_CELL
        out.loc[mask, col] = pd.NA
        out["suppressed"] = mask
    return out


def sql_in(values: list[str]) -> str:
    return ",".join(q(x) for x in values) or "''"


def name_condition(name_expr: str, terms: list[str]) -> str:
    return "(" + " OR ".join(f"strpos({name_expr}, {q(t)}) > 0" for t in terms) + ")"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("data_root", type=Path)
    ap.add_argument("--output-dir", type=Path, required=True)
    args = ap.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    root = args.data_root / "PCORnet" / "parquet"
    prescribing = root / "prescribing.parquet"
    med_admin = root / "med_admin.parquet"
    if not prescribing.exists() or not med_admin.exists():
        raise FileNotFoundError("Required PSU prescribing or med_admin parquet missing")

    con = duckdb.connect()
    con.execute("PRAGMA threads=4")
    con.execute(f"CREATE VIEW p AS SELECT * FROM read_parquet({q(str(prescribing))})")
    con.execute(f"CREATE VIEW m AS SELECT * FROM read_parquet({q(str(med_admin))})")

    specs = [
        ("PRESCRIBING", "p", "raw_rx_med_name", "rxnorm_cui", "rx_route", "encounterid"),
        ("MED_ADMIN", "m", "raw_medadmin_med_name", "medadmin_code", "medadmin_route", "encounterid"),
    ]

    rows = []
    encounter_rows = []
    sanity_rows = []
    code_rows = []
    for agent, spec in AGENTS.items():
        for code in spec["include_codes"]:
            code_rows.append({"agent": agent, "mapping_role": "include", "legacy_code": code})
        for code in spec["exclude_codes"]:
            code_rows.append({"agent": agent, "mapping_role": "exclude", "legacy_code": code})

    for source, view, name_col, code_col, route_col, enc_col in specs:
        total_rows = int(con.execute(f"SELECT count(*) FROM {view}").fetchone()[0])
        for agent, spec in AGENTS.items():
            name_expr = f"lower(coalesce(cast({name_col} AS VARCHAR),''))"
            route_expr = f"upper(trim(coalesce(cast({route_col} AS VARCHAR),'')))"
            code_expr = f"trim(coalesce(cast({code_col} AS VARCHAR),''))"
            text_match = name_condition(name_expr, spec["name_terms"])
            include_codes = sql_in(spec["include_codes"])
            exclude_codes = sql_in(spec["exclude_codes"])
            code_match = f"{code_expr} IN ({include_codes})"
            code_excluded = f"{code_expr} IN ({exclude_codes})" if spec["exclude_codes"] else "FALSE"
            base = f"(({text_match}) OR ({code_match})) AND NOT ({code_excluded})"

            inhaled = f"(strpos({name_expr}, 'inhal')>0 OR strpos({name_expr}, 'nebul')>0 OR strpos({name_expr}, 'cayston')>0 OR {route_expr} IN ('RESPIRATORY_TRACT','INHALATION'))"
            oral = f"(strpos({name_expr}, 'oral')>0 OR strpos({name_expr}, 'capsule')>0 OR strpos({name_expr}, 'tablet')>0 OR strpos({name_expr}, 'by mouth')>0 OR strpos({name_expr}, 'oral solution')>0 OR {route_expr}='ORAL')"
            injectable = f"(strpos({name_expr}, 'inject')>0 OR strpos({name_expr}, 'intravenous')>0 OR strpos({name_expr}, ' vial')>0 OR strpos({name_expr}, 'infusion')>0 OR strpos({name_expr}, 'premix')>0 OR strpos({name_expr}, 'piggyback')>0 OR {route_expr}='INTRAVENOUS')"
            cls = f"CASE WHEN {inhaled} THEN 'inhaled' WHEN {oral} THEN 'oral' WHEN {injectable} THEN 'injectable_or_iv' ELSE 'unspecified' END"

            counts = con.execute(f"""
                SELECT {q(source)} AS source, {q(agent)} AS agent,
                       {cls} AS formulation_class, count(*)::BIGINT AS count
                FROM {view}
                WHERE {base}
                GROUP BY 1,2,3
                ORDER BY count DESC
            """).fetchdf()
            rows.append(suppress(counts))

            enc = con.execute(f"""
                SELECT {q(source)} AS source, {q(agent)} AS agent,
                       count(DISTINCT CASE WHEN {base} THEN {enc_col} END)::BIGINT AS any_agent_encounters,
                       count(DISTINCT CASE WHEN {base} AND ({cls})='injectable_or_iv' THEN {enc_col} END)::BIGINT AS injectable_or_iv_encounters,
                       count(DISTINCT CASE WHEN {base} AND ({cls})='oral' THEN {enc_col} END)::BIGINT AS oral_encounters,
                       count(DISTINCT CASE WHEN {base} AND ({cls})='inhaled' THEN {enc_col} END)::BIGINT AS inhaled_encounters,
                       count(DISTINCT CASE WHEN {base} AND ({cls})='unspecified' THEN {enc_col} END)::BIGINT AS unspecified_encounters
                FROM {view}
            """).fetchdf()
            encounter_rows.append(enc)

            sanity = con.execute(f"""
                SELECT
                  count(*) FILTER (WHERE {text_match})::BIGINT AS text_match_rows,
                  count(*) FILTER (WHERE {code_match})::BIGINT AS include_code_rows,
                  count(*) FILTER (WHERE {code_excluded})::BIGINT AS explicit_excluded_code_rows,
                  count(*) FILTER (WHERE {base})::BIGINT AS final_candidate_rows
                FROM {view}
            """).fetchdf().iloc[0].to_dict()
            sanity_rows.append({"source": source, "agent": agent, "source_total_rows": total_rows, **{k:int(v) for k,v in sanity.items()}})

    pd.concat(rows, ignore_index=True).to_csv(args.output_dir / "ambiguous_agent_formulation_counts.csv", index=False)
    pd.concat(encounter_rows, ignore_index=True).to_csv(args.output_dir / "ambiguous_agent_encounter_counts.csv", index=False)
    pd.DataFrame(code_rows).to_csv(args.output_dir / "ambiguous_agent_legacy_codes.csv", index=False)
    pd.DataFrame(sanity_rows).to_csv(args.output_dir / "ambiguous_agent_sanity_counts.csv", index=False)

    summary = {
        "privacy_mode": "aggregate_only_no_ids_no_patient_rows_no_free_text_export",
        "minimum_reported_cell": MIN_CELL,
        "agents": list(AGENTS),
        "correction": "Prior run produced implausibly large counts and an empty legacy-code artifact. This run uses explicit reviewed include/exclude codes and strpos-based drug-name matching with source-total sanity counts.",
        "classification_priority": ["inhaled", "oral", "injectable_or_iv", "unspecified"],
        "guardrail": "Diagnostic only; freeze the final PSU antibiotic phenotype only if candidate totals and formulation distributions are clinically coherent.",
    }
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
