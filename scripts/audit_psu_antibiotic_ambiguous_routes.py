#!/usr/bin/env python3
"""Aggregate-only audit of route ambiguity for vancomycin, linezolid, and aztreonam.

Uses local medication names plus corrected legacy RxNorm inclusion codes to classify
ambiguous agents as injectable/IV, oral, inhaled, or unspecified. Exports aggregate
counts only; no patient identifiers, rows, or free-text values are exported.
"""
from __future__ import annotations

import argparse
import ast
import json
import re
from pathlib import Path

import duckdb
import pandas as pd

MIN_CELL = 11
AMBIGUOUS = {
    "vancomycin": r"vancomycin|vancocin",
    "linezolid": r"linezolid|zyvox",
    "aztreonam": r"aztreonam|azactam|cayston",
}
FORMULATION_PATTERNS = {
    "injectable_or_iv": r"inject|intravenous|\biv\b|vial|infusion|premix|piggyback|solution for injection|injection solution",
    "oral": r"oral|capsule|tablet|po\b|by mouth|suspension|oral solution",
    "inhaled": r"inhal|nebul|cayston",
}


def q(value: str) -> str:
    return "'" + value.replace("'", "''") + "'"


def parse_legacy_codes(path: Path) -> tuple[dict[str, set[str]], set[str]]:
    text = path.read_text(encoding="utf-8", errors="ignore")
    tree = ast.parse(text)
    assignments: dict[str, object] = {}
    for node in tree.body:
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name):
                    try:
                        assignments[target.id] = ast.literal_eval(node.value)
                    except Exception:
                        pass
    include: dict[str, set[str]] = {}
    for key, value in assignments.items():
        if not isinstance(value, (list, tuple)):
            continue
        codes: set[str] = set()
        for item in value:
            if isinstance(item, (list, tuple)) and item:
                token = str(item[0]).strip()
                if re.fullmatch(r"\d{4,9}", token):
                    codes.add(token)
        if codes:
            include[key] = codes
    excluded: set[str] = set()
    for key, value in assignments.items():
        if "exclude" in key.lower() and isinstance(value, (set, list, tuple)):
            for item in value:
                token = str(item).strip()
                if re.fullmatch(r"\d{4,9}", token):
                    excluded.add(token)
    return include, excluded


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
    legacy = root / "code" / "config" / "codes_antibiotics.py"
    if not prescribing.exists() or not med_admin.exists() or not legacy.exists():
        raise FileNotFoundError("Required PSU prescribing, med_admin, or legacy antibiotic code missing")

    include_map, excluded = parse_legacy_codes(legacy)
    agent_codes: dict[str, set[str]] = {}
    for agent in AMBIGUOUS:
        codes = set()
        for key, vals in include_map.items():
            if agent in key.lower():
                codes |= vals
        agent_codes[agent] = codes - excluded

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
    for source, view, name_col, code_col, route_col, enc_col in specs:
        for agent, text_pat in AMBIGUOUS.items():
            codes = sorted(agent_codes.get(agent, set()))
            code_sql = ",".join(q(c) for c in codes) or "''"
            name_expr = f"lower(coalesce(cast({name_col} AS VARCHAR),''))"
            route_expr = f"upper(coalesce(cast({route_col} AS VARCHAR),''))"
            base = (
                f"(cast({code_col} AS VARCHAR) IN ({code_sql}) OR "
                f"regexp_matches({name_expr}, {q(text_pat)}))"
            )
            cls = (
                f"CASE "
                f"WHEN regexp_matches({name_expr}, {q(FORMULATION_PATTERNS['inhaled'])}) OR {route_expr} IN ('RESPIRATORY_TRACT','INHALATION') THEN 'inhaled' "
                f"WHEN regexp_matches({name_expr}, {q(FORMULATION_PATTERNS['oral'])}) OR {route_expr}='ORAL' THEN 'oral' "
                f"WHEN regexp_matches({name_expr}, {q(FORMULATION_PATTERNS['injectable_or_iv'])}) OR {route_expr}='INTRAVENOUS' THEN 'injectable_or_iv' "
                f"ELSE 'unspecified' END"
            )
            df = con.execute(f"""
                SELECT {q(source)} AS source, {q(agent)} AS agent,
                       {cls} AS formulation_class,
                       count(*)::BIGINT AS count
                FROM {view}
                WHERE {base}
                GROUP BY 1,2,3
                ORDER BY count DESC
            """).fetchdf()
            rows.append(suppress(df))

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

    pd.concat(rows, ignore_index=True).to_csv(args.output_dir / "ambiguous_agent_formulation_counts.csv", index=False)
    pd.concat(encounter_rows, ignore_index=True).to_csv(args.output_dir / "ambiguous_agent_encounter_counts.csv", index=False)

    code_rows = []
    for agent, codes in agent_codes.items():
        for code in sorted(codes):
            code_rows.append({"agent": agent, "legacy_include_code": code})
    pd.DataFrame(code_rows).to_csv(args.output_dir / "ambiguous_agent_legacy_codes.csv", index=False)

    summary = {
        "privacy_mode": "aggregate_only_no_ids_no_patient_rows_no_free_text_export",
        "minimum_reported_cell": MIN_CELL,
        "agents": list(AMBIGUOUS),
        "purpose": "Resolve route/formulation ambiguity for vancomycin, linezolid, and aztreonam before freezing the PSU broad-spectrum antibiotic phenotype.",
        "classification_priority": ["inhaled", "oral", "injectable_or_iv", "unspecified"],
        "guardrail": "Diagnostic only; do not freeze the final PSU antibiotic phenotype unless ambiguous-agent formulation results are clinically coherent.",
    }
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
