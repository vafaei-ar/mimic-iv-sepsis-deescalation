#!/usr/bin/env python3
"""Aggregate-only PSU antibiotic mapping/source audit.

Compares the frozen MIMIC broad-spectrum name phenotype with the existing PSU
antibiotic code source and with PCORnet PRESCRIBING versus MED_ADMIN. Patient
identifiers are used only for aggregate encounter-set comparison. No patient rows,
identifiers, or free-text medication records are exported.
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
BROAD_AGENTS = {
    "vancomycin": r"vancomycin",
    "piperacillin_tazobactam": r"piperacillin|tazobactam|zosyn",
    "cefepime": r"cefepime",
    "ceftazidime": r"ceftazidime",
    "meropenem": r"meropenem",
    "imipenem": r"imipenem",
    "aztreonam": r"aztreonam",
    "linezolid": r"linezolid",
    "daptomycin": r"daptomycin",
    "ceftolozane": r"ceftolozane",
    "avibactam": r"avibactam",
}
BROAD_PATTERN = "|".join(BROAD_AGENTS.values())
IV_ROUTE_PATTERN = r"(^|[^a-z])(iv|intravenous)([^a-z]|$)"


def qident(name: str) -> str:
    return '"' + name.replace('"', '""') + '"'


def qstr(value: str) -> str:
    return "'" + value.replace("'", "''") + "'"


def first_present(cols: set[str], candidates: list[str]) -> str | None:
    upper = {c.upper(): c for c in cols}
    for cand in candidates:
        if cand.upper() in upper:
            return upper[cand.upper()]
    return None


def find_parquet(root: Path, stem: str) -> Path:
    candidates = sorted(root.glob(f"PCORnet/parquet/**/{stem}*.parquet"))
    if not candidates:
        candidates = sorted(root.glob(f"PCORnet/parquet/**/{stem.replace('_','')}*.parquet"))
    if not candidates:
        raise FileNotFoundError(f"No parquet found for {stem}")
    # Prefer exact basename if present, otherwise the largest likely consolidated file.
    exact = [p for p in candidates if p.stem.lower() == stem.lower()]
    if exact:
        return exact[0]
    return max(candidates, key=lambda p: p.stat().st_size)


def safe_count(v: int) -> int | None:
    return int(v) if int(v) >= MIN_CELL else None


def old_code_evidence(path: Path) -> dict:
    out = {"path": str(path), "exists": path.exists(), "agent_mentions": {}, "rxnorm_literal_count": 0}
    if not path.exists():
        return out
    text = path.read_text(encoding="utf-8", errors="ignore")
    low = text.lower()
    out["agent_mentions"] = {agent: bool(re.search(pattern, low)) for agent, pattern in BROAD_AGENTS.items()}
    try:
        tree = ast.parse(text)
        numeric = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Constant):
                if isinstance(node.value, int) and 1000 <= node.value <= 999999999:
                    numeric.add(str(node.value))
                elif isinstance(node.value, str):
                    for token in re.findall(r"\b\d{4,9}\b", node.value):
                        numeric.add(token)
        out["rxnorm_literal_count"] = len(numeric)
    except SyntaxError:
        out["rxnorm_literal_count"] = len(set(re.findall(r"\b\d{4,9}\b", text)))
    out["mentions_rxnorm"] = "rxnorm" in low
    out["mentions_med_admin"] = "med_admin" in low or "medadmin" in low
    out["mentions_prescribing"] = "prescribing" in low
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("data_root", type=Path)
    ap.add_argument("--output-dir", type=Path, required=True)
    args = ap.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    prescribing = find_parquet(args.data_root, "prescribing")
    med_admin = find_parquet(args.data_root, "med_admin")
    old_code = args.data_root / "PCORnet" / "code" / "config" / "codes_antibiotics.py"

    con = duckdb.connect()
    con.execute("PRAGMA threads=4")
    con.execute(f"CREATE VIEW prescribing AS SELECT * FROM read_parquet({qstr(str(prescribing))})")
    con.execute(f"CREATE VIEW med_admin AS SELECT * FROM read_parquet({qstr(str(med_admin))})")

    rows = []
    agent_rows = []
    source_specs = {
        "PRESCRIBING": {
            "view": "prescribing",
            "text_candidates": ["RAW_RX_MED_NAME", "RX_MED_NAME", "RAW_MED_NAME", "MEDICATION_NAME"],
            "rxnorm_candidates": ["RXNORM_CUI", "RAW_RXNORM_CUI", "RXNORM"],
            "route_candidates": ["RX_ROUTE", "RAW_RX_ROUTE", "ROUTE"],
            "enc_candidates": ["ENCOUNTERID", "ENCOUNTER_ID"],
        },
        "MED_ADMIN": {
            "view": "med_admin",
            "text_candidates": ["RAW_MEDADMIN_MED_NAME", "MEDADMIN_MED_NAME", "RAW_MED_NAME", "MEDICATION_NAME"],
            "rxnorm_candidates": ["MEDADMIN_CODE", "RAW_MEDADMIN_CODE", "RXNORM_CUI", "RAW_RXNORM_CUI"],
            "route_candidates": ["MEDADMIN_ROUTE", "RAW_MEDADMIN_ROUTE", "ROUTE"],
            "enc_candidates": ["ENCOUNTERID", "ENCOUNTER_ID"],
        },
    }

    encounter_views = []
    schema_rows = []
    for source, spec in source_specs.items():
        view = spec["view"]
        cols = set(con.execute(f"DESCRIBE {view}").fetchdf()["column_name"].astype(str))
        text_col = first_present(cols, spec["text_candidates"])
        rxnorm_col = first_present(cols, spec["rxnorm_candidates"])
        route_col = first_present(cols, spec["route_candidates"])
        enc_col = first_present(cols, spec["enc_candidates"])
        schema_rows.append({
            "source": source,
            "text_column": text_col or "<NONE>",
            "rxnorm_or_code_column": rxnorm_col or "<NONE>",
            "route_column": route_col or "<NONE>",
            "encounter_column": enc_col or "<NONE>",
        })
        text_expr = f"lower(coalesce(cast({qident(text_col)} as varchar),''))" if text_col else "''"
        route_expr = f"lower(trim(coalesce(cast({qident(route_col)} as varchar),'')))" if route_col else "''"
        rx_expr = f"trim(coalesce(cast({qident(rxnorm_col)} as varchar),''))" if rxnorm_col else "''"
        total = int(con.execute(f"SELECT count(*) FROM {view}").fetchone()[0])
        q = f"""
        SELECT
          sum(CASE WHEN {text_expr}<>'' THEN 1 ELSE 0 END)::BIGINT AS text_present,
          sum(CASE WHEN {rx_expr}<>'' THEN 1 ELSE 0 END)::BIGINT AS code_present,
          sum(CASE WHEN {route_expr}<>'' THEN 1 ELSE 0 END)::BIGINT AS route_present,
          sum(CASE WHEN regexp_matches({text_expr}, {qstr(BROAD_PATTERN)}) THEN 1 ELSE 0 END)::BIGINT AS broad_text_rows,
          sum(CASE WHEN regexp_matches({text_expr}, {qstr(BROAD_PATTERN)}) AND regexp_matches({route_expr}, {qstr(IV_ROUTE_PATTERN)}) THEN 1 ELSE 0 END)::BIGINT AS broad_text_iv_rows,
          sum(CASE WHEN regexp_matches({text_expr}, {qstr(BROAD_PATTERN)}) AND {route_expr}='' THEN 1 ELSE 0 END)::BIGINT AS broad_text_route_missing_rows
        FROM {view}
        """
        r = con.execute(q).fetchdf().iloc[0].to_dict()
        rows.append({"source": source, "total_rows": total, **{k: int(v) for k, v in r.items()}})

        for agent, pattern in BROAD_AGENTS.items():
            aq = f"""
            SELECT
              sum(CASE WHEN regexp_matches({text_expr}, {qstr(pattern)}) THEN 1 ELSE 0 END)::BIGINT AS rows,
              sum(CASE WHEN regexp_matches({text_expr}, {qstr(pattern)}) AND regexp_matches({route_expr}, {qstr(IV_ROUTE_PATTERN)}) THEN 1 ELSE 0 END)::BIGINT AS iv_rows
            FROM {view}
            """
            a = con.execute(aq).fetchone()
            agent_rows.append({"source": source, "agent": agent, "rows": safe_count(a[0]), "iv_rows": safe_count(a[1])})

        if enc_col:
            ev = f"broad_enc_{source.lower()}"
            con.execute(f"""
              CREATE TEMP TABLE {ev} AS
              SELECT DISTINCT cast({qident(enc_col)} as varchar) AS encounter_id
              FROM {view}
              WHERE {qident(enc_col)} IS NOT NULL
                AND regexp_matches({text_expr}, {qstr(BROAD_PATTERN)})
                AND regexp_matches({route_expr}, {qstr(IV_ROUTE_PATTERN)})
            """)
            encounter_views.append((source, ev))

    pd.DataFrame(schema_rows).to_csv(args.output_dir / "source_field_map.csv", index=False)
    pd.DataFrame(rows).to_csv(args.output_dir / "source_mapping_summary.csv", index=False)
    pd.DataFrame(agent_rows).to_csv(args.output_dir / "broad_agent_counts.csv", index=False)

    overlap = []
    if len(encounter_views) == 2:
        a_name, a_view = encounter_views[0]
        b_name, b_view = encounter_views[1]
        a_n = int(con.execute(f"SELECT count(*) FROM {a_view}").fetchone()[0])
        b_n = int(con.execute(f"SELECT count(*) FROM {b_view}").fetchone()[0])
        both = int(con.execute(f"SELECT count(*) FROM {a_view} a INNER JOIN {b_view} b USING(encounter_id)").fetchone()[0])
        overlap = [
            {"category": f"{a_name}_any_iv_broad_encounters", "count": safe_count(a_n)},
            {"category": f"{b_name}_any_iv_broad_encounters", "count": safe_count(b_n)},
            {"category": "both_sources_any_iv_broad_encounters", "count": safe_count(both)},
            {"category": f"{a_name}_only_any_iv_broad_encounters", "count": safe_count(a_n - both)},
            {"category": f"{b_name}_only_any_iv_broad_encounters", "count": safe_count(b_n - both)},
        ]
    pd.DataFrame(overlap).to_csv(args.output_dir / "encounter_source_overlap.csv", index=False)

    code = old_code_evidence(old_code)
    (args.output_dir / "old_antibiotic_code_evidence.json").write_text(json.dumps(code, indent=2), encoding="utf-8")

    summary = {
        "privacy_mode": "aggregate_only_no_ids_no_patient_rows_no_free_text_export",
        "minimum_reported_cell": MIN_CELL,
        "mimic_frozen_broad_pattern": BROAD_PATTERN,
        "mimic_primary_route_requirement": "IV/intravenous required",
        "prescribing_path": str(prescribing),
        "med_admin_path": str(med_admin),
        "old_antibiotic_code_path": str(old_code),
        "old_code_exists": old_code.exists(),
        "purpose": "Quantify mapping coverage and PRESCRIBING-versus-MED_ADMIN disagreement before freezing the PSU antibiotic source and broad-spectrum phenotype.",
        "guardrail": "No PSU source is selected as primary from this diagnostic alone; clinical/code mapping review follows aggregate results.",
    }
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
