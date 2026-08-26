#!/usr/bin/env python3
"""Aggregate-only audit of route ambiguity for PSU broad-spectrum antibiotics.

Focuses on agents where ingredient/brand RxNorm codes may span oral, inhaled, or IV
formulations despite route fields being predominantly UN. Patient identifiers and raw
medication names remain local; only predefined aggregate categories are exported.
"""
from __future__ import annotations

import argparse
import ast
import json
from pathlib import Path

import duckdb
import pandas as pd

MIN_CELL = 11
AMBIGUOUS = {
    "vancomycin": r"vancomycin|vancocin",
    "linezolid": r"linezolid|zyvox",
    "aztreonam": r"aztreonam|azactam|cayston",
}
ORAL_PAT = r"oral|capsule|tablet|solution|suspension|po\b|by mouth"
INJECT_PAT = r"inject|injection|intravenous|\biv\b|vial|premix|piggyback"
INHALED_PAT = r"inhal|nebul|cayston|respiratory"


def q(s: str) -> str:
    return "'" + s.replace("'", "''") + "'"


def suppress(df: pd.DataFrame, col: str = "count") -> pd.DataFrame:
    out = df.copy()
    if col in out.columns:
        mask = out[col].fillna(0) < MIN_CELL
        out.loc[mask, col] = pd.NA
        out["suppressed"] = mask
    return out


def parse_legacy(path: Path) -> tuple[set[str], set[str]]:
    tree = ast.parse(path.read_text(encoding="utf-8", errors="ignore"))
    includes: set[str] = set()
    excludes: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign):
            names = [t.id for t in node.targets if isinstance(t, ast.Name)]
            if not names:
                continue
            vals: set[str] = set()
            for sub in ast.walk(node.value):
                if isinstance(sub, ast.Constant) and isinstance(sub.value, str) and sub.value.isdigit():
                    vals.add(sub.value)
            if any("exclude" in n.lower() for n in names):
                excludes |= vals
            elif any("antibiotic" in n.lower() or "rxnorm" in n.lower() or "broad" in n.lower() for n in names):
                includes |= vals
    # Fallback to all digit strings in tuple/list assignments if parser names are not informative.
    if not includes:
        for node in ast.walk(tree):
            if isinstance(node, (ast.List, ast.Tuple, ast.Dict)):
                for sub in ast.walk(node):
                    if isinstance(sub, ast.Constant) and isinstance(sub.value, str) and sub.value.isdigit():
                        includes.add(sub.value)
    return includes - excludes, excludes


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("data_root", type=Path)
    ap.add_argument("--output-dir", type=Path, required=True)
    args = ap.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    root = args.data_root / "PCORnet"
    pfile = root / "parquet" / "prescribing.parquet"
    mfile = root / "parquet" / "med_admin.parquet"
    legacy = root / "code" / "config" / "codes_antibiotics.py"
    if not pfile.exists() or not mfile.exists() or not legacy.exists():
        raise FileNotFoundError("Required PSU files missing")

    include_codes, exclude_codes = parse_legacy(legacy)
    inc_sql = ",".join(q(x) for x in sorted(include_codes)) or "''"

    con = duckdb.connect()
    con.execute("PRAGMA threads=4")
    con.execute(f"CREATE VIEW p AS SELECT * FROM read_parquet({q(str(pfile))})")
    con.execute(f"CREATE VIEW m AS SELECT * FROM read_parquet({q(str(mfile))})")

    specs = [
        ("PRESCRIBING", "p", "raw_rx_med_name", "rxnorm_cui", "rx_route"),
        ("MED_ADMIN", "m", "raw_medadmin_med_name", "medadmin_code", "medadmin_route"),
    ]
    rows = []
    code_rows = []
    for source, view, text_col, code_col, route_col in specs:
        for agent, agent_pat in AMBIGUOUS.items():
            base = f"cast({code_col} AS VARCHAR) IN ({inc_sql}) AND regexp_matches(lower(coalesce(cast({text_col} AS VARCHAR),'')), {q(agent_pat)})"
            df = con.execute(f"""
                SELECT {q(source)} AS source, {q(agent)} AS agent,
                       CASE
                         WHEN regexp_matches(lower(coalesce(cast({text_col} AS VARCHAR),'')), {q(ORAL_PAT)}) THEN 'oral_name_signal'
                         WHEN regexp_matches(lower(coalesce(cast({text_col} AS VARCHAR),'')), {q(INHALED_PAT)}) THEN 'inhaled_name_signal'
                         WHEN regexp_matches(lower(coalesce(cast({text_col} AS VARCHAR),'')), {q(INJECT_PAT)}) THEN 'injectable_name_signal'
                         ELSE 'no_formulation_name_signal'
                       END AS formulation_signal,
                       coalesce(cast({route_col} AS VARCHAR), '<MISSING>') AS route,
                       count(*)::BIGINT AS count
                FROM {view}
                WHERE {base}
                GROUP BY 1,2,3,4
                ORDER BY 2, count DESC
            """).fetchdf()
            rows.append(suppress(df))
            cdf = con.execute(f"""
                SELECT {q(source)} AS source, {q(agent)} AS agent,
                       cast({code_col} AS VARCHAR) AS code,
                       count(*)::BIGINT AS count,
                       sum(regexp_matches(lower(coalesce(cast({text_col} AS VARCHAR),'')), {q(ORAL_PAT)}))::BIGINT AS oral_name_rows,
                       sum(regexp_matches(lower(coalesce(cast({text_col} AS VARCHAR),'')), {q(INHALED_PAT)}))::BIGINT AS inhaled_name_rows,
                       sum(regexp_matches(lower(coalesce(cast({text_col} AS VARCHAR),'')), {q(INJECT_PAT)}))::BIGINT AS injectable_name_rows
                FROM {view}
                WHERE {base}
                GROUP BY 1,2,3 ORDER BY count DESC
            """).fetchdf()
            code_rows.append(suppress(cdf))

    pd.concat(rows, ignore_index=True).to_csv(args.output_dir / "ambiguous_agent_formulation_counts.csv", index=False)
    pd.concat(code_rows, ignore_index=True).to_csv(args.output_dir / "ambiguous_agent_code_counts.csv", index=False)

    summary = {
        "privacy_mode": "aggregate_only_no_ids_no_patient_rows_no_raw_names_export",
        "minimum_reported_cell": MIN_CELL,
        "corrected_legacy_inclusion_code_count": len(include_codes),
        "explicit_exclusion_code_count": len(exclude_codes),
        "agents_audited": list(AMBIGUOUS),
        "purpose": "Check whether corrected legacy RxNorm broad-spectrum mapping still admits oral or inhaled formulations for route-ambiguous agents before freezing PSU antibiotic phenotype.",
        "guardrail": "Diagnostic only; raw medication names are used locally only for predefined formulation-pattern classification and are never exported.",
    }
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
