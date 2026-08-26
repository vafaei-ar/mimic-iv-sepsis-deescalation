#!/usr/bin/env python3
"""Parse prior PSU antibiotic mapping code without conflating include/exclude RxNorm codes."""
from __future__ import annotations

import argparse
import ast
import json
import re
from pathlib import Path

import duckdb
import pandas as pd

MIN_CELL = 11
BROAD_PATTERN = (
    "vancomycin|piperacillin|tazobactam|zosyn|cefepime|ceftazidime|meropenem|"
    "imipenem|aztreonam|linezolid|daptomycin|ceftolozane|avibactam"
)


def q(s: str) -> str:
    return "'" + s.replace("'", "''") + "'"


def literal_eval_assignments(path: Path) -> dict[str, object]:
    tree = ast.parse(path.read_text(encoding="utf-8", errors="ignore"))
    out: dict[str, object] = {}
    for node in tree.body:
        if isinstance(node, ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name):
            try:
                out[node.targets[0].id] = ast.literal_eval(node.value)
            except Exception:
                pass
    return out


def collect_numeric_strings(obj: object) -> set[str]:
    vals: set[str] = set()
    if isinstance(obj, str) and re.fullmatch(r"\d{3,9}", obj):
        vals.add(obj)
    elif isinstance(obj, dict):
        for k, v in obj.items():
            vals |= collect_numeric_strings(k)
            vals |= collect_numeric_strings(v)
    elif isinstance(obj, (list, tuple, set)):
        for x in obj:
            vals |= collect_numeric_strings(x)
    return vals


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("data_root", type=Path)
    ap.add_argument("--output-dir", type=Path, required=True)
    args = ap.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    root = args.data_root / "PCORnet"
    code = root / "code" / "config" / "codes_antibiotics.py"
    pfile = root / "parquet" / "prescribing.parquet"
    mfile = root / "parquet" / "med_admin.parquet"
    if not code.exists():
        raise FileNotFoundError(code)

    assigns = literal_eval_assignments(code)
    rows = []
    include_codes: set[str] = set()
    exclude_codes: set[str] = set()
    for name, obj in assigns.items():
        nums = collect_numeric_strings(obj)
        if not nums:
            continue
        lname = name.lower()
        role = "exclude" if "exclude" in lname else ("include" if any(k in lname for k in ["antibi", "broad", "rxnorm", "code", "drug"]) else "other")
        if role == "include":
            include_codes |= nums
        elif role == "exclude":
            exclude_codes |= nums
        rows.append({"variable": name, "role": role, "numeric_code_count": len(nums)})

    # Conservative parsed include set: remove any explicitly excluded codes.
    parsed_include = sorted(include_codes - exclude_codes)
    parsed_exclude = sorted(exclude_codes)
    pd.DataFrame(rows).sort_values(["role", "variable"]).to_csv(args.output_dir / "mapping_object_summary.csv", index=False)
    pd.DataFrame({"code": parsed_include}).to_csv(args.output_dir / "parsed_include_codes.csv", index=False)
    pd.DataFrame({"code": parsed_exclude}).to_csv(args.output_dir / "parsed_exclude_codes.csv", index=False)

    con = duckdb.connect()
    con.execute(f"CREATE VIEW p AS SELECT * FROM read_parquet({q(str(pfile))})")
    con.execute(f"CREATE VIEW m AS SELECT * FROM read_parquet({q(str(mfile))})")
    inc_sql = ",".join(q(x) for x in parsed_include) or "''"
    exc_sql = ",".join(q(x) for x in parsed_exclude) or "''"
    specs = [
        ("PRESCRIBING", "p", "raw_rx_med_name", "rxnorm_cui", "rx_route", "encounterid"),
        ("MED_ADMIN", "m", "raw_medadmin_med_name", "medadmin_code", "medadmin_route", "encounterid"),
    ]
    out = []
    for source, view, text_col, code_col, route_col, enc_col in specs:
        df = con.execute(f"""
        SELECT {q(source)} AS source,
          sum(cast({code_col} AS VARCHAR) IN ({inc_sql}))::BIGINT AS parsed_include_rows,
          sum(cast({code_col} AS VARCHAR) IN ({exc_sql}))::BIGINT AS explicit_exclude_rows,
          sum(cast({code_col} AS VARCHAR) IN ({inc_sql}) AND regexp_matches(lower(coalesce(cast({text_col} AS VARCHAR),'')), {q(BROAD_PATTERN)}))::BIGINT AS include_and_broad_text_rows,
          count(DISTINCT CASE WHEN cast({code_col} AS VARCHAR) IN ({inc_sql}) THEN {enc_col} END)::BIGINT AS parsed_include_encounters,
          sum(cast({code_col} AS VARCHAR) IN ({inc_sql}) AND upper(coalesce(cast({route_col} AS VARCHAR),''))='ORAL')::BIGINT AS include_oral_route_rows,
          sum(cast({code_col} AS VARCHAR) IN ({inc_sql}) AND upper(coalesce(cast({route_col} AS VARCHAR),'')) IN ('INTRAVENOUS','IV'))::BIGINT AS include_iv_route_rows
        FROM {view}
        """).fetchdf()
        out.append(df)
    pd.concat(out, ignore_index=True).to_csv(args.output_dir / "parsed_mapping_coverage.csv", index=False)

    summary = {
        "privacy_mode": "aggregate_only_source_code_config_and_counts_no_ids_no_patient_rows",
        "minimum_reported_cell": MIN_CELL,
        "parsed_assignment_count": len(assigns),
        "parsed_include_code_count": len(parsed_include),
        "parsed_exclude_code_count": len(parsed_exclude),
        "important_guardrail": "Previous regex extraction mixed numeric literals from include and exclude sections. This audit separates them using AST-assigned objects and removes explicit exclusions before coverage assessment.",
    }
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

if __name__ == "__main__":
    main()
