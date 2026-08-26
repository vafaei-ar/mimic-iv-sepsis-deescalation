#!/usr/bin/env python3
"""Aggregate-only audit of legacy PSU RxNorm specificity and source concordance.

This corrected version separates inclusion and exclusion RxNorm codes from the
legacy mapping instead of regex-collecting every numeric literal in the file.
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
BROAD_PATTERN = (
    "vancomycin|piperacillin|tazobactam|zosyn|cefepime|ceftazidime|meropenem|"
    "imipenem|aztreonam|linezolid|daptomycin|ceftolozane|avibactam"
)
AGENT_PATTERNS = {
    "vancomycin": "vancomycin",
    "piperacillin_tazobactam": "piperacillin|tazobactam|zosyn",
    "cefepime": "cefepime",
    "ceftazidime": "ceftazidime",
    "meropenem": "meropenem",
    "imipenem": "imipenem",
    "aztreonam": "aztreonam",
    "linezolid": "linezolid",
    "daptomycin": "daptomycin",
    "ceftolozane": "ceftolozane",
    "avibactam": "avibactam",
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


def parse_assignments(path: Path) -> dict[str, object]:
    tree = ast.parse(path.read_text(encoding="utf-8", errors="ignore"))
    out: dict[str, object] = {}
    for node in tree.body:
        if isinstance(node, ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name):
            try:
                out[node.targets[0].id] = ast.literal_eval(node.value)
            except Exception:
                pass
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
    old_code = root / "code" / "config" / "codes_antibiotics.py"
    if not prescribing.exists() or not med_admin.exists() or not old_code.exists():
        raise FileNotFoundError("Required PSU prescribing, med_admin, or legacy antibiotic code missing")

    old_text = old_code.read_text(encoding="utf-8", errors="ignore")
    assignments = parse_assignments(old_code)

    include_codes: set[str] = set()
    exclude_codes: set[str] = set()
    object_rows = []
    for name, obj in assignments.items():
        nums = collect_numeric_strings(obj)
        if not nums:
            continue
        lname = name.lower()
        if "exclude" in lname:
            role = "exclude"
            exclude_codes |= nums
        elif any(k in lname for k in ("antibi", "broad", "rxnorm", "drug", "code")):
            role = "include"
            include_codes |= nums
        else:
            role = "other"
        object_rows.append({"variable": name, "role": role, "numeric_code_count": len(nums)})

    # Explicit exclusions override inclusion. This is critical for oral vancomycin.
    rxnorm_codes = sorted(include_codes - exclude_codes)
    excluded_codes = sorted(exclude_codes)
    rxnorm_sql = ",".join(q(x) for x in rxnorm_codes) or "''"
    excluded_sql = ",".join(q(x) for x in excluded_codes) or "''"

    pd.DataFrame(object_rows).to_csv(args.output_dir / "legacy_mapping_objects.csv", index=False)

    # Export source-code context for parsed include/exclude codes only; no patient data.
    lines = old_text.splitlines()
    context_rows = []
    for code in sorted(set(rxnorm_codes) | set(excluded_codes)):
        hits = [i for i, line in enumerate(lines) if re.search(rf"(?<!\d){re.escape(code)}(?!\d)", line)]
        for i in hits[:3]:
            start, stop = max(0, i - 1), min(len(lines), i + 2)
            context = " | ".join(x.strip() for x in lines[start:stop] if x.strip())
            context_rows.append({
                "code": code,
                "mapping_role": "exclude" if code in exclude_codes else "include",
                "legacy_code_context": context[:500],
            })
    pd.DataFrame(context_rows).drop_duplicates().to_csv(args.output_dir / "legacy_rxnorm_code_context.csv", index=False)

    con = duckdb.connect()
    con.execute("PRAGMA threads=4")
    con.execute(f"CREATE VIEW p AS SELECT * FROM read_parquet({q(str(prescribing))})")
    con.execute(f"CREATE VIEW m AS SELECT * FROM read_parquet({q(str(med_admin))})")

    specs = [
        ("PRESCRIBING", "p", "raw_rx_med_name", "rxnorm_cui", "rx_route", "encounterid"),
        ("MED_ADMIN", "m", "raw_medadmin_med_name", "medadmin_code", "medadmin_route", "encounterid"),
    ]

    route_frames = []
    agreement_rows = []
    code_agent_frames = []
    for source, view, text_col, code_col, route_col, enc_col in specs:
        route = con.execute(f"""
            SELECT {q(source)} AS source,
                   coalesce(cast({route_col} AS VARCHAR), '<MISSING>') AS route,
                   count(*)::BIGINT AS count
            FROM {view}
            WHERE cast({code_col} AS VARCHAR) IN ({rxnorm_sql})
            GROUP BY 1,2 ORDER BY count DESC
        """).fetchdf()
        route_frames.append(suppress(route))

        agreement = con.execute(f"""
            SELECT {q(source)} AS source,
              count(*) FILTER (WHERE regexp_matches(lower(coalesce(cast({text_col} AS VARCHAR),'')), {q(BROAD_PATTERN)})
                               AND cast({code_col} AS VARCHAR) IN ({rxnorm_sql}))::BIGINT AS text_and_legacy_rows,
              count(*) FILTER (WHERE regexp_matches(lower(coalesce(cast({text_col} AS VARCHAR),'')), {q(BROAD_PATTERN)})
                               AND NOT cast({code_col} AS VARCHAR) IN ({rxnorm_sql}))::BIGINT AS text_only_rows,
              count(*) FILTER (WHERE NOT regexp_matches(lower(coalesce(cast({text_col} AS VARCHAR),'')), {q(BROAD_PATTERN)})
                               AND cast({code_col} AS VARCHAR) IN ({rxnorm_sql}))::BIGINT AS legacy_only_rows,
              count(*) FILTER (WHERE cast({code_col} AS VARCHAR) IN ({excluded_sql}))::BIGINT AS explicit_excluded_rows,
              count(DISTINCT CASE WHEN cast({code_col} AS VARCHAR) IN ({rxnorm_sql}) THEN {enc_col} END)::BIGINT AS legacy_encounters
            FROM {view}
        """).fetchdf().iloc[0].to_dict()
        agreement_rows.append(agreement)

        parts = []
        for agent, pat in AGENT_PATTERNS.items():
            parts.append(
                f"SELECT {q(source)} AS source, cast({code_col} AS VARCHAR) AS code, {q(agent)} AS agent, "
                f"count(*)::BIGINT AS count FROM {view} WHERE cast({code_col} AS VARCHAR) IN ({rxnorm_sql}) "
                f"AND regexp_matches(lower(coalesce(cast({text_col} AS VARCHAR),'')), {q(pat)}) GROUP BY 1,2,3"
            )
        code_agent = con.execute(" UNION ALL ".join(parts)).fetchdf()
        code_agent_frames.append(suppress(code_agent))

    pd.concat(route_frames, ignore_index=True).to_csv(args.output_dir / "legacy_rxnorm_route_counts.csv", index=False)
    pd.DataFrame(agreement_rows).to_csv(args.output_dir / "text_rxnorm_agreement.csv", index=False)
    pd.concat(code_agent_frames, ignore_index=True).to_csv(args.output_dir / "legacy_code_agent_counts.csv", index=False)

    overlap = con.execute(f"""
        WITH pe AS (
          SELECT DISTINCT cast(encounterid AS VARCHAR) AS encounterid
          FROM p WHERE cast(rxnorm_cui AS VARCHAR) IN ({rxnorm_sql}) AND encounterid IS NOT NULL
        ),
        me AS (
          SELECT DISTINCT cast(encounterid AS VARCHAR) AS encounterid
          FROM m WHERE cast(medadmin_code AS VARCHAR) IN ({rxnorm_sql}) AND encounterid IS NOT NULL
        ),
        u AS (SELECT encounterid FROM pe UNION SELECT encounterid FROM me)
        SELECT
          count(*) FILTER (WHERE pe.encounterid IS NOT NULL)::BIGINT AS prescribing_legacy_encounters,
          count(*) FILTER (WHERE me.encounterid IS NOT NULL)::BIGINT AS med_admin_legacy_encounters,
          count(*) FILTER (WHERE pe.encounterid IS NOT NULL AND me.encounterid IS NOT NULL)::BIGINT AS both_legacy_encounters,
          count(*) FILTER (WHERE pe.encounterid IS NOT NULL AND me.encounterid IS NULL)::BIGINT AS prescribing_only_legacy_encounters,
          count(*) FILTER (WHERE pe.encounterid IS NULL AND me.encounterid IS NOT NULL)::BIGINT AS med_admin_only_legacy_encounters
        FROM u LEFT JOIN pe USING (encounterid) LEFT JOIN me USING (encounterid)
    """).fetchdf()
    overlap.to_csv(args.output_dir / "legacy_encounter_source_overlap.csv", index=False)

    summary = {
        "privacy_mode": "aggregate_only_no_ids_no_patient_rows_no_free_text_export",
        "minimum_reported_cell": MIN_CELL,
        "parsed_assignment_count": len(assignments),
        "parsed_include_code_count_before_exclusions": len(include_codes),
        "parsed_explicit_exclude_code_count": len(exclude_codes),
        "final_legacy_rxnorm_code_count": len(rxnorm_codes),
        "correction": "Prior audit regex-collected all numeric literals, including explicit oral vancomycin exclusion codes. This run parses literal assignment objects and removes explicit exclusions before assessing specificity/concordance.",
        "guardrail": "Diagnostic only; do not freeze source or IV phenotype unless corrected mapping and timing concordance support it.",
    }
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
