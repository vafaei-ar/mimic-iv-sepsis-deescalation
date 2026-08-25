#!/usr/bin/env python3
import argparse
import ast
import csv
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

KEYWORDS = [
    "sepsis", "icu", "intime", "admit_time", "admit_date", "discharge_time",
    "encounter", "med_admin", "prescribing", "rxnorm", "antibiotic", "broad",
    "culture", "micro", "specimen", "result_time", "result_date", "loinc",
    "sofa", "vasopressor", "death", "mortality", "lab_result_cm", "midnight",
]
SKIP_PARTS = {"site-packages", ".venv", "venv", "env", "__pycache__", ".git", "los"}


def safe_rel(path: Path, root: Path) -> str:
    try:
        return str(path.relative_to(root))
    except Exception:
        return path.name


def should_skip(path: Path) -> bool:
    return any(part in SKIP_PARTS for part in path.parts)


def ast_summary(text: str):
    out = {"functions": [], "classes": [], "imports": [], "assignments": []}
    try:
        tree = ast.parse(text)
    except Exception:
        return out
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            out["functions"].append(node.name)
        elif isinstance(node, ast.ClassDef):
            out["classes"].append(node.name)
        elif isinstance(node, ast.Import):
            out["imports"].extend(a.name for a in node.names)
        elif isinstance(node, ast.ImportFrom):
            mod = node.module or ""
            out["imports"].append(mod)
        elif isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name):
                    out["assignments"].append(target.id)
    for k in out:
        out[k] = sorted(set(out[k]))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("data_root")
    ap.add_argument("--output-dir", required=True)
    args = ap.parse_args()

    root = Path(args.data_root).resolve()
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    py_files = []
    for p in root.rglob("*.py"):
        if should_skip(p):
            continue
        # prioritize project analysis code, but still allow root-level Python helpers
        py_files.append(p)
    py_files = sorted(py_files)

    inventory = []
    evidence = []
    keyword_files = defaultdict(set)
    keyword_counts = Counter()
    pipeline_signals = defaultdict(list)

    for p in py_files:
        rel = safe_rel(p, root)
        try:
            text = p.read_text(encoding="utf-8", errors="replace")
        except Exception as e:
            inventory.append({"file": rel, "read_status": f"error:{type(e).__name__}", "n_lines": None,
                              "functions": "", "classes": "", "imports": "", "assignments": ""})
            continue

        lines = text.splitlines()
        s = ast_summary(text)
        inventory.append({
            "file": rel,
            "read_status": "ok",
            "n_lines": len(lines),
            "functions": ";".join(s["functions"][:100]),
            "classes": ";".join(s["classes"][:50]),
            "imports": ";".join(s["imports"][:100]),
            "assignments": ";".join(s["assignments"][:150]),
        })

        lower_lines = [ln.lower() for ln in lines]
        matched_lines = set()
        for kw in KEYWORDS:
            for i, low in enumerate(lower_lines):
                if kw in low:
                    keyword_counts[kw] += 1
                    keyword_files[kw].add(rel)
                    if i not in matched_lines:
                        # sanitized code-only evidence, no execution/data values
                        code = lines[i].strip()
                        if len(code) > 500:
                            code = code[:500]
                        evidence.append({"file": rel, "line": i + 1, "keyword": kw, "code": code})
                        matched_lines.add(i)

        # high-level pipeline signals inferred only from code text
        tests = {
            "uses_sepsis_encounter": r"sepsis_encounter",
            "uses_med_admin": r"med_admin",
            "uses_prescribing": r"prescribing",
            "uses_lab_result_cm": r"lab_result_cm|lab_reduced",
            "uses_rxnorm": r"rxnorm",
            "uses_sofa": r"\bsofa\b",
            "uses_icu_term": r"\bicu\b|intensive care",
            "uses_micro_term": r"culture|microbiology|specimen",
            "uses_midnight_literal": r"midnight|00:00|time\(0\s*,\s*0",
            "uses_admit_time": r"admit_time",
            "uses_discharge_time": r"discharge_time",
            "uses_result_time": r"result_time",
            "uses_vasopressor": r"vasopressor|norepinephrine|epinephrine|vasopressin|phenylephrine",
        }
        lowtxt = text.lower()
        for name, pat in tests.items():
            if re.search(pat, lowtxt, flags=re.I):
                pipeline_signals[name].append(rel)

    with (outdir / "python_inventory.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["file", "read_status", "n_lines", "functions", "classes", "imports", "assignments"])
        w.writeheader(); w.writerows(inventory)

    # cap evidence to avoid oversized artifact while preserving all important files/keywords
    evidence = evidence[:5000]
    with (outdir / "relevant_code_evidence.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["file", "line", "keyword", "code"])
        w.writeheader(); w.writerows(evidence)

    pipeline = {k: sorted(set(v)) for k, v in pipeline_signals.items()}
    (outdir / "inferred_pipeline_map.json").write_text(json.dumps(pipeline, indent=2), encoding="utf-8")

    summary = {
        "privacy_mode": "source_code_only_no_patient_rows_no_identifiers_no_data_values",
        "data_root": str(root),
        "python_files_scanned": len(py_files),
        "python_files_readable": sum(1 for r in inventory if r["read_status"] == "ok"),
        "keyword_counts": dict(keyword_counts),
        "keyword_file_counts": {k: len(v) for k, v in keyword_files.items()},
        "important_signals": pipeline,
        "note": "This audit reads Python source only. It does not read patient rows or export patient-level data.",
    }
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps({"status": "ok", "python_files_scanned": len(py_files), "output_dir": str(outdir)}, indent=2))

if __name__ == "__main__":
    main()
