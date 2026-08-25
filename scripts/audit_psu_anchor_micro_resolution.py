#!/usr/bin/env python3
import argparse
import csv
import json
import re
from collections import Counter
from pathlib import Path

import pyarrow.compute as pc
import pyarrow.parquet as pq

MIN_CELL = 11


def safe_count(n: int):
    return n if n >= MIN_CELL else None


def write_csv(path: Path, fieldnames, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)


def candidate_paths(root: Path):
    terms = ("icu", "intensive", "critical", "adt", "transfer", "location", "encounter", "micro", "culture", "organism")
    rows = []
    for p in root.rglob("*"):
        rel = str(p.relative_to(root))
        low = rel.lower()
        if not any(t in low for t in terms):
            continue
        try:
            is_file = p.is_file()
            size = p.stat().st_size if is_file else None
        except OSError:
            is_file = False
            size = None
        rows.append({"relative_path": rel, "kind": "file" if is_file else "directory", "size_bytes": size})
        if len(rows) >= 500:
            break
    return rows


def scan_micro(lab_path: Path):
    keywords = {
        "culture": r"culture|cult\b",
        "microbiology": r"microbiol",
        "organism": r"organism",
        "susceptibility": r"suscept|sensitivity",
        "gram_stain": r"gram\s*stain|gramstain",
        "blood": r"blood",
        "urine": r"urine",
        "respiratory": r"respir|sputum|bronch",
        "wound": r"wound",
        "fungal": r"fung|yeast|mold",
        "afb": r"acid[- ]?fast|\bafb\b|mycobacter",
        "body_fluid": r"body\s*fluid|sterile\s*site|pleural|peritoneal|csf",
    }
    keyword_counts = Counter()
    loinc_counts = Counter()
    specimen_counts = Counter()
    result_qual_counts = Counter()
    total_candidates = 0
    pf = pq.ParquetFile(lab_path)
    cols = ["RAW_LAB_NAME", "LAB_LOINC", "SPECIMEN_SOURCE", "RESULT_QUAL"]
    for rg in range(pf.num_row_groups):
        tbl = pf.read_row_group(rg, columns=cols)
        names = pc.fill_null(tbl["RAW_LAB_NAME"], "")
        names_lower = pc.utf8_lower(names)
        any_mask = None
        for label, pattern in keywords.items():
            mask = pc.match_substring_regex(names_lower, pattern)
            n = int(pc.sum(pc.cast(mask, "int64")).as_py() or 0)
            keyword_counts[label] += n
            any_mask = mask if any_mask is None else pc.or_(any_mask, mask)
        n_any = int(pc.sum(pc.cast(any_mask, "int64")).as_py() or 0)
        total_candidates += n_any
        if n_any:
            sub = tbl.filter(any_mask)
            for value in sub["LAB_LOINC"].to_pylist():
                if value not in (None, ""):
                    loinc_counts[str(value)] += 1
            for value in sub["SPECIMEN_SOURCE"].to_pylist():
                if value not in (None, ""):
                    specimen_counts[str(value)] += 1
            for value in sub["RESULT_QUAL"].to_pylist():
                if value not in (None, ""):
                    result_qual_counts[str(value)] += 1
    return {
        "total_candidate_rows": total_candidates,
        "keyword_counts": keyword_counts,
        "loinc_counts": loinc_counts,
        "specimen_counts": specimen_counts,
        "result_qual_counts": result_qual_counts,
    }


def scan_obs_icu(obs_path: Path):
    pattern = r"icu|intensive|critical\s*care|transfer|patient\s*location|care\s*unit|unit\s*location|adt"
    pf = pq.ParquetFile(obs_path)
    counts = Counter()
    matched_rows = 0
    cols = ["obsclin_type", "obsclin_code", "raw_obsclin_name", "raw_obsclin_code", "raw_obsclin_type"]
    for rg in range(pf.num_row_groups):
        tbl = pf.read_row_group(rg, columns=cols)
        any_mask = None
        for col in cols:
            arr = pc.utf8_lower(pc.fill_null(tbl[col], ""))
            mask = pc.match_substring_regex(arr, pattern)
            any_mask = mask if any_mask is None else pc.or_(any_mask, mask)
        n_any = int(pc.sum(pc.cast(any_mask, "int64")).as_py() or 0)
        matched_rows += n_any
        if n_any:
            sub = tbl.filter(any_mask)
            for col in ("obsclin_type", "obsclin_code", "raw_obsclin_type"):
                for value in sub[col].to_pylist():
                    if value not in (None, ""):
                        counts[(col, str(value))] += 1
    return matched_rows, counts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("root")
    ap.add_argument("--output-dir", required=True)
    args = ap.parse_args()
    root = Path(args.root)
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)

    path_rows = candidate_paths(root)
    write_csv(out / "candidate_path_inventory.csv", ["relative_path", "kind", "size_bytes"], path_rows)

    lab = root / "PCORnet/parquet/lab/lab_reduced.parquet"
    micro = scan_micro(lab)
    micro_rows = []
    for label, n in sorted(micro["keyword_counts"].items()):
        micro_rows.append({"category": "keyword", "value": label, "count": safe_count(n), "suppressed": n < MIN_CELL})
    for value, n in micro["loinc_counts"].most_common(100):
        if n >= MIN_CELL:
            micro_rows.append({"category": "loinc", "value": value, "count": n, "suppressed": False})
    for value, n in micro["specimen_counts"].most_common(50):
        if n >= MIN_CELL:
            micro_rows.append({"category": "specimen_source", "value": value, "count": n, "suppressed": False})
    for value, n in micro["result_qual_counts"].most_common(50):
        if n >= MIN_CELL:
            micro_rows.append({"category": "result_qual", "value": value, "count": n, "suppressed": False})
    write_csv(out / "microbiology_aggregate_counts.csv", ["category", "value", "count", "suppressed"], micro_rows)

    obs = root / "PCORnet/parquet/obs_clin.parquet"
    icu_matched, icu_counts = scan_obs_icu(obs)
    icu_rows = []
    for (col, value), n in icu_counts.most_common(100):
        if n >= MIN_CELL:
            icu_rows.append({"column": col, "value": value, "count": n, "suppressed": False})
    write_csv(out / "icu_observation_aggregate_counts.csv", ["column", "value", "count", "suppressed"], icu_rows)

    summary = {
        "privacy_mode": "aggregate_only_no_ids_no_rows_no_result_free_text_export",
        "minimum_reported_cell": MIN_CELL,
        "candidate_path_count": len(path_rows),
        "microbiology_candidate_rows": safe_count(micro["total_candidate_rows"]),
        "microbiology_candidate_rows_suppressed": micro["total_candidate_rows"] < MIN_CELL,
        "icu_observation_candidate_rows": safe_count(icu_matched),
        "icu_observation_candidate_rows_suppressed": icu_matched < MIN_CELL,
        "assessment": {
            "icu_anchor": "Exact ICU intime is supportable only if a local ADT/location source or a validated ICU observation signal is identified; hospital admit_date alone is insufficient for strict MIMIC replication.",
            "microbiology": "Culture/result-availability phenotype can be frozen only after aggregate lab-name/LOINC signals identify clinical culture tests and positivity semantics are validated.",
        },
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": "ok", "output_dir": str(out), "candidate_paths": len(path_rows), "micro_candidates": micro["total_candidate_rows"], "icu_obs_candidates": icu_matched}, indent=2))


if __name__ == "__main__":
    main()
