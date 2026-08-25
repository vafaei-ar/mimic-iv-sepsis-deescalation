#!/usr/bin/env python
from __future__ import annotations

import argparse
import ast
import json
import re
from collections import Counter
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq

EXPECTED_ROOT = Path("/home/asadr/Depts/PHS/PATH_CDM/Hwang_Bonavia").resolve()
MIN_CELL = 11

TABLES = {
    "sepsis_encounter": "PCORnet/parquet/sepsis_encounter.parquet",
    "prescribing": "PCORnet/parquet/prescribing.parquet",
    "med_admin": "PCORnet/parquet/med_admin.parquet",
    "lab_reduced": "PCORnet/parquet/lab/lab_reduced.parquet",
    "death": "PCORnet/parquet/death.parquet",
    "obs_clin": "PCORnet/parquet/obs_clin.parquet",
    "sepsis_vital": "PCORnet/parquet/sepsis_vital.parquet",
}

CODE_FILES = {
    "cohort_extractor": "PCORnet/code/extract/01_cohort.py",
    "antibiotic_codes": "PCORnet/code/config/codes_antibiotics.py",
    "settings": "PCORnet/code/config/settings.py",
}

MICRO_NAME_RE = re.compile(
    r"(?:blood|urine|respiratory|sputum|wound|sterile|body fluid|csf|cerebrospinal|bronch|stool|fungal|anaerobic|aerobic)?\\s*"
    r"(?:culture|cultures)|microbiol|organism|susceptib|gram stain",
    re.IGNORECASE,
)
POS_RE = re.compile(r"\\b(?:positive|detected|growth|grew|isolated)\\b", re.IGNORECASE)
NEG_RE = re.compile(r"\\b(?:negative|no growth|not detected|none isolated)\\b", re.IGNORECASE)


def safe_root(root: Path) -> Path:
    root = root.expanduser().resolve()
    if root != EXPECTED_ROOT:
        raise ValueError(f"Refusing unexpected data root: {root}")
    return root


def null_count_from_metadata(path: Path, column: str) -> tuple[int | None, int]:
    pf = pq.ParquetFile(path)
    names = pf.schema_arrow.names
    if column not in names:
        return None, int(pf.metadata.num_rows)
    idx = names.index(column)
    total_null = 0
    for rg in range(pf.metadata.num_row_groups):
        stats = pf.metadata.row_group(rg).column(idx).statistics
        if stats is None or stats.null_count is None:
            return None, int(pf.metadata.num_rows)
        total_null += int(stats.null_count)
    return total_null, int(pf.metadata.num_rows)


def timing_row(root: Path, table: str, column: str, role: str) -> dict:
    path = root / TABLES[table]
    pf = pq.ParquetFile(path)
    present = column in pf.schema_arrow.names
    null_count, n = null_count_from_metadata(path, column) if present else (None, int(pf.metadata.num_rows))
    return {
        "domain": role,
        "table": table,
        "column": column,
        "present": present,
        "row_count": n,
        "null_count_metadata": null_count,
        "nonmissing_fraction_metadata": (1 - null_count / n) if null_count is not None and n else None,
    }


def aggregate_categories(path: Path, columns: list[str], batch_size: int = 250_000) -> list[dict]:
    pf = pq.ParquetFile(path)
    present = [c for c in columns if c in pf.schema_arrow.names]
    counters = {c: Counter() for c in present}
    missing = Counter()
    totals = Counter()
    for batch in pf.iter_batches(columns=present, batch_size=batch_size):
        d = batch.to_pandas()
        for c in present:
            s = d[c]
            totals[c] += len(s)
            missing[c] += int(s.isna().sum())
            vals = s.dropna().astype(str).str.strip()
            vals = vals[vals != ""]
            counters[c].update(vals.tolist())
    rows: list[dict] = []
    for c in present:
        rows.append({"column": c, "value": "<MISSING>", "count": int(missing[c]), "suppressed": missing[c] < MIN_CELL})
        for value, count in counters[c].most_common(25):
            rows.append({
                "column": c,
                "value": value if count >= MIN_CELL else "<SUPPRESSED>",
                "count": int(count) if count >= MIN_CELL else None,
                "suppressed": count < MIN_CELL,
            })
    return rows


def scan_microbiology(path: Path, batch_size: int = 150_000) -> dict:
    pf = pq.ParquetFile(path)
    needed = [
        "RAW_LAB_NAME", "RAW_PANEL", "LAB_LOINC", "SPECIMEN_SOURCE",
        "SPECIMEN_DATE", "SPECIMEN_TIME", "RESULT_DATE", "RESULT_TIME",
        "RESULT_QUAL", "RESULT_SNOMED", "RAW_RESULT",
    ]
    present = [c for c in needed if c in pf.schema_arrow.names]
    counts = Counter()
    loincs = Counter()
    specimens = Counter()
    for batch in pf.iter_batches(columns=present, batch_size=batch_size):
        d = batch.to_pandas()
        name = d.get("RAW_LAB_NAME", pd.Series("", index=d.index)).fillna("").astype(str)
        panel = d.get("RAW_PANEL", pd.Series("", index=d.index)).fillna("").astype(str)
        mask = (name + " " + panel).str.contains(MICRO_NAME_RE, regex=True, na=False)
        if not bool(mask.any()):
            continue
        m = d.loc[mask]
        counts["candidate_micro_rows"] += len(m)
        for c in ["SPECIMEN_DATE", "SPECIMEN_TIME", "RESULT_DATE", "RESULT_TIME", "RESULT_QUAL", "RESULT_SNOMED", "RAW_RESULT"]:
            if c in m.columns:
                counts[f"nonmissing_{c.lower()}"] += int(m[c].notna().sum())
        if "RAW_RESULT" in m.columns:
            raw = m["RAW_RESULT"].fillna("").astype(str)
            counts["raw_result_positive_term"] += int(raw.str.contains(POS_RE, regex=True, na=False).sum())
            counts["raw_result_negative_term"] += int(raw.str.contains(NEG_RE, regex=True, na=False).sum())
        if "RESULT_QUAL" in m.columns:
            qual = m["RESULT_QUAL"].fillna("").astype(str)
            counts["result_qual_positive_term"] += int(qual.str.contains(POS_RE, regex=True, na=False).sum())
            counts["result_qual_negative_term"] += int(qual.str.contains(NEG_RE, regex=True, na=False).sum())
        if "LAB_LOINC" in m.columns:
            loincs.update(m["LAB_LOINC"].dropna().astype(str).str.strip().tolist())
        if "SPECIMEN_SOURCE" in m.columns:
            specimens.update(m["SPECIMEN_SOURCE"].dropna().astype(str).str.strip().tolist())

    n = counts["candidate_micro_rows"]
    return {
        "candidate_definition": "RAW_LAB_NAME or RAW_PANEL contains culture/microbiology/organism/susceptibility/gram-stain terms",
        "candidate_micro_rows": int(n),
        "timestamp_coverage": {
            k: (int(v) / n if n else None)
            for k, v in counts.items()
            if k.startswith("nonmissing_")
        },
        "positivity_signal_counts": {
            k: int(v) if v >= MIN_CELL else None
            for k, v in counts.items()
            if "positive_term" in k or "negative_term" in k
        },
        "distinct_loinc_count": len([x for x in loincs if x]),
        "top_loinc_counts": [
            {"loinc": v, "count": int(c)} for v, c in loincs.most_common(20) if c >= MIN_CELL and v
        ],
        "top_specimen_source_counts": [
            {"specimen_source": v, "count": int(c)} for v, c in specimens.most_common(20) if c >= MIN_CELL and v
        ],
        "warning": "Keyword-based microbiology identification and positivity terms are diagnostic only and must not be treated as the frozen culture phenotype without validation.",
    }


def inspect_code(path: Path, label: str) -> dict:
    text = path.read_text(encoding="utf-8", errors="replace")
    lower = text.lower()
    patterns = {
        "references_sepsis_encounter": "sepsis_encounter" in lower,
        "references_med_admin": "med_admin" in lower,
        "references_prescribing": "prescribing" in lower,
        "references_rxnorm": "rxnorm" in lower,
        "references_icu": "icu" in lower,
        "references_admit_time": "admit_time" in lower,
        "references_discharge_time": "discharge_time" in lower,
        "references_lab_result_cm": "lab_result_cm" in lower,
        "references_specimen": "specimen" in lower,
        "references_result_time": "result_time" in lower,
        "references_vasopressor": "vasopressor" in lower,
        "midnight_literal_present": bool(re.search(r"00:00(?::00)?|midnight", lower)),
    }
    out = {"logical_file": label, "relative_path": str(path.relative_to(EXPECTED_ROOT)), **patterns}
    if label == "antibiotic_codes":
        try:
            tree = ast.parse(text)
            objects = []
            for node in tree.body:
                if isinstance(node, (ast.Assign, ast.AnnAssign)):
                    targets = node.targets if isinstance(node, ast.Assign) else [node.target]
                    value = node.value
                    names = [t.id for t in targets if isinstance(t, ast.Name)]
                    for name in names:
                        if any(k in name.lower() for k in ["abx", "antibi", "broad", "rxnorm"]):
                            length = None
                            if isinstance(value, (ast.List, ast.Tuple, ast.Set, ast.Dict)):
                                length = len(value.elts) if hasattr(value, "elts") else len(value.keys)
                            objects.append({"name": name, "literal_length": length})
            out["candidate_mapping_objects"] = objects
        except SyntaxError:
            out["candidate_mapping_objects"] = []
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description="Aggregate-only PSU-to-MIMIC target-trial crosswalk audit.")
    ap.add_argument("data_root", type=Path)
    ap.add_argument("--output-dir", type=Path, default=Path("outputs/psu_target_trial_crosswalk/latest"))
    args = ap.parse_args()
    root = safe_root(args.data_root)
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)

    timing_specs = [
        ("sepsis_encounter", "admit_date", "hospital/encounter anchor"),
        ("sepsis_encounter", "discharge_date", "96h landmark/discharge"),
        ("prescribing", "rx_order_date", "antibiotic order date"),
        ("prescribing", "rx_order_time", "antibiotic order time"),
        ("prescribing", "rx_start_date", "antibiotic start date"),
        ("prescribing", "rx_end_date", "antibiotic end date"),
        ("med_admin", "medadmin_start_date", "administration start date"),
        ("med_admin", "medadmin_start_time", "administration start time"),
        ("med_admin", "medadmin_stop_date", "administration stop date"),
        ("med_admin", "medadmin_stop_time", "administration stop time"),
        ("lab_reduced", "SPECIMEN_DATE", "culture specimen date"),
        ("lab_reduced", "SPECIMEN_TIME", "culture specimen time"),
        ("lab_reduced", "RESULT_DATE", "culture result-availability date"),
        ("lab_reduced", "RESULT_TIME", "culture result-availability time"),
        ("death", "death_date", "30-day mortality"),
        ("sepsis_vital", "measure_date", "trajectory vital date"),
        ("sepsis_vital", "measure_time", "trajectory vital time"),
        ("obs_clin", "obsclin_start_date", "clinical observation date"),
        ("obs_clin", "obsclin_start_time", "clinical observation time"),
    ]
    timing = [timing_row(root, *spec) for spec in timing_specs]
    pd.DataFrame(timing).to_csv(out / "timestamp_crosswalk.csv", index=False)

    category_rows = []
    for table, cols in {
        "sepsis_encounter": ["enc_type", "discharge_status"],
        "prescribing": ["rx_route", "rx_source", "rx_basis"],
        "med_admin": ["medadmin_route", "medadmin_source", "medadmin_type"],
    }.items():
        for row in aggregate_categories(root / TABLES[table], cols):
            category_rows.append({"table": table, **row})
    pd.DataFrame(category_rows).to_csv(out / "aggregate_category_counts.csv", index=False)

    code_summary = [inspect_code(root / rel, label) for label, rel in CODE_FILES.items() if (root / rel).exists()]
    (out / "local_code_logic_summary.json").write_text(json.dumps(code_summary, indent=2), encoding="utf-8")

    micro = scan_microbiology(root / TABLES["lab_reduced"])
    (out / "microbiology_mapping_summary.json").write_text(json.dumps(micro, indent=2), encoding="utf-8")

    encounter_cols = set(pq.ParquetFile(root / TABLES["sepsis_encounter"]).schema_arrow.names)
    summary = {
        "privacy_mode": "aggregate_only_no_ids_no_rows_no_free_text_export",
        "minimum_reported_cell": MIN_CELL,
        "icu_timing": {
            "sepsis_encounter_has_admit_date": "admit_date" in encounter_cols,
            "sepsis_encounter_has_admit_time": "admit_time" in encounter_cols,
            "sepsis_encounter_has_icu_intime": any("icu" in c.lower() and "time" in c.lower() for c in encounter_cols),
            "assessment": "Prebuilt sepsis_encounter alone does not provide an exact ICU intime if both admit_time and ICU-specific time fields are absent.",
        },
        "cohort_source": {
            "prebuilt_sepsis_encounter_present": True,
            "status": "not_frozen",
            "reason": "The target trial is antibiotic/suspected-infection anchored; a prefiltered coded-sepsis cohort must not be assumed equivalent without provenance validation.",
        },
        "antibiotic_source": {
            "prescribing_primary_candidate": True,
            "med_admin_sensitivity_candidate": True,
            "status": "not_frozen_pending_route_timing_and_mapping_review",
        },
        "microbiology": {
            "result_time_fields_present": True,
            "diagnostic_candidate_rows": micro["candidate_micro_rows"],
            "status": "not_frozen_pending_validation_of_positive-clinical-culture representation",
        },
        "next_decision": "Freeze the PSU-to-MIMIC field crosswalk only after exact ICU timing/provenance and microbiology positivity semantics are resolved; do not construct the final external target-trial cohort before that.",
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps({"status": "ok", "privacy_mode": summary["privacy_mode"], "candidate_micro_rows": micro["candidate_micro_rows"]}, indent=2))


if __name__ == "__main__":
    main()
