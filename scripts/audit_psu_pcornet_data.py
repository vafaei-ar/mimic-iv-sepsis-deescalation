#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq

EXPECTED_ROOT = Path("/home/asadr/Depts/PHS/PATH_CDM/Hwang_Bonavia").resolve()

PARQUET_CANDIDATES = {
    "sepsis_encounter": "PCORnet/parquet/sepsis_encounter.parquet",
    "sepsis_demographic": "PCORnet/parquet/sepsis_demographic.parquet",
    "sepsis_diagnosis": "PCORnet/parquet/sepsis_diagnosis.parquet",
    "sepsis_vital": "PCORnet/parquet/sepsis_vital.parquet",
    "condition": "PCORnet/parquet/condition.parquet",
    "death": "PCORnet/parquet/death.parquet",
    "prescribing": "PCORnet/parquet/prescribing.parquet",
    "med_admin": "PCORnet/parquet/med_admin.parquet",
    "obs_clin": "PCORnet/parquet/obs_clin.parquet",
    "lab_reduced": "PCORnet/parquet/lab/lab_reduced.parquet",
}

CODE_FILES = {
    "cohort_extractor": "PCORnet/code/extract/01_cohort.py",
    "antibiotic_codes": "PCORnet/code/config/codes_antibiotics.py",
    "settings": "PCORnet/code/config/settings.py",
}

CODE_KEYWORDS = [
    "sepsis_encounter",
    "prescribing",
    "med_admin",
    "rxnorm",
    "encounterid",
    "patid",
    "admit_date",
    "admit_time",
    "discharge_date",
    "discharge_time",
    "lab_result_cm",
    "specimen",
    "result_date",
    "result_time",
    "vasopressor",
    "icu",
]

IDENTIFIER_LIKE = {
    "patid",
    "patient_id",
    "encounterid",
    "encounter_id",
    "providerid",
    "provider_id",
}


def _safe_root(root: Path) -> Path:
    root = root.expanduser().resolve()
    if root != EXPECTED_ROOT:
        raise ValueError(f"Refusing unexpected data root: {root}")
    return root


def _parquet_metadata(path: Path, logical_name: str) -> tuple[dict, list[dict]]:
    pf = pq.ParquetFile(path)
    metadata = pf.metadata
    schema = pf.schema_arrow
    total_rows = int(metadata.num_rows)
    rows: list[dict] = []

    for idx, field in enumerate(schema):
        null_count = 0
        null_known = True
        for rg_idx in range(metadata.num_row_groups):
            col_meta = metadata.row_group(rg_idx).column(idx)
            stats = col_meta.statistics
            if stats is None or stats.null_count is None:
                null_known = False
                break
            null_count += int(stats.null_count)
        name_lower = field.name.lower()
        rows.append(
            {
                "logical_table": logical_name,
                "file": str(path.relative_to(EXPECTED_ROOT)),
                "row_count": total_rows,
                "column": field.name,
                "dtype": str(field.type),
                "identifier_like": name_lower in IDENTIFIER_LIKE or name_lower.endswith("_id"),
                "null_count_metadata": null_count if null_known else None,
                "null_fraction_metadata": (null_count / total_rows) if null_known and total_rows else None,
            }
        )

    table_summary = {
        "logical_table": logical_name,
        "file": str(path.relative_to(EXPECTED_ROOT)),
        "row_count": total_rows,
        "n_columns": len(schema),
        "n_row_groups": int(metadata.num_row_groups),
        "columns": [field.name for field in schema],
    }
    return table_summary, rows


def _scan_code(path: Path, logical_name: str) -> dict:
    text = path.read_text(encoding="utf-8", errors="replace")
    lower = text.lower()
    hits = {keyword: lower.count(keyword.lower()) for keyword in CODE_KEYWORDS if keyword.lower() in lower}
    return {
        "logical_file": logical_name,
        "file": str(path.relative_to(EXPECTED_ROOT)),
        "exists": True,
        "keyword_counts": hits,
        "note": "Counts only; source text is intentionally not exported.",
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Read-only metadata audit of local PSU/PCORnet sepsis data. No row-level data are exported."
    )
    parser.add_argument("data_root", type=Path)
    parser.add_argument("--output-dir", type=Path, default=Path("outputs/psu_pcornet_audit/latest"))
    args = parser.parse_args()

    root = _safe_root(args.data_root)
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)

    inventory_rows: list[dict] = []
    table_summaries: list[dict] = []
    schema_rows: list[dict] = []

    for logical_name, rel in PARQUET_CANDIDATES.items():
        path = root / rel
        inventory_rows.append(
            {
                "kind": "parquet",
                "logical_name": logical_name,
                "relative_path": rel,
                "exists": path.exists(),
                "size_bytes": path.stat().st_size if path.exists() else None,
            }
        )
        if path.exists():
            summary, rows = _parquet_metadata(path, logical_name)
            table_summaries.append(summary)
            schema_rows.extend(rows)

    code_summaries: list[dict] = []
    for logical_name, rel in CODE_FILES.items():
        path = root / rel
        inventory_rows.append(
            {
                "kind": "code",
                "logical_name": logical_name,
                "relative_path": rel,
                "exists": path.exists(),
                "size_bytes": path.stat().st_size if path.exists() else None,
            }
        )
        if path.exists():
            code_summaries.append(_scan_code(path, logical_name))

    # Presence-only checks for additional local files that may matter later.
    presence_only = [
        "PCORnet/sepsis_procedures.sas7bdat",
        "PCORnet/Full/lab_result_cm.sas7bdat",
        "PCORnet/Full/prescribing.sas7bdat",
        "PCORnet/Full/med_admin.sas7bdat",
        "PCORnet/Full/obs_clin.sas7bdat",
        "PCORnet/Full/obs_gen.sas7bdat",
        "PCORnet/PCORnet-Common-Data-Model-v61-2023_04_031.pdf",
    ]
    for rel in presence_only:
        path = root / rel
        inventory_rows.append(
            {
                "kind": "presence_only",
                "logical_name": Path(rel).stem,
                "relative_path": rel,
                "exists": path.exists(),
                "size_bytes": path.stat().st_size if path.exists() else None,
            }
        )

    pd.DataFrame(inventory_rows).to_csv(out / "file_inventory.csv", index=False)
    pd.DataFrame(schema_rows).to_csv(out / "parquet_schema_summary.csv", index=False)
    (out / "code_reference_summary.json").write_text(json.dumps(code_summaries, indent=2), encoding="utf-8")

    columns_by_table = {x["logical_table"]: {c.lower() for c in x["columns"]} for x in table_summaries}

    def has(table: str, *tokens: str) -> bool:
        cols = columns_by_table.get(table, set())
        return all(any(token in col for col in cols) for token in tokens)

    candidate_rows = [
        {
            "domain": "encounter / landmark timing",
            "preferred_local_source": "sepsis_encounter",
            "available": "sepsis_encounter" in columns_by_table,
            "mapping_signal": "encounter-like date/time columns present" if has("sepsis_encounter", "admit") else "needs manual mapping audit",
        },
        {
            "domain": "antibiotic orders",
            "preferred_local_source": "prescribing",
            "available": "prescribing" in columns_by_table,
            "mapping_signal": "Rx/order table available; route/timing columns require audit",
        },
        {
            "domain": "antibiotic administrations",
            "preferred_local_source": "med_admin",
            "available": "med_admin" in columns_by_table,
            "mapping_signal": "administration table available for sensitivity/reclassification",
        },
        {
            "domain": "mortality",
            "preferred_local_source": "death",
            "available": "death" in columns_by_table,
            "mapping_signal": "death table available",
        },
        {
            "domain": "vitals",
            "preferred_local_source": "sepsis_vital / obs_clin",
            "available": ("sepsis_vital" in columns_by_table) or ("obs_clin" in columns_by_table),
            "mapping_signal": "vital/clinical observation tables available",
        },
        {
            "domain": "laboratory covariates",
            "preferred_local_source": "lab_reduced / LAB_RESULT_CM",
            "available": "lab_reduced" in columns_by_table or (root / "PCORnet/Full/lab_result_cm.sas7bdat").exists(),
            "mapping_signal": "lab source available; specimen/result timing must be verified",
        },
        {
            "domain": "microbiology culture result availability",
            "preferred_local_source": "LAB_RESULT_CM or local microbiology extension",
            "available": (root / "PCORnet/Full/lab_result_cm.sas7bdat").exists(),
            "mapping_signal": "presence confirmed only; organism positivity and true result-availability timestamps remain unresolved",
        },
    ]
    pd.DataFrame(candidate_rows).to_csv(out / "candidate_sources.csv", index=False)

    summary = {
        "data_root": str(root),
        "privacy_mode": "metadata_only_no_row_level_export",
        "n_parquet_candidates_found": sum(1 for x in table_summaries),
        "tables": table_summaries,
        "code_files_scanned": code_summaries,
        "publication_replication_status": {
            "cohort_source_not_yet_frozen": True,
            "microbiology_result_time_mapping_not_yet_frozen": True,
            "icu_timing_mapping_not_yet_frozen": True,
            "prescribing_vs_med_admin_primary_not_yet_finalized": True,
        },
        "next_step": "Use this aggregate audit to freeze a PSU-to-MIMIC field crosswalk before any patient-level cohort construction.",
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print(json.dumps({
        "status": "ok",
        "privacy_mode": summary["privacy_mode"],
        "n_parquet_candidates_found": summary["n_parquet_candidates_found"],
        "output_dir": str(out),
    }, indent=2))


if __name__ == "__main__":
    main()
