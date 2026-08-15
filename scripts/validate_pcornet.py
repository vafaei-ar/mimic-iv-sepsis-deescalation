#!/usr/bin/env python
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from sepsis_deescalation.config import load_config
from sepsis_deescalation.pcornet_io import columns_for, resolve_path

CORE_REQUIREMENTS = {
    "demographic": ["patid", "birth_date", "sex", "race"],
    "encounter": ["patid", "encounterid", "admit_date", "enc_type"],
    "prescribing": ["patid", "encounterid"],
    "med_admin": ["patid", "encounterid", "medadmin_start_date"],
    "lab_result_cm": ["patid", "specimen_date", "result_date"],
    "death": ["patid", "death_date"],
}


def _check_mapping(name: str, path: Path, mapping: dict[str, str], required_keys: list[str]) -> tuple[bool, list[str]]:
    cols = {c.lower() for c in columns_for(path)}
    errors = []
    for key in required_keys:
        col = mapping.get(key)
        if not col:
            errors.append(f"mapping '{key}' not configured")
        elif col.lower() not in cols:
            errors.append(f"mapped column '{col}' for '{key}' not present")
    return not errors, errors


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate PSU/PCORnet external-validation inputs.")
    parser.add_argument("--config", default="config/pcornet_psu.local.yaml")
    args = parser.parse_args()
    cfg = load_config(args.config)
    failed = False

    print("PCORnet/PSU preflight")
    for key, required in CORE_REQUIREMENTS.items():
        try:
            path = resolve_path(cfg, key)
        except KeyError as exc:
            print(f"FAIL {key}: {exc}")
            failed = True
            continue
        if not path.exists():
            print(f"FAIL {key}: missing file {path}")
            failed = True
            continue
        cols = {c.lower() for c in columns_for(path)}
        missing = [c for c in required if c not in cols]
        if missing:
            print(f"WARN {key}: standard fields not all present: {missing}. Local mappings may still resolve this.")
        else:
            print(f"OK   {key}: {path}")

    # ICU timestamps are scientifically required for the harmonized primary cohort.
    try:
        icu_path = resolve_path(cfg, "icu_stays")
        if not icu_path.exists():
            raise FileNotFoundError(icu_path)
        ok, errors = _check_mapping(
            "icu_stays", icu_path, cfg.get("icu", {}), ["patid", "encounterid", "start_datetime", "end_datetime"]
        )
        if ok:
            print("OK   local ICU table and timestamp mappings")
        else:
            print("FAIL local ICU mappings: " + "; ".join(errors)); failed = True
    except Exception as exc:
        print(f"FAIL local ICU table: {exc}")
        failed = True

    meds = cfg.get("medications", {})
    primary = meds.get("primary_source", "prescribing")
    if primary != "prescribing":
        print("WARN primary medication source is not PRESCRIBING; this changes comparability with MIMIC order/prescription exposure.")
    try:
        p = resolve_path(cfg, "prescribing")
        ok, errors = _check_mapping(
            "prescribing", p, meds.get("prescribing_columns", {}),
            ["patid", "encounterid", "raw_name", "start_date", "end_date"]
        )
        if not ok:
            print("FAIL prescribing mappings: " + "; ".join(errors)); failed = True
        else:
            print("OK   prescribing mappings")
    except Exception as exc:
        print(f"FAIL prescribing mapping check: {exc}"); failed = True

    micro = cfg.get("microbiology", {})
    micro_source = micro.get("source", "local")
    if micro_source == "local":
        try:
            p = resolve_path(cfg, "microbiology_local")
            ok, errors = _check_mapping(
                "microbiology_local", p, micro.get("local_columns", {}),
                ["patid", "encounterid", "specimen_datetime", "result_datetime", "test_name", "organism"]
            )
            if not ok:
                print("FAIL local microbiology mappings: " + "; ".join(errors)); failed = True
            else:
                print("OK   local microbiology mapping")
        except Exception as exc:
            print(f"FAIL local microbiology source: {exc}"); failed = True
    else:
        print("WARN using LAB_RESULT_CM for microbiology. Confirm PSU ETL semantics for RESULT_DATE/TIME before treating them as result availability.")
        try:
            p = resolve_path(cfg, "lab_result_cm")
            ok, errors = _check_mapping(
                "lab_result_cm", p, micro.get("lab_result_cm_columns", {}),
                ["patid", "specimen_date", "result_date", "test_name", "result_qual", "raw_result"]
            )
            if not ok:
                print("FAIL LAB_RESULT_CM microbiology mappings: " + "; ".join(errors)); failed = True
        except Exception as exc:
            print(f"FAIL LAB_RESULT_CM microbiology mapping: {exc}"); failed = True

    if failed:
        print("\nPreflight FAILED. Correct mappings before running the external validation.")
        sys.exit(2)
    print("\nPreflight passed. The data structure is sufficient to begin the harmonized external run.")


if __name__ == "__main__":
    main()
