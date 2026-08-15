#!/usr/bin/env python
from __future__ import annotations

import argparse
import sys

from sepsis_deescalation.config import load_config, resolve_mimic_source
from sepsis_deescalation.mimic_io import get_columns, table_exists

REQUIRED_TABLES = {
    "hosp/admissions.csv.gz": ["subject_id", "hadm_id", "admittime", "dischtime", "deathtime", "race"],
    "hosp/patients.csv.gz": ["subject_id", "gender", "anchor_age", "dod"],
    "hosp/prescriptions.csv.gz": ["subject_id", "hadm_id", "starttime", "stoptime", "drug", "route"],
    "hosp/microbiologyevents.csv.gz": ["subject_id", "hadm_id", "org_name", "spec_type_desc", "test_name"],
    "hosp/diagnoses_icd.csv.gz": ["subject_id", "hadm_id", "icd_code", "icd_version"],
    "hosp/d_icd_diagnoses.csv.gz": ["icd_code", "icd_version", "long_title"],
    "hosp/procedures_icd.csv.gz": ["subject_id", "hadm_id", "chartdate", "icd_code", "icd_version"],
    "hosp/d_icd_procedures.csv.gz": ["icd_code", "icd_version", "long_title"],
    "hosp/d_labitems.csv.gz": ["itemid", "label"],
    "hosp/labevents.csv.gz": ["subject_id", "hadm_id", "itemid", "charttime", "valuenum"],
    "icu/icustays.csv.gz": ["subject_id", "hadm_id", "stay_id", "first_careunit", "intime", "outtime"],
    "icu/d_items.csv.gz": ["itemid", "label"],
    "icu/chartevents.csv.gz": ["subject_id", "hadm_id", "stay_id", "charttime", "itemid", "valuenum"],
    "icu/inputevents.csv.gz": ["subject_id", "hadm_id", "stay_id", "starttime", "endtime", "itemid"],
    "icu/outputevents.csv.gz": ["subject_id", "hadm_id", "stay_id", "charttime", "itemid", "value"],
}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/mimic.yaml")
    args = parser.parse_args()
    cfg = load_config(args.config)
    source = resolve_mimic_source(cfg)
    print(f"MIMIC source: {source}")
    failed = False
    for rel, required in REQUIRED_TABLES.items():
        if not table_exists(source, rel):
            print(f"FAIL missing table: {rel}")
            failed = True
            continue
        cols = set(get_columns(source, rel))
        missing = [c for c in required if c not in cols]
        if missing:
            print(f"FAIL {rel}: missing columns {missing}")
            failed = True
        else:
            print(f"OK   {rel}")
    micro_cols = set(get_columns(source, "hosp/microbiologyevents.csv.gz"))
    if not ({"storetime", "storedate"} & micro_cols):
        print("FAIL microbiologyevents lacks both storetime and storedate; primary result-availability phenotype cannot run.")
        failed = True
    else:
        print("OK   microbiology result-availability timestamp present")
    if failed:
        sys.exit(2)
    print("Preflight passed.")


if __name__ == "__main__":
    main()
