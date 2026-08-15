from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import pandas as pd


@dataclass(frozen=True)
class HarmonizationStatus:
    construct: str
    status: str
    mimic_definition: str
    psu_definition: str
    source: str
    limitation: str = ""


def harmonization_table(rows: list[HarmonizationStatus]) -> pd.DataFrame:
    return pd.DataFrame([r.__dict__ for r in rows])


def default_crosswalk() -> pd.DataFrame:
    rows = [
        HarmonizationStatus(
            "Trial time zero",
            "direct",
            "First qualifying systemic IV broad-spectrum prescription coverage start",
            "First qualifying broad-spectrum order start in PRESCRIBING",
            "PRESCRIBING",
        ),
        HarmonizationStatus(
            "Medication administration sensitivity",
            "PSU enrichment",
            "Not available from MIMIC prescription phenotype",
            "First/continued broad-spectrum medication actually administered",
            "MED_ADMIN",
            "Different measurement construct; use as a measurement-validation sensitivity, not silently as the primary replication.",
        ),
        HarmonizationStatus(
            "Microbiology specimen time",
            "approximate",
            "microbiologyevents charttime/chartdate",
            "LAB_RESULT_CM SPECIMEN_DATE/TIME or validated local microbiology table",
            "LAB_RESULT_CM/local microbiology",
        ),
        HarmonizationStatus(
            "Microbiology result availability",
            "approximate",
            "microbiologyevents storetime/storedate",
            "LAB_RESULT_CM RESULT_DATE/TIME or validated local final-result timestamp",
            "LAB_RESULT_CM/local microbiology",
            "PCORnet permits RESULT_DATE substitution when true result date is unavailable. PSU ETL semantics must be audited before interpreting it as result availability.",
        ),
        HarmonizationStatus(
            "ICU admission/exit",
            "local extension required",
            "MIMIC icustays intime/outtime",
            "PSU local ICU/ADT location table with timestamps",
            "local ICU table",
            "Core PCORnet ENCOUNTER does not by itself guarantee ICU-location timestamps.",
        ),
        HarmonizationStatus(
            "Vital signs",
            "approximate",
            "MIMIC chartevents",
            "VITAL and/or OBS_CLIN with encounter and time",
            "VITAL/OBS_CLIN",
        ),
        HarmonizationStatus(
            "Death",
            "direct/approximate",
            "MIMIC deathtime/DOD",
            "DEATH.DEATH_DATE plus encounter death status when appropriate",
            "DEATH/ENCOUNTER",
        ),
    ]
    return harmonization_table(rows)


def require_columns(df: pd.DataFrame, columns: list[str], table: str) -> None:
    missing = [c for c in columns if c not in df.columns]
    if missing:
        raise ValueError(f"{table} missing required columns: {missing}")


def first_present(mapping: dict[str, Any], candidates: list[str]) -> str | None:
    for key in candidates:
        value = mapping.get(key)
        if value:
            return str(value)
    return None
