#!/usr/bin/env python3
"""Shared data contracts for manuscript-facing publication figures.

This module centralizes the manuscript-facing values and validation logic that must
remain identical across the submission, reviewer-support, reviewer-final, and Nature
figure builders. It does not fit models or read restricted data by itself.
"""
from __future__ import annotations

from collections.abc import Iterable

import numpy as np
import pandas as pd

PUBLICATION_SECONDARY_OVERRIDES = {
    ("PSU", "Antibiotic-free days"): {
        "estimate": 3.16,
        "ci95_low": 2.85,
        "ci95_high": 3.47,
    },
}

PRETTY_LABELS = {
    "broad_abx_hours_pre72": "Broad-spectrum antibiotic hours, pre-72 h",
    "systemic_abx_hours_pre72": "Systemic antibiotic hours, pre-72 h",
    "antipseudomonal_pre72": "Antipseudomonal therapy, pre-72 h",
    "broad_abx_agents_pre72": "Broad-spectrum agents, pre-72 h",
    "micro_records_pre72": "Microbiology records, pre-72 h",
    "anaerobic_coverage_pre72": "Anaerobic coverage, pre-72 h",
    "strict_culture_records_pre72": "Culture records, pre-72 h",
    "carbapenem_pre72": "Carbapenem use, pre-72 h",
    "distinct_specimen_types_pre72": "Distinct specimen types, pre-72 h",
    "cardiac_icu": "Cardiac ICU",
    "respiratory_culture_pre72": "Respiratory culture, pre-72 h",
    "temperature_48_72h": "Temperature, 48-72 h",
    "vent_proc": "Mechanical ventilation procedure",
    "sterile_fluid_culture_pre72": "Sterile-fluid culture, pre-72 h",
    "sofa_like_change_pre72": "SOFA-like change, pre-72 h",
    "sofa_like_48_72h": "SOFA-like score, 48-72 h",
    "fever_last12h_pre72": "Fever in prior 12 h, pre-72 h",
    "repeat_micro_48_72h": "Repeat microbiology, 48-72 h",
    "blood_culture_pre72": "Blood culture, pre-72 h",
    "hr_max_pre72": "Maximum heart rate, pre-72 h",
    "platelet_late_worst_48_72h": "Platelet count, worst 48-72 h",
    "sofa_like_improved_pre72": "SOFA-like improvement, pre-72 h",
    "lactate_last_pre72": "Last lactate, pre-72 h",
    "bilirubin_rising_pre72": "Rising bilirubin, pre-72 h",
    "wbc_late_last_48_72h": "Last WBC count, 48-72 h",
    "sicu": "Surgical ICU",
    "systemic_abx_agents_pre72": "Systemic antibiotic agents, pre-72 h",
    "lactate_rising_pre72": "Rising lactate, pre-72 h",
    "hours_admit_to_icu": "Admission-to-ICU interval, h",
    "micu": "Medical ICU",
    "white_blood_cells_last_pre72": "Last WBC count, pre-72 h",
    "wbc_rising_pre72": "Rising WBC count, pre-72 h",
    "severity_pre72": "Severity index, pre-72 h",
    "bilirubin_late_worst_48_72h": "Bilirubin, worst 48-72 h",
    "vasopressor_stopped_before_72h": "Vasopressor stopped before 72 h",
}

ABBREVIATIONS = {
    "icu": "ICU",
    "wbc": "WBC",
    "sofa": "SOFA",
    "iv": "IV",
    "bmi": "BMI",
    "spo2": "SpO2",
    "map": "MAP",
    "rr": "RR",
    "hr": "HR",
}


def apply_publication_secondary_overrides(sec: pd.DataFrame) -> pd.DataFrame:
    """Return manuscript-facing secondary outcomes with frozen reconciliations."""
    out = sec.copy()
    for (dataset, outcome), values in PUBLICATION_SECONDARY_OVERRIDES.items():
        mask = (out["dataset"] == dataset) & (out["outcome"] == outcome)
        if int(mask.sum()) != 1:
            raise RuntimeError(
                f"Expected exactly one row for publication override: {dataset}, {outcome}"
            )
        for column, value in values.items():
            out.loc[mask, column] = value
    return out


def prepare_progressive_mortality(
    mortality: pd.DataFrame,
    progressive: pd.DataFrame,
    *,
    tolerance: float = 1e-9,
) -> pd.DataFrame:
    """Apply the designated primary M4 CI only after point-estimate parity passes."""
    primary_rows = mortality.loc[
        mortality["dataset_analysis"].astype(str).str.startswith("MIMIC-IV primary")
    ]
    if len(primary_rows) != 1:
        raise RuntimeError(
            f"Expected exactly one MIMIC-IV primary mortality row, found {len(primary_rows)}"
        )
    if progressive.empty:
        raise RuntimeError("Progressive mortality table is empty")

    primary = primary_rows.iloc[0]
    out = progressive.copy()
    last = out.index[-1]
    progressive_point = float(out.loc[last, "risk_difference"])
    primary_point = float(primary["mortality_rd"])
    if abs(progressive_point - primary_point) > tolerance:
        raise RuntimeError(
            "M4 point estimate disagrees with the harmonized primary result: "
            f"progressive={progressive_point}, primary={primary_point}"
        )
    out.loc[last, "rd_lower_95"] = float(primary["rd_ci95_low"])
    out.loc[last, "rd_upper_95"] = float(primary["rd_ci95_high"])
    return out


def pretty_label(name: str) -> str:
    """Convert internal covariate names to manuscript-facing labels."""
    if name in PRETTY_LABELS:
        return PRETTY_LABELS[name]
    words = [ABBREVIATIONS.get(w.lower(), w) for w in str(name).replace("_", " ").split()]
    text = " ".join(words)
    return text[:1].upper() + text[1:]


def shared_histogram_edges(groups: Iterable[Iterable[float]], bins: int = 30) -> np.ndarray:
    """Build one numeric bin grid shared by every overlaid histogram group."""
    arrays = []
    for values in groups:
        x = np.asarray(values, dtype=float)
        x = x[np.isfinite(x)]
        if x.size:
            arrays.append(x)
    if not arrays:
        raise ValueError("Cannot build histogram edges from empty/non-finite data")
    pooled = np.concatenate(arrays)
    lo = float(pooled.min())
    hi = float(pooled.max())
    if hi <= lo:
        delta = max(abs(lo) * 0.01, 0.5)
        lo -= delta
        hi += delta
    return np.linspace(lo, hi, int(bins) + 1)
