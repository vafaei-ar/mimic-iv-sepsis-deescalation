"""Public-CI guardrails for manuscript-facing figure data contracts."""
from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
COMMON = ROOT / "scripts" / "publication_figure_common.py"

spec = importlib.util.spec_from_file_location("publication_figure_common", COMMON)
common = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(common)


def test_psu_antibiotic_free_days_override_is_publication_locked():
    sec = pd.DataFrame(
        [
            {
                "dataset": "PSU",
                "outcome": "Antibiotic-free days",
                "estimate": 3.157555,
                "ci95_low": 2.837738,
                "ci95_high": 3.467932,
            },
            {
                "dataset": "MIMIC-IV",
                "outcome": "Antibiotic-free days",
                "estimate": 1.75,
                "ci95_low": 0.33,
                "ci95_high": 3.34,
            },
        ]
    )
    out = common.apply_publication_secondary_overrides(sec)
    row = out.loc[(out["dataset"] == "PSU") & (out["outcome"] == "Antibiotic-free days")].iloc[0]
    assert row["estimate"] == 3.16
    assert row["ci95_low"] == 2.85
    assert row["ci95_high"] == 3.47


def test_progressive_m4_requires_point_estimate_parity_before_ci_graft():
    mortality = pd.DataFrame(
        [{
            "dataset_analysis": "MIMIC-IV primary",
            "mortality_rd": 0.0083758268262657,
            "rd_ci95_low": -0.0439033352268441,
            "rd_ci95_high": 0.056713649179014,
        }]
    )
    progressive = pd.DataFrame(
        {
            "risk_difference": [-0.0414, -0.0319, -0.0201, 0.0083758268262657],
            "rd_lower_95": [-0.0597, -0.0501, -0.0395, -0.0551],
            "rd_upper_95": [-0.0229, -0.0139, -0.0010, 0.0563],
        }
    )
    out = common.prepare_progressive_mortality(mortality, progressive)
    assert out.iloc[-1]["rd_lower_95"] == mortality.iloc[0]["rd_ci95_low"]
    assert out.iloc[-1]["rd_upper_95"] == mortality.iloc[0]["rd_ci95_high"]

    bad = progressive.copy()
    bad.loc[bad.index[-1], "risk_difference"] = 0.01
    with pytest.raises(RuntimeError, match="M4 point estimate disagrees"):
        common.prepare_progressive_mortality(mortality, bad)


def test_shared_histogram_edges_use_one_grid_for_all_groups():
    g1 = np.array([0.2, 0.5, 1.0, 4.0])
    g0 = np.array([0.1, 0.3, 0.8, 30.0])
    edges = common.shared_histogram_edges([g1, g0], bins=30)
    assert len(edges) == 31
    assert edges[0] == 0.1
    assert edges[-1] == 30.0
    assert np.all(np.diff(edges) > 0)


def test_pretty_label_hides_internal_covariate_names():
    assert common.pretty_label("broad_abx_hours_pre72") == "Broad-spectrum antibiotic hours, pre-72 h"
    assert common.pretty_label("sofa_like_48_72h").startswith("SOFA-like")
