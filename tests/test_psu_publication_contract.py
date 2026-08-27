"""Static guardrails for the frozen PSU publication pipeline.

These tests do not access restricted PSU data and therefore can run in public CI. Their
purpose is to catch accidental software drift in the publication entry points before a
reviewer or collaborator runs the real-data pipeline on the approved machine.
"""
import importlib.util
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"


def _text(name: str) -> str:
    return (SCRIPTS / name).read_text(encoding="utf-8")


def _load_script_module(name: str):
    path = SCRIPTS / name
    spec = importlib.util.spec_from_file_location(path.stem, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_final_psu_entry_points_are_present():
    required = {
        "audit_psu_final_covariate_freeze.py",
        "audit_psu_ps_balance.py",
        "audit_psu_final_outcome_freeze.py",
        "run_psu_point_estimates.py",
        "run_psu_bootstrap_inference.py",
        "run_psu_prespecified_robustness.py",
        "run_psu_prespecified_robustness_bootstrap.py",
    }
    missing = sorted(name for name in required if not (SCRIPTS / name).exists())
    assert not missing, f"Missing frozen PSU entry points: {missing}"


def test_point_estimate_insertion_marker_is_unique():
    ps = _text("audit_psu_ps_balance.py")
    marker = (
        "    max_pre=float(bdf.abs_pre_smd.max()); max_post=float(bdf.abs_post_smd.max()); "
        "worst=str(bdf.iloc[0].variable)\n"
    )
    assert ps.count(marker) == 1, "Frozen PS insertion marker must appear exactly once"


def test_bootstrap_publication_settings_are_fixed():
    boot = _text("run_psu_bootstrap_inference.py")
    assert "REPS = 1000" in boot
    assert "SEED = 20260826" in boot
    assert "fit_regularized(alpha=0.001" in boot
    assert "0.001, 0.999" in boot


def test_bootstrap_refits_propensity_score_inside_worker():
    boot = _text("run_psu_bootstrap_inference.py")
    worker = boot.split("def _worker", 1)[1].split("def run_bootstrap", 1)[0]
    assert "sm.GLM" in worker
    assert "model.fit" in worker
    assert "sw = np.where" in worker


def test_only_prespecified_psu_robustness_variants_are_named():
    rob = _text("run_psu_prespecified_robustness.py")
    assert '"medadmin_exposure"' in rob
    assert '"lenient_landmark"' in rob
    assert "MEDADMIN_OVERRIDE" in rob
    assert "STRICT_COHORT" in rob and "LENIENT_COHORT" in rob


def test_reviewer_documentation_exists_and_uses_correct_terminology():
    walkthrough = (ROOT / "docs" / "psu_analysis_walkthrough.md").read_text(encoding="utf-8")
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    combined = walkthrough + "\n" + readme
    assert "modified external replication" in combined
    assert "ordered systemic broad-spectrum antibiotic proxy" in combined
    assert "not an exact" in combined.lower()
    assert "refit" in walkthrough.lower() and "bootstrap" in walkthrough.lower()


def test_obsolete_execution_placeholders_are_removed():
    assert not list((ROOT / "docs").glob(".placeholder_psu_point_estimates*"))


def test_parity_report_rows_are_json_safe_for_numpy_scalars():
    parity = _load_script_module("check_psu_publication_parity.py")
    checks = []
    parity._record(checks, "integer", np.int64(19841), 19841)
    parity._record(checks, "float", np.float64(-0.0256), -0.0256, 1e-6)
    json.dumps({"checks": checks})
    assert checks[0]["observed"] == 19841
    assert checks[1]["passed"] is True


def test_parity_tolerances_are_narrow_and_frozen():
    """Prevent a future edit from silently weakening the real-data parity contract."""
    parity = _load_script_module("check_psu_publication_parity.py")
    assert parity.POINT_ATOL == 1e-5
    assert parity.BOOTSTRAP_RD_ATOL == 0.001
    assert parity.BOOTSTRAP_RR_ATOL == 0.005
