"""Static guardrails for the corrected MIMIC publication workflow.

These tests require no credentialed MIMIC data. They protect the documented final
v5.7 analysis contract from accidental software/documentation drift.
"""
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def _text(path: str) -> str:
    return (ROOT / path).read_text(encoding="utf-8")


def test_final_vital_correction_and_inference_entry_points_exist():
    assert (ROOT / "scripts" / "repair_v57_vital_covariates.py").exists()
    assert (ROOT / "scripts" / "rerun_inference.py").exists()


def test_primary_ps_excludes_audited_direct_gcs_and_fio2_terms():
    spec = _text("src/sepsis_deescalation/specification.py")
    assert '"gcs_total_48_72h"' in spec
    assert '"fio2_48_72h"' in spec
    assert "PS_EXCLUDED_VITAL_DIRECT_TERMS" in spec
    assert "CANDIDATE_PS_VARS" in spec


def test_ps_preparation_keeps_frozen_mean_fill_and_probability_clipping():
    stats = _text("src/sepsis_deescalation/stats.py")
    assert "x.fillna(mu)" in stats
    assert ".clip(lower=-8, upper=8)" in stats
    assert "x.fillna(0)" in stats
    assert "0.001, 0.999" in stats
    assert "fit_regularized(alpha=0.001" in stats


def test_bootstrap_refits_ps_within_resample():
    stats = _text("src/sepsis_deescalation/stats.py")
    boot = stats.split("def bootstrap_iptw_ci", 1)[1]
    assert "sample = df.iloc" in boot
    assert "fit_stabilized_iptw(sample" in boot


def test_readme_identifies_corrected_inference_as_manuscript_source():
    readme = _text("README.md")
    assert "repair_v57_vital_covariates.py" in readme
    assert "analysis_cohort_vital_corrected.csv" in readme
    assert "vital_corrected_final" in readme
    assert "manuscript primary/secondary" in readme.lower()


def test_mimic_review_records_primary_ci_reporting_rule():
    review = _text("docs/mimic_v57_freeze_review.md")
    assert "primary_secondary_outcomes.csv" in review
    assert "Do not choose between those CIs based on favorability" in review
    assert "mean" in review.lower() and "PS preparation" in review
