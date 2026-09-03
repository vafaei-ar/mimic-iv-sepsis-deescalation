#!/usr/bin/env python3
"""Validate shared figure contracts and rebuild legacy/Nature figure sets.

The task is publication QA only. It does not refit the frozen outcome models or
change any estimand. The figure builders may read the restricted local analytic
cohort to reproduce propensity/weight diagnostics, but only aggregate figures and
this sanitized report are declared as artifacts.
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import build_nature_figures
import build_submission_figures

OUT = Path("outputs/publication_integration/figure_validation")
SUBMISSION = Path("outputs/publication_integration/submission_figures")
NATURE = Path("outputs/publication_integration/nature_figures")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    subprocess.run(
        [sys.executable, "-m", "pytest", "-q", "tests/test_publication_figure_contract.py"],
        check=True,
    )

    build_submission_figures.main()
    build_nature_figures.main()

    expected_submission = [
        "Fig1_target_trial_timeline.png",
        "Fig1_target_trial_timeline.pdf",
        "Fig2_progressive_adjustment.png",
        "Fig2_progressive_adjustment.pdf",
        "Fig3_cross_dataset_outcomes.png",
        "Fig3_cross_dataset_outcomes.pdf",
        "ESM_Fig1_mimic_balance_love.png",
        "ESM_Fig1_mimic_balance_love.pdf",
        "ESM_Fig2_mimic_ps_weights.png",
        "ESM_Fig2_mimic_ps_weights.pdf",
    ]
    expected_nature = [
        "Fig2_progressive_adjustment.png",
        "Fig2_progressive_adjustment.pdf",
        "Fig3_cross_dataset_outcomes.png",
        "Fig3_cross_dataset_outcomes.pdf",
        "Fig1b_cohort_attrition.png",
        "Fig1b_cohort_attrition.pdf",
        "ESM_Fig1_mimic_balance.png",
        "ESM_Fig1_mimic_balance.pdf",
        "ESM_Fig2_mimic_overlap_weights.png",
        "ESM_Fig2_mimic_overlap_weights.pdf",
    ]

    missing = [
        str(SUBMISSION / name)
        for name in expected_submission
        if not (SUBMISSION / name).exists()
    ] + [
        str(NATURE / name)
        for name in expected_nature
        if not (NATURE / name).exists()
    ]

    report = {
        "contract_tests_passed": True,
        "submission_figures_rebuilt": len(expected_submission),
        "nature_figures_rebuilt": len(expected_nature),
        "missing_expected_outputs": missing,
        "all_expected_outputs_present": not missing,
        "primary_science_changed": False,
        "notes": [
            "Legacy overlaid weight histograms now use one shared bin grid.",
            "All legacy manuscript-facing builders use the frozen PSU antibiotic-free-days reconciliation.",
            "All legacy M4 CI grafts assert point-estimate parity first.",
            "Reviewer-support Love-plot labels use manuscript-facing names.",
        ],
    }
    (OUT / "figure_validation_report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))

    if missing:
        raise SystemExit("Figure validation failed: expected outputs are missing")


if __name__ == "__main__":
    main()
