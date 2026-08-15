#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REQUIRED = [
    "run_manifest.json",
    "run_summary.json",
    "tables/cohort_flow.csv",
    "tables/table1_unweighted.csv",
    "tables/table1_weighted.csv",
    "tables/primary_secondary_outcomes.csv",
    "tables/progressive_adjustment.csv",
    "tables/mortality_sensitivity_analyses.csv",
    "diagnostics/balance_before_after.csv",
    "diagnostics/weight_summary.csv",
    "diagnostics/propensity_score_distribution.csv",
    "audits/microbiology_counts.csv",
    "figure_data/figure2_progressive_adjustment.csv",
    "figure_data/figure3_mortality_sensitivities.csv",
]


def main() -> None:
    parser = argparse.ArgumentParser(description="Audit a completed MIMIC scripted run before manuscript use.")
    parser.add_argument("run_dir")
    args = parser.parse_args()
    root = Path(args.run_dir)
    errors, warnings = [], []
    for rel in REQUIRED:
        if not (root / rel).exists():
            errors.append(f"missing required output: {rel}")

    if errors:
        for e in errors: print("FAIL", e)
        sys.exit(2)

    summary = json.loads((root / "run_summary.json").read_text())
    if summary.get("n", 0) < 1000: warnings.append(f"analytic cohort unexpectedly small: n={summary.get('n')}")
    if summary.get("deescalated_stopped", 0) < 100: errors.append("too few de-escalated/stopped patients for stable primary analysis")
    if summary.get("continued", 0) < 100: errors.append("too few continued-broad patients for stable primary analysis")
    if summary.get("max_post_smd") is not None and summary["max_post_smd"] > 0.10: errors.append(f"post-weighting max absolute SMD exceeds 0.10: {summary['max_post_smd']:.3f}")

    weights = pd.read_csv(root / "diagnostics/weight_summary.csv")
    if not np.isfinite(weights.select_dtypes(include="number").to_numpy()).all(): errors.append("non-finite values in weight diagnostics")
    overall = weights.loc[weights["group"] == "overall"]
    if len(overall) and float(overall.iloc[0]["max_weight"]) > 50: warnings.append(f"large maximum stabilized weight: {float(overall.iloc[0]['max_weight']):.2f}")

    outcomes = pd.read_csv(root / "tables/primary_secondary_outcomes.csv")
    if "bootstrap_success" in outcomes:
        bad = outcomes.loc[outcomes["bootstrap_success"] < 0.95 * outcomes["bootstrap_success"].max()]
        if len(bad): warnings.append("some outcome bootstraps had materially fewer successful resamples; inspect bootstrap diagnostics")

    micro = pd.read_csv(root / "audits/microbiology_counts.csv")
    if not (micro["metric"] == "positive_missing_result_time").any(): warnings.append("microbiology audit lacks positive-missing-result-time count")

    print(json.dumps(summary, indent=2))
    for w in warnings: print("WARN", w)
    for e in errors: print("FAIL", e)
    if errors: sys.exit(2)
    print("PASS completed run passed structural and primary diagnostic checks")


if __name__ == "__main__":
    main()
