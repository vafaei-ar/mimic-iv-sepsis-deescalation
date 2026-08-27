#!/usr/bin/env python3
"""Bootstrap inference for the two prespecified PSU robustness analyses.

This script reuses the already-frozen robustness variants and adds uncertainty
estimation only. It does not introduce any new exposure, cohort, covariate, outcome,
or weighting definitions.

For each variant, the propensity score is refit inside each of 1,000 encounter-level
nonparametric bootstrap samples. This matters because keeping the original fitted
weights fixed would understate uncertainty in a propensity-score analysis.

Variants
--------
1. Strict cohort with MED_ADMIN day-3 exposure reclassification.
2. Lenient 96-hour landmark with PRESCRIBING exposure.

Only aggregate confidence intervals and diagnostics are written.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

import run_psu_point_estimates as pe
import run_psu_prespecified_robustness as rob
from run_psu_bootstrap_inference import run_bootstrap

BOOT_INJECT = r'''
    run_bootstrap(df, binary, continuous, effect_outcomes, args.output_dir)
'''


def run_variant_bootstrap(source: str, variant: str, data_root: Path, outdir: Path) -> dict:
    s = source
    if variant == "medadmin_exposure":
        if s.count(rob.DF_MARKER) != 1:
            raise RuntimeError("Frozen PS dataframe marker missing or ambiguous")
        s = s.replace(rob.DF_MARKER, rob.DF_MARKER + rob.MEDADMIN_OVERRIDE + "\n", 1)
    elif variant == "lenient_landmark":
        if s.count(rob.STRICT_COHORT) != 1:
            raise RuntimeError("Frozen strict landmark clause missing or ambiguous")
        s = s.replace(rob.STRICT_COHORT, rob.LENIENT_COHORT, 1)
    else:
        raise ValueError(variant)

    if s.count(rob.BAL_MARKER) != 1:
        raise RuntimeError("Frozen PS balance marker missing or ambiguous")
    # Fail closed if the frozen source no longer has the expected insertion point.
    s = s.replace(rob.BAL_MARKER, pe.INJECT + "\n" + BOOT_INJECT + "\n" + rob.BAL_MARKER, 1)

    outdir.mkdir(parents=True, exist_ok=True)
    old_argv = sys.argv[:]
    try:
        sys.argv = ["audit_psu_ps_balance.py", str(data_root), "--output-dir", str(outdir)]
        ns = {
            "__name__": "__main__",
            "__file__": str(Path(__file__).resolve().parent / "audit_psu_ps_balance.py"),
            "run_bootstrap": run_bootstrap,
        }
        exec(compile(s, ns["__file__"], "exec"), ns, ns)
    finally:
        sys.argv = old_argv

    pts = pd.read_csv(outdir / "point_estimates.csv")
    death = pts[(pts["method"] == "stabilized_ate_iptw") & (pts["outcome"] == "death_30d")].iloc[0]
    ci = pd.read_csv(outdir / "bootstrap_ci.csv")
    rd = ci[
        (ci["method"] == "stabilized_ate_iptw")
        & (ci["outcome"] == "death_30d")
        & (ci["estimand"] == "difference_A1_minus_A0")
    ].iloc[0]
    rr = ci[
        (ci["method"] == "stabilized_ate_iptw")
        & (ci["outcome"] == "death_30d")
        & (ci["estimand"] == "risk_ratio_A1_over_A0")
    ].iloc[0]
    diag = json.loads((outdir / "bootstrap_diagnostics.json").read_text())
    summary = json.loads((outdir / "effect_summary.json").read_text())
    return {
        "sensitivity": variant,
        "cohort_n": int(summary["strict_cohort_n"]),
        "deescalated_n": int(summary["deescalated_n"]),
        "continued_n": int(summary["continued_n"]),
        "death_risk_difference": float(death["difference_A1_minus_A0"]),
        "death_rd_lower_95": float(rd["lower_95"]),
        "death_rd_upper_95": float(rd["upper_95"]),
        "death_risk_ratio": float(death["risk_ratio_A1_over_A0"]),
        "death_rr_lower_95": float(rr["lower_95"]),
        "death_rr_upper_95": float(rr["upper_95"]),
        "bootstrap_successful": int(diag["n_successful_replicates"]),
        "bootstrap_failed": int(diag["n_failed_replicates"]),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("data_root", type=Path)
    ap.add_argument("--output-dir", type=Path, required=True)
    args = ap.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    source_path = Path(__file__).resolve().parent / "audit_psu_ps_balance.py"
    source = source_path.read_text(encoding="utf-8")

    rows = [
        run_variant_bootstrap(source, "medadmin_exposure", args.data_root, args.output_dir / "medadmin_exposure"),
        run_variant_bootstrap(source, "lenient_landmark", args.data_root, args.output_dir / "lenient_landmark"),
    ]
    pd.DataFrame(rows).to_csv(args.output_dir / "robustness_bootstrap_summary.csv", index=False)
    meta = {
        "privacy_mode": "aggregate_only",
        "bootstrap_replicates_per_sensitivity": 1000,
        "sensitivities": ["medadmin_exposure", "lenient_landmark"],
        "guardrail": "Only prespecified PSU robustness variants are bootstrapped. Exposure variant/landmark rule changes are exactly those frozen before outcome inspection; covariate, PS, and outcome definitions are otherwise unchanged.",
    }
    (args.output_dir / "robustness_bootstrap_metadata.json").write_text(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()
