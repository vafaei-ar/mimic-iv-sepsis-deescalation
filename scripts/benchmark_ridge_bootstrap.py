#!/usr/bin/env python
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from sepsis_deescalation.ridge_benchmark import benchmark_ridge_vs_auto
from sepsis_deescalation.specification import CANDIDATE_PS_VARS


def main() -> None:
    parser = argparse.ArgumentParser(description="Benchmark auto GLM/fallback versus direct ridge PS fitting on identical bootstrap samples.")
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("--reps", type=int, default=20)
    parser.add_argument("--seed", type=int, default=20260426)
    args = parser.parse_args()

    cohort_path = args.run_dir / "analysis_cohort_weighted.csv"
    if not cohort_path.exists():
        raise FileNotFoundError(f"Missing {cohort_path}")
    cohort = pd.read_csv(cohort_path, low_memory=False)
    out = args.run_dir / "audits" / "ridge_speed_benchmark"
    out.mkdir(parents=True, exist_ok=True)
    detail = benchmark_ridge_vs_auto(cohort, CANDIDATE_PS_VARS, reps=args.reps, seed=args.seed)
    detail.to_csv(out / "ridge_vs_auto_detail.csv", index=False)
    summary = pd.DataFrame([{
        "reps": len(detail),
        "median_auto_seconds": detail["auto_seconds"].median(),
        "median_ridge_seconds": detail["ridge_seconds"].median(),
        "median_speedup": detail["speedup"].median(),
        "max_abs_rd_difference": detail["rd_difference"].abs().max(),
        "median_abs_rd_difference": detail["rd_difference"].abs().median(),
        "max_abs_rr_difference": detail["rr_difference"].abs().max(),
        "max_abs_ps_difference": detail["max_abs_ps_difference"].max(),
        "max_abs_weight_difference": detail["max_abs_weight_difference"].max(),
        "n_auto_glm": int((detail["auto_method"] == "glm").sum()),
        "n_auto_ridge_fallback": int((detail["auto_method"] == "ridge_fallback").sum()),
    }])
    summary.to_csv(out / "ridge_vs_auto_summary.csv", index=False)
    print(summary.to_string(index=False))
    print(f"Outputs: {out}")


if __name__ == "__main__":
    main()
