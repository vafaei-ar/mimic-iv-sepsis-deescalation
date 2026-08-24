#!/usr/bin/env python
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from sepsis_deescalation.capped_glm_benchmark import benchmark_capped_glm
from sepsis_deescalation.specification import CANDIDATE_PS_VARS


def main() -> None:
    parser = argparse.ArgumentParser(description="Benchmark capped GLM attempts against the historical 200-iteration fallback strategy.")
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("--reps", type=int, default=20)
    parser.add_argument("--seed", type=int, default=20260426)
    parser.add_argument("--caps", nargs="+", type=int, default=[10, 20, 30, 50])
    args = parser.parse_args()

    cohort_path = args.run_dir / "analysis_cohort_weighted.csv"
    if not cohort_path.exists():
        raise FileNotFoundError(f"Missing {cohort_path}")
    cohort = pd.read_csv(cohort_path, low_memory=False)

    detail, summary = benchmark_capped_glm(
        cohort,
        CANDIDATE_PS_VARS,
        reps=args.reps,
        seed=args.seed,
        caps=args.caps,
    )
    out = args.run_dir / "audits" / "capped_glm_benchmark"
    out.mkdir(parents=True, exist_ok=True)
    detail.to_csv(out / "capped_glm_detail.csv", index=False)
    summary.to_csv(out / "capped_glm_summary.csv", index=False)

    print(summary.to_string(index=False))
    print(f"Outputs: {out}")


if __name__ == "__main__":
    main()
