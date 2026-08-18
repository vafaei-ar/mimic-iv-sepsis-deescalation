from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from sepsis_deescalation.specification import CANDIDATE_PS_VARS
from sepsis_deescalation.weighting_v57 import run_v57_weighting_sensitivities


def main() -> None:
    parser = argparse.ArgumentParser(description="Run v5.7 weighting sensitivities from a completed MIMIC analytic cohort.")
    parser.add_argument("run_dir", type=Path, help="Completed MIMIC run directory containing analysis_cohort_weighted.csv")
    parser.add_argument("--reps", type=int, default=1000, help="Bootstrap replicates per weighting sensitivity")
    parser.add_argument("--seed", type=int, default=20260426)
    parser.add_argument("--out-dir", type=Path, default=None)
    args = parser.parse_args()

    cohort_path = args.run_dir / "analysis_cohort_weighted.csv"
    if not cohort_path.exists():
        raise FileNotFoundError(f"Missing {cohort_path}")

    out = args.out_dir or (args.run_dir / "audits" / "v57_weighting")
    cohort = pd.read_csv(cohort_path, low_memory=False)
    result = run_v57_weighting_sensitivities(
        cohort,
        CANDIDATE_PS_VARS,
        out_dir=out,
        reps=args.reps,
        seed=args.seed,
    )

    print(f"v5.7 weighting sensitivity outputs written to: {out}")
    print("\nPoint estimates:")
    print(result["summary"].to_string(index=False))
    print("\nBootstrap confidence intervals:")
    print(result["ci"].to_string(index=False))


if __name__ == "__main__":
    main()
