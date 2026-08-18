from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from sepsis_deescalation.specification import CANDIDATE_PS_VARS
from sepsis_deescalation.weighting_audit import design_matrix_summary, run_weighting_audit


def main() -> None:
    parser = argparse.ArgumentParser(description="Audit weighting/positivity for a completed MIMIC run.")
    parser.add_argument("run_dir", type=Path, help="Completed MIMIC run directory containing analysis_cohort_weighted.csv")
    parser.add_argument("--out-dir", type=Path, default=None, help="Optional output directory; default is <run_dir>/audits/weighting")
    args = parser.parse_args()

    cohort_path = args.run_dir / "analysis_cohort_weighted.csv"
    if not cohort_path.exists():
        raise FileNotFoundError(f"Missing {cohort_path}")

    out_dir = args.out_dir or (args.run_dir / "audits" / "weighting")
    out_dir.mkdir(parents=True, exist_ok=True)
    cohort = pd.read_csv(cohort_path, low_memory=False)
    result = run_weighting_audit(cohort, CANDIDATE_PS_VARS, out_dir=out_dir)

    std_path = args.run_dir / "diagnostics" / "ps_standardization.csv"
    design = pd.DataFrame()
    if std_path.exists():
        std = pd.read_csv(std_path)
        ps_vars = std["ps_variable"].dropna().astype(str).tolist()
        design = design_matrix_summary(cohort, ps_vars)
        design.to_csv(out_dir / "ps_design_matrix_summary.csv", index=False)

    print("Weighting audit written to:", out_dir)
    print(result["summary"].to_string(index=False))
    if len(design):
        print("\nPropensity-score design matrix:")
        print(design.to_string(index=False))
    if len(result["duplicates"]):
        print("\nExact duplicate covariates detected:")
        print(result["duplicates"].to_string(index=False))
    print("\nTop residual imbalances under primary IPTW:")
    print(result["top_imbalanced"].head(10).to_string(index=False))


if __name__ == "__main__":
    main()
