#!/usr/bin/env python
from __future__ import annotations

import argparse
import logging

import pandas as pd

from sepsis_deescalation.config import load_config
from sepsis_deescalation.mimic_pipeline import run_mimic
from sepsis_deescalation.provenance import zip_run
from sepsis_deescalation.specification import CANDIDATE_PS_VARS
from sepsis_deescalation.weighting_final import run_final_weighting_sensitivities


def main() -> None:
    parser = argparse.ArgumentParser(description="Run MIMIC-IV day-3 antibiotic de-escalation analysis.")
    parser.add_argument("--config", default="config/mimic.yaml")
    parser.add_argument("--log-level", default="INFO")
    args = parser.parse_args()
    logging.basicConfig(level=getattr(logging, args.log_level.upper()), format="%(asctime)s %(levelname)s %(name)s: %(message)s")

    cfg = load_config(args.config)
    run_dir = run_mimic(args.config)

    cohort_path = run_dir / "analysis_cohort_weighted.csv"
    if cohort_path.exists():
        cohort = pd.read_csv(cohort_path, low_memory=False)
        weighting_dir = run_dir / "audits" / "final_weighting"
        run_final_weighting_sensitivities(
            cohort,
            CANDIDATE_PS_VARS,
            out_dir=weighting_dir,
            truncation_percentiles=cfg.get("weighting", {}).get(
                "truncation_percentiles", [[1.0, 99.0], [2.5, 97.5]]
            ),
            reps=int(cfg.get("bootstrap", {}).get("weighting_sensitivity_reps", 1000)),
            seed=int(cfg.get("bootstrap", {}).get("seed", 20260426)) + 700,
        )
        # run_mimic creates its archive before this post-processing step. Rebuild
        # the archive so the weighting sensitivity outputs are included.
        zip_run(run_dir)

    print(f"Completed MIMIC analysis: {run_dir}")
    print(f"ZIP: {run_dir}.zip")


if __name__ == "__main__":
    main()
