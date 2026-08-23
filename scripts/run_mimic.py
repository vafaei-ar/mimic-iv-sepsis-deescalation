#!/usr/bin/env python
from __future__ import annotations

import argparse
import logging
import time

import pandas as pd

import sepsis_deescalation.mimic_pipeline as mimic_pipeline
from sepsis_deescalation.cache import save_inference_cache
from sepsis_deescalation.config import load_config
from sepsis_deescalation.provenance import zip_run
from sepsis_deescalation.runtime_optimized import install_runtime_optimizations, write_timings
from sepsis_deescalation.specification import CANDIDATE_PS_VARS
from sepsis_deescalation.weighting_final import run_final_weighting_sensitivities


def main() -> None:
    parser = argparse.ArgumentParser(description="Run MIMIC-IV day-3 antibiotic de-escalation analysis.")
    parser.add_argument("--config", default="config/mimic.yaml")
    parser.add_argument("--log-level", default="INFO")
    parser.add_argument(
        "--mode",
        choices=["fast", "final"],
        default="final",
        help="fast uses a reduced bootstrap count for development; final uses configured publication counts",
    )
    parser.add_argument(
        "--jobs",
        default="auto",
        help="Bootstrap worker processes (integer or 'auto'; auto uses up to 8 cores)",
    )
    parser.add_argument(
        "--no-cache",
        action="store_true",
        help="Do not write the local Parquet inference checkpoint under outputs/cache/",
    )
    args = parser.parse_args()
    logging.basicConfig(level=getattr(logging, args.log_level.upper()), format="%(asctime)s %(levelname)s %(name)s: %(message)s")

    cfg = load_config(args.config)
    jobs = install_runtime_optimizations(mimic_pipeline, jobs=args.jobs, mode=args.mode)

    started = time.perf_counter()
    run_dir = mimic_pipeline.run_mimic(args.config)
    core_seconds = time.perf_counter() - started

    cohort_path = run_dir / "analysis_cohort_weighted.csv"
    if cohort_path.exists():
        cohort = pd.read_csv(cohort_path, low_memory=False)

        if not args.no_cache:
            cache_path = save_inference_cache(
                cohort,
                run_dir,
                analysis_version=str(cfg.get("analysis_version", "unknown")),
                metadata={"mode": args.mode, "jobs": jobs},
            )
            logging.info("Inference cache: %s", cache_path)

        weighting_dir = run_dir / "audits" / "final_weighting"
        weighting_reps = int(cfg.get("bootstrap", {}).get("weighting_sensitivity_reps", 1000))
        if args.mode == "fast":
            weighting_reps = min(weighting_reps, int(cfg.get("execution", {}).get("fast_bootstrap_reps", 100)))
        weighting_started = time.perf_counter()
        run_final_weighting_sensitivities(
            cohort,
            CANDIDATE_PS_VARS,
            out_dir=weighting_dir,
            truncation_percentiles=cfg.get("weighting", {}).get(
                "truncation_percentiles", [[1.0, 99.0], [2.5, 97.5]]
            ),
            reps=weighting_reps,
            seed=int(cfg.get("bootstrap", {}).get("seed", 20260426)) + 700,
            jobs=jobs,
        )
        weighting_seconds = time.perf_counter() - weighting_started

        timings = run_dir / "logs" / "runtime_timings.csv"
        write_timings(timings)
        # Add the two outer stages without coupling timing code to cohort-building internals.
        outer = pd.DataFrame([
            {"stage": "core_pipeline_total", "seconds": core_seconds, "mode": args.mode, "jobs": jobs},
            {"stage": "final_weighting_total", "seconds": weighting_seconds, "mode": args.mode, "jobs": jobs},
        ])
        if timings.exists():
            inner = pd.read_csv(timings)
            pd.concat([outer, inner], ignore_index=True, sort=False).to_csv(timings, index=False)
        else:
            outer.to_csv(timings, index=False)

        # run_mimic creates its archive before post-processing. Rebuild so the
        # weighting, cache metadata/timing diagnostics are reflected in the ZIP.
        zip_run(run_dir)

    print(f"Completed MIMIC analysis: {run_dir}")
    print(f"ZIP: {run_dir}.zip")
    print(f"Mode: {args.mode}; bootstrap jobs: {jobs}")


if __name__ == "__main__":
    main()
