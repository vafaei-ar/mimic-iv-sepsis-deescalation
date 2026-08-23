#!/usr/bin/env python
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from pathlib import Path
from types import SimpleNamespace

import pandas as pd

from sepsis_deescalation.config import load_config
from sepsis_deescalation.runtime_optimized import (
    install_runtime_optimizations,
    primary_and_secondary_results_fast,
    progressive_adjustment_fast,
    write_timings,
)
from sepsis_deescalation.specification import CANDIDATE_PS_VARS
from sepsis_deescalation.stats import fit_stabilized_iptw
from sepsis_deescalation.weighting_final import run_final_weighting_sensitivities


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Resume/re-run statistical inference from an existing patient-level analytic cohort without rereading MIMIC."
    )
    parser.add_argument("run_dir", type=Path, help="Existing completed MIMIC run directory")
    parser.add_argument("--config", default="config/mimic.yaml")
    parser.add_argument("--mode", choices=["fast", "final"], default="fast")
    parser.add_argument("--jobs", default="auto")
    args = parser.parse_args()

    cfg = load_config(args.config)
    csv_path = args.run_dir / "analysis_cohort_weighted.csv"
    parquet_path = Path(cfg.get("execution", {}).get("cache_root", "outputs/cache/mimic")) / str(cfg.get("analysis_version", "unknown")) / "analysis_cohort_weighted.parquet"
    if parquet_path.exists():
        cohort = pd.read_parquet(parquet_path)
        source = parquet_path
    elif csv_path.exists():
        cohort = pd.read_csv(csv_path, low_memory=False)
        source = csv_path
    else:
        raise FileNotFoundError(f"No cached Parquet or {csv_path} found")

    # Configure global runtime settings without altering cohort construction.
    import sepsis_deescalation.mimic_pipeline as mimic_pipeline

    jobs = install_runtime_optimizations(mimic_pipeline, jobs=args.jobs, mode=args.mode)
    rerun_dir = args.run_dir / "inference_reruns" / f"{args.mode}_{_timestamp()}"
    tables = rerun_dir / "tables"
    diagnostics = rerun_dir / "diagnostics"
    logs = rerun_dir / "logs"
    for p in (rerun_dir, tables, diagnostics, logs):
        p.mkdir(parents=True, exist_ok=True)
    paths = SimpleNamespace(run_dir=rerun_dir, tables=tables, diagnostics=diagnostics, logs=logs)

    # Refit the point-estimate primary PS so this command does not depend on the
    # saved SW_A column being present or current.
    cohort_w, _, _ = fit_stabilized_iptw(cohort, CANDIDATE_PS_VARS)
    primary = primary_and_secondary_results_fast(cohort_w, cfg, paths)["outcomes"]
    primary.to_csv(tables / "primary_secondary_outcomes.csv", index=False)

    progressive = progressive_adjustment_fast(cohort, cfg, paths)
    progressive.to_csv(tables / "progressive_adjustment.csv", index=False)

    weighting_reps = int(cfg.get("bootstrap", {}).get("weighting_sensitivity_reps", 1000))
    if args.mode == "fast":
        weighting_reps = min(weighting_reps, int(cfg.get("execution", {}).get("fast_bootstrap_reps", 100)))
    run_final_weighting_sensitivities(
        cohort,
        CANDIDATE_PS_VARS,
        out_dir=rerun_dir / "final_weighting",
        truncation_percentiles=cfg.get("weighting", {}).get("truncation_percentiles", [[1.0, 99.0], [2.5, 97.5]]),
        reps=weighting_reps,
        seed=int(cfg.get("bootstrap", {}).get("seed", 20260426)) + 700,
        jobs=jobs,
    )
    write_timings(logs / "runtime_timings.csv")
    (rerun_dir / "README.txt").write_text(
        "Inference-only rerun from: " + str(source) + "\n"
        "This avoids MIMIC extraction/feature construction. It recomputes primary/secondary, progressive, and weighting inference.\n"
        "It does not replace a complete final MIMIC package because microbiology-membership and missing-stop-time sensitivities require upstream source data.\n",
        encoding="utf-8",
    )
    print(f"Inference-only rerun completed: {rerun_dir}")
    print(f"Mode: {args.mode}; bootstrap jobs: {jobs}")


if __name__ == "__main__":
    main()
