from __future__ import annotations

import json
import shutil
from pathlib import Path
from typing import Any

import pandas as pd


def save_inference_cache(
    cohort: pd.DataFrame,
    run_dir: str | Path,
    analysis_version: str,
    metadata: dict[str, Any] | None = None,
    cache_root: str | Path = "outputs/cache/mimic",
) -> Path:
    """Persist a local patient-level Parquet checkpoint for inference reruns.

    The cache lives under outputs/, is gitignored, and must never be committed.
    It intentionally requires explicit reuse; the main MIMIC extraction never
    silently trusts an old cache.
    """
    run_dir = Path(run_dir)
    cache_dir = Path(cache_root) / analysis_version
    cache_dir.mkdir(parents=True, exist_ok=True)
    cohort_path = cache_dir / "analysis_cohort_weighted.parquet"
    cohort.to_parquet(cohort_path, index=False)
    manifest = {
        "analysis_version": analysis_version,
        "source_run_dir": str(run_dir),
        "n_rows": int(len(cohort)),
    }
    if metadata:
        manifest.update(metadata)
    (cache_dir / "cache_manifest.json").write_text(json.dumps(manifest, indent=2, default=str), encoding="utf-8")
    return cohort_path


def copy_cached_inference_artifacts(source_run_dir: str | Path, target_run_dir: str | Path) -> None:
    """Copy non-recomputed deterministic artifacts into an inference-only rerun."""
    source_run_dir = Path(source_run_dir)
    target_run_dir = Path(target_run_dir)
    for rel in [
        "tables/cohort_flow.csv",
        "tables/stop_time_assumption_sensitivity.csv",
        "audits/microbiology_counts.csv",
        "audits/antibiotic_definition.csv",
    ]:
        src = source_run_dir / rel
        dst = target_run_dir / rel
        if src.exists():
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, dst)
