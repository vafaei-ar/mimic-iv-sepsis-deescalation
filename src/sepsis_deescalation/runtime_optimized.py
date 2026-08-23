from __future__ import annotations

import logging
import time
from pathlib import Path

import numpy as np
import pandas as pd

from .fast_bootstrap import OutcomeSpec, bootstrap_ci_long, bootstrap_multi_outcome_iptw, resolve_jobs
from .specification import CANDIDATE_PS_VARS, PROGRESSIVE_MODELS
from .stats import balance_table, fit_stabilized_iptw, mean_difference, risks

LOG = logging.getLogger(__name__)
_RUNTIME_JOBS: int = 1
_RUNTIME_MODE: str = "final"
_TIMINGS: list[dict] = []


def _reps(cfg: dict, key: str) -> int:
    requested = int(cfg["bootstrap"][key])
    if _RUNTIME_MODE == "fast":
        return min(requested, int(cfg.get("execution", {}).get("fast_bootstrap_reps", 100)))
    return requested


def _record(stage: str, started: float, **extra) -> None:
    row = {"stage": stage, "seconds": time.perf_counter() - started, "mode": _RUNTIME_MODE, "jobs": _RUNTIME_JOBS}
    row.update(extra)
    _TIMINGS.append(row)
    LOG.info("Timing %-38s %.1f s", stage, row["seconds"])


def reset_timings() -> None:
    _TIMINGS.clear()


def write_timings(path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(_TIMINGS).to_csv(path, index=False)


def _results_row(name: str, values: tuple, kind: str = "risk") -> pd.DataFrame:
    if kind == "risk":
        rt, rc, rd, rr = values
        return pd.DataFrame([{"analysis": name, "risk_deescalated_stopped": rt, "risk_continued": rc, "risk_difference": rd, "risk_ratio": rr}])
    mt, mc, md = values
    return pd.DataFrame([{"analysis": name, "mean_deescalated_stopped": mt, "mean_continued": mc, "mean_difference": md}])


def _write(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)


def primary_and_secondary_results_fast(cohort_w: pd.DataFrame, cfg: dict, paths) -> dict[str, pd.DataFrame]:
    started = time.perf_counter()
    reps = _reps(cfg, "primary_reps")
    seed = int(cfg["bootstrap"]["seed"])
    outcomes = [
        OutcomeSpec("30-day mortality", "death_by_horizon", "risk"),
        OutcomeSpec("hospital-free days", "hospital_free_days", "mean"),
        OutcomeSpec("antibiotic-free days", "antibiotic_free_days", "mean"),
        OutcomeSpec("normalized systemic antibiotic exposure", "normalized_antibiotic_exposure_30d", "mean"),
        OutcomeSpec("normalized broad-spectrum exposure", "normalized_broad_antibiotic_exposure_30d", "mean"),
        OutcomeSpec("late recurrent/persistent antibiotic-course use", "late_recurrent_or_persistent_abx_course_30d", "risk"),
    ]

    # One PS fit per bootstrap replicate is shared by all six outcomes.
    boot, boot_diag = bootstrap_multi_outcome_iptw(
        cohort_w,
        CANDIDATE_PS_VARS,
        outcomes,
        reps,
        seed,
        jobs=_RUNTIME_JOBS,
    )
    _write(boot_diag, paths.diagnostics / "bootstrap_primary_shared_fit_diagnostics.csv")
    ci_all = bootstrap_ci_long(boot, reps)

    rows = []
    for spec in outcomes:
        point = risks(cohort_w, spec.column, "SW_A") if spec.kind == "risk" else mean_difference(cohort_w, spec.column, "SW_A")
        result = _results_row(spec.label, point, spec.kind)
        result["kind"] = spec.kind
        b = boot.loc[boot["analysis"] == spec.label].copy() if len(boot) else pd.DataFrame()
        c = ci_all.loc[ci_all["analysis"] == spec.label].copy() if len(ci_all) else pd.DataFrame()
        _write(b, paths.diagnostics / f"bootstrap_{spec.column}.csv")
        _write(c, paths.tables / f"ci_{spec.column}.csv")
        primary_est = "risk_difference" if spec.kind == "risk" else "mean_difference"
        hit = c.loc[c["estimand"] == primary_est] if len(c) else pd.DataFrame()
        if len(hit):
            result["lower_95"] = float(hit.iloc[0]["lower_95"])
            result["upper_95"] = float(hit.iloc[0]["upper_95"])
            result["bootstrap_success"] = int(hit.iloc[0]["n_success"])
        rows.append(result)
    out = {"outcomes": pd.concat(rows, ignore_index=True)}
    _record("bootstrap_primary_secondary_shared", started, reps=reps, ps_fits=reps)
    return out


def progressive_adjustment_fast(cohort: pd.DataFrame, cfg: dict, paths) -> pd.DataFrame:
    started_all = time.perf_counter()
    reps = _reps(cfg, "progressive_reps")
    seed = int(cfg["bootstrap"]["seed"])
    rows = []
    for i, spec in enumerate(PROGRESSIVE_MODELS):
        started = time.perf_counter()
        w, _, diag = fit_stabilized_iptw(cohort, spec["vars"])
        rt, rc, rd, rr = risks(w, "death_by_horizon", "SW_A")
        bal = balance_table(w, spec["vars"])
        outcomes = [OutcomeSpec(spec["model"], "death_by_horizon", "risk")]
        boot, boot_diag = bootstrap_multi_outcome_iptw(
            cohort,
            spec["vars"],
            outcomes,
            reps,
            seed + 100 + i,
            jobs=_RUNTIME_JOBS,
        )
        _write(boot, paths.diagnostics / f"progressive_bootstrap_M{i+1}.csv")
        _write(boot_diag, paths.diagnostics / f"progressive_bootstrap_M{i+1}_diagnostics.csv")
        ci = bootstrap_ci_long(boot, reps)
        rd_ci = ci.loc[ci["estimand"] == "risk_difference"] if len(ci) else pd.DataFrame()
        rr_ci = ci.loc[ci["estimand"] == "risk_ratio"] if len(ci) else pd.DataFrame()
        rows.append({
            "model": spec["model"],
            "n_vars_requested": len(spec["vars"]),
            "n_vars_used": len(diag["used_vars"]),
            "risk_deescalated_stopped": rt,
            "risk_continued": rc,
            "risk_difference": rd,
            "rd_lower_95": float(rd_ci.iloc[0]["lower_95"]) if len(rd_ci) else np.nan,
            "rd_upper_95": float(rd_ci.iloc[0]["upper_95"]) if len(rd_ci) else np.nan,
            "risk_ratio": rr,
            "rr_lower_95": float(rr_ci.iloc[0]["lower_95"]) if len(rr_ci) else np.nan,
            "rr_upper_95": float(rr_ci.iloc[0]["upper_95"]) if len(rr_ci) else np.nan,
            "max_post_smd": float(bal["after"].max()) if len(bal) else np.nan,
        })
        _record(f"progressive_M{i+1}", started, reps=reps, ps_fits=reps)
    _record("progressive_all", started_all, reps_per_model=reps, ps_fits=reps * len(PROGRESSIVE_MODELS))
    return pd.DataFrame(rows)


def simple_sensitivity_fast(
    df: pd.DataFrame,
    label: str,
    cfg: dict,
    reps: int | None = None,
    seed_offset: int = 0,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    started = time.perf_counter()
    if len(df) < 100 or df["A"].nunique() != 2:
        return pd.DataFrame([{"analysis": label, "n": len(df), "status": "not_run"}]), pd.DataFrame()
    w, _, _ = fit_stabilized_iptw(df, CANDIDATE_PS_VARS)
    rt, rc, rd, rr = risks(w, "death_by_horizon", "SW_A")
    out = pd.DataFrame([{
        "analysis": label,
        "n": len(w),
        "treated": int(w["A"].sum()),
        "control": int((w["A"] == 0).sum()),
        "risk_treated": rt,
        "risk_control": rc,
        "risk_difference": rd,
        "risk_ratio": rr,
        "status": "ok",
    }])
    reps = _reps(cfg, "primary_reps") if reps is None else (min(int(reps), int(cfg.get("execution", {}).get("fast_bootstrap_reps", 100))) if _RUNTIME_MODE == "fast" else int(reps))
    outcomes = [OutcomeSpec(label, "death_by_horizon", "risk")]
    boot, _ = bootstrap_multi_outcome_iptw(
        df,
        CANDIDATE_PS_VARS,
        outcomes,
        reps,
        int(cfg["bootstrap"]["seed"]) + seed_offset,
        jobs=_RUNTIME_JOBS,
    )
    ci = bootstrap_ci_long(boot, reps)
    rd_ci = ci.loc[ci["estimand"] == "risk_difference"] if len(ci) else pd.DataFrame()
    rr_ci = ci.loc[ci["estimand"] == "risk_ratio"] if len(ci) else pd.DataFrame()
    if len(rd_ci):
        out["rd_lower_95"] = float(rd_ci.iloc[0]["lower_95"])
        out["rd_upper_95"] = float(rd_ci.iloc[0]["upper_95"])
    if len(rr_ci):
        out["rr_lower_95"] = float(rr_ci.iloc[0]["lower_95"])
        out["rr_upper_95"] = float(rr_ci.iloc[0]["upper_95"])
    _record(f"sensitivity:{label[:45]}", started, reps=reps, ps_fits=reps)
    return out, boot


def install_runtime_optimizations(mimic_pipeline_module, jobs: int | str | None = "auto", mode: str = "final") -> int:
    """Install speed-only replacements for the expensive bootstrap helpers.

    The cohort construction, phenotype definitions, point-estimate PS model and
    outcome definitions remain in mimic_pipeline unchanged. Only bootstrap
    execution is replaced by shared-fit/matrix/parallel implementations.
    """
    global _RUNTIME_JOBS, _RUNTIME_MODE
    if mode not in {"fast", "final"}:
        raise ValueError("mode must be 'fast' or 'final'")
    _RUNTIME_MODE = mode
    _RUNTIME_JOBS = resolve_jobs(jobs)
    reset_timings()
    mimic_pipeline_module._primary_and_secondary_results = primary_and_secondary_results_fast
    mimic_pipeline_module._progressive_adjustment = progressive_adjustment_fast
    mimic_pipeline_module._simple_sensitivity = simple_sensitivity_fast
    LOG.info("Installed optimized inference runtime: mode=%s jobs=%s", _RUNTIME_MODE, _RUNTIME_JOBS)
    return _RUNTIME_JOBS
