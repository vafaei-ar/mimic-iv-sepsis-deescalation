from __future__ import annotations

from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd

from .stats import balance_table, effective_sample_size, fit_stabilized_iptw, risks


def add_overlap_weights(df: pd.DataFrame, ps_col: str = "ps_den", treatment_col: str = "A") -> pd.DataFrame:
    d = df.copy()
    ps = pd.to_numeric(d[ps_col], errors="coerce").clip(0.001, 0.999)
    a = pd.to_numeric(d[treatment_col], errors="coerce")
    d["OW_A"] = np.where(a == 1, 1.0 - ps, ps)
    return d


def _strategy_summary(df: pd.DataFrame, vars_: Sequence[str], weight_col: str, label: str, estimand: str) -> dict:
    rt, rc, rd, rr = risks(df, "death_by_horizon", weight_col)
    bal = balance_table(df, vars_, weight_col=weight_col)
    t = df["A"] == 1
    c = df["A"] == 0
    return {
        "analysis": label,
        "estimand": estimand,
        "risk_deescalated_stopped": rt,
        "risk_continued": rc,
        "risk_difference": rd,
        "risk_ratio": rr,
        "max_post_smd": float(bal["after"].max()) if len(bal) else np.nan,
        "worst_balanced_variable": str(bal.sort_values("after", ascending=False).iloc[0]["variable"]) if len(bal) else "",
        "ess_deescalated_stopped": effective_sample_size(df.loc[t, weight_col]),
        "ess_continued": effective_sample_size(df.loc[c, weight_col]),
        "max_weight": float(pd.to_numeric(df[weight_col], errors="coerce").max()),
    }


def _bootstrap_all(
    cohort: pd.DataFrame,
    variables: Sequence[str],
    truncation_percentiles: Sequence[Sequence[float]],
    reps: int,
    seed: int,
) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    rows: list[dict] = []
    n = len(cohort)
    for rep in range(reps):
        s = cohort.iloc[rng.integers(0, n, size=n)].copy()
        if s["A"].nunique() != 2:
            continue
        try:
            w, _, _ = fit_stabilized_iptw(s, variables)
        except Exception:
            continue

        ow = add_overlap_weights(w)
        rt, rc, rd, rr = risks(ow, "death_by_horizon", "OW_A")
        rows.append({
            "rep": rep,
            "analysis": "Overlap weighting",
            "risk_treated": rt,
            "risk_control": rc,
            "risk_difference": rd,
            "risk_ratio": rr,
        })

        for low, high in truncation_percentiles:
            low = float(low)
            high = float(high)
            col = "SW_A_truncated"
            tmp = w.copy()
            lo = float(tmp["SW_A"].quantile(low / 100.0))
            hi = float(tmp["SW_A"].quantile(high / 100.0))
            tmp[col] = tmp["SW_A"].clip(lo, hi)
            rt, rc, rd, rr = risks(tmp, "death_by_horizon", col)
            rows.append({
                "rep": rep,
                "analysis": f"IPTW truncated {low:g}/{high:g}",
                "risk_treated": rt,
                "risk_control": rc,
                "risk_difference": rd,
                "risk_ratio": rr,
            })
    return pd.DataFrame(rows)


def _bootstrap_ci(boot: pd.DataFrame, requested: int) -> pd.DataFrame:
    rows: list[dict] = []
    if boot.empty:
        return pd.DataFrame(columns=["analysis", "estimand", "lower_95", "upper_95", "n_success", "n_requested"])
    for label, sub in boot.groupby("analysis", sort=False):
        for est in ["risk_difference", "risk_ratio"]:
            vals = pd.to_numeric(sub[est], errors="coerce").replace([np.inf, -np.inf], np.nan).dropna()
            rows.append({
                "analysis": label,
                "estimand": est,
                "lower_95": vals.quantile(0.025),
                "upper_95": vals.quantile(0.975),
                "n_success": len(vals),
                "n_requested": requested,
            })
    return pd.DataFrame(rows)


def run_final_weighting_sensitivities(
    cohort: pd.DataFrame,
    candidate_vars: Sequence[str],
    out_dir: str | Path,
    truncation_percentiles: Sequence[Sequence[float]] = ((1.0, 99.0), (2.5, 97.5)),
    reps: int = 1000,
    seed: int = 20260426,
) -> dict[str, pd.DataFrame]:
    """Run prespecified final weighting sensitivities.

    Primary inference remains stabilized IPTW for the ATE. Overlap weighting is
    reported as a separate overlap-population estimand. Weight truncation remains
    an ATE sensitivity and is not promoted solely because it changes the point estimate.
    """
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    base, _, _ = fit_stabilized_iptw(cohort, candidate_vars)
    rows = [_strategy_summary(base, candidate_vars, "SW_A", "Primary stabilized IPTW", "ATE")]

    overlap = add_overlap_weights(base)
    rows.append(_strategy_summary(overlap, candidate_vars, "OW_A", "Overlap weighting", "overlap population"))
    balance_table(overlap, candidate_vars, weight_col="OW_A").to_csv(out / "balance_overlap.csv", index=False)

    for low, high in truncation_percentiles:
        low = float(low)
        high = float(high)
        tmp = base.copy()
        col = f"SW_A_trunc_{low:g}_{high:g}"
        tmp[col] = tmp["SW_A"].clip(tmp["SW_A"].quantile(low / 100.0), tmp["SW_A"].quantile(high / 100.0))
        rows.append(_strategy_summary(tmp, candidate_vars, col, f"IPTW truncated {low:g}/{high:g}", "ATE, truncated weights"))
        balance_table(tmp, candidate_vars, weight_col=col).to_csv(out / f"balance_truncated_{low:g}_{high:g}.csv", index=False)

    summary = pd.DataFrame(rows)
    summary.to_csv(out / "final_weighting_point_estimates.csv", index=False)

    boot = _bootstrap_all(cohort, candidate_vars, truncation_percentiles, reps, seed)
    boot.to_csv(out / "final_weighting_bootstrap_replicates.csv", index=False)
    ci = _bootstrap_ci(boot, reps)
    ci.to_csv(out / "final_weighting_bootstrap_ci.csv", index=False)
    return {"summary": summary, "bootstrap": boot, "ci": ci}
