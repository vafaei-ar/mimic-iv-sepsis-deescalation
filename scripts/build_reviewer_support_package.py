#!/usr/bin/env python3
"""Build reviewer-support aggregates and revised publication figures.

This task deliberately preserves the frozen primary MIMIC-IV and PSU analyses. It
uses the final vital-corrected MIMIC analytic cohort only for targeted reviewer
support analyses that were not part of the primary estimand: treatment-subtype
composition/descriptive outcomes and subtype-vs-continuation mortality
sensitivities. All exported artifacts are aggregate/sanitized.
"""
from __future__ import annotations

from pathlib import Path
import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from publication_figure_common import (
    apply_publication_secondary_overrides,
    prepare_progressive_mortality,
    pretty_label,
    shared_histogram_edges,
)
from sepsis_deescalation.specification import CANDIDATE_PS_VARS
from sepsis_deescalation.stats import (
    balance_table,
    bootstrap_iptw_ci,
    effective_sample_size,
    fit_stabilized_iptw,
    risks,
)

BASE_RUN = Path("outputs/mimic/mimic_iv_v5_7_final_20260820T003506Z")
COHORT = BASE_RUN / "audits/vital_repair/analysis_cohort_vital_corrected.csv"
FLOW = BASE_RUN / "tables/cohort_flow.csv"
HARM = Path("outputs/publication_integration/harmonized")
OUT = Path("outputs/publication_integration/reviewer_support")

EXPECTED_N = 9589
EXPECTED_DEESC = 1863
EXPECTED_CONT = 7726
EXPECTED_PRIMARY_RD = 0.0083758268262657


def savefig(fig: plt.Figure, stem: str) -> None:
    fig.tight_layout()
    fig.savefig(OUT / f"{stem}.png", dpi=400, bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)


def validate_primary(d: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    if len(d) != EXPECTED_N:
        raise RuntimeError(f"Corrected cohort n={len(d)}, expected {EXPECTED_N}")
    if int((d["A"] == 1).sum()) != EXPECTED_DEESC or int((d["A"] == 0).sum()) != EXPECTED_CONT:
        raise RuntimeError("Corrected treatment counts do not match the frozen publication cohort")
    w, _, _ = fit_stabilized_iptw(d, CANDIDATE_PS_VARS)
    _, _, rd, _ = risks(w, "death_by_horizon", "SW_A")
    if abs(rd - EXPECTED_PRIMARY_RD) > 1e-4:
        raise RuntimeError(f"Corrected primary RD parity failed: observed {rd}, expected {EXPECTED_PRIMARY_RD}")
    bal = balance_table(w, CANDIDATE_PS_VARS)
    return w, bal


def build_flow() -> pd.DataFrame:
    f = pd.read_csv(FLOW)
    count_col = next((c for c in f.columns if c.lower() in {"n", "count", "cohort_n"}), None)
    stage_col = next((c for c in f.columns if c != count_col), None)
    if count_col is None or stage_col is None:
        raise RuntimeError(f"Unexpected cohort_flow columns: {list(f.columns)}")
    f = f[[stage_col, count_col]].rename(columns={stage_col: "stage", count_col: "n"})
    f["removed_since_prior"] = np.nan
    prev = None
    for i, row in f.iterrows():
        n = int(row["n"])
        if prev is not None and n <= prev:
            f.loc[i, "removed_since_prior"] = prev - n
        prev = n
    f["percent_of_first_stage"] = 100 * f["n"] / float(f.iloc[0]["n"])
    f.to_csv(OUT / "mimic_cohort_flow_reviewer.csv", index=False)
    return f


def subtype_descriptives(d: pd.DataFrame) -> pd.DataFrame:
    outcomes = [
        "death_by_horizon",
        "hospital_free_days",
        "antibiotic_free_days",
        "normalized_antibiotic_exposure_30d",
        "normalized_broad_antibiotic_exposure_30d",
    ]
    rows = []
    for typ, g in d.groupby("deescalation_type", dropna=False):
        row = {
            "deescalation_type": typ,
            "n": len(g),
            "percent_of_total": 100 * len(g) / len(d),
            "percent_within_deescalated": 100 * len(g) / EXPECTED_DEESC if typ != "continued_broad" else np.nan,
        }
        for o in outcomes:
            if o in g:
                row[f"{o}_mean"] = pd.to_numeric(g[o], errors="coerce").mean()
        rows.append(row)
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "mimic_deescalation_subtype_descriptives.csv", index=False)
    return out


def subtype_sensitivity(d: pd.DataFrame, reps: int = 1000) -> pd.DataFrame:
    contrasts = [
        ("narrowed_or_non_broad_only", "Narrowed/non-broad only vs continued broad"),
        ("stopped_all_observed_systemic_antibiotics", "Stopped all observed systemic antibiotics vs continued broad"),
    ]
    rows = []
    for i, (typ, label) in enumerate(contrasts):
        s = d.loc[d["deescalation_type"].isin([typ, "continued_broad"])].copy()
        s["A"] = (s["deescalation_type"] == typ).astype(int)
        w, _, _ = fit_stabilized_iptw(s, CANDIDATE_PS_VARS)
        rt, rc, rd, rr = risks(w, "death_by_horizon", "SW_A")
        bal = balance_table(w, CANDIDATE_PS_VARS)
        ci, _ = bootstrap_iptw_ci(
            s,
            CANDIDATE_PS_VARS,
            "death_by_horizon",
            "risk",
            reps,
            20260827 + i,
        )
        rdrow = ci.loc[ci["estimand"] == "risk_difference"].iloc[0]
        rrrow = ci.loc[ci["estimand"] == "risk_ratio"].iloc[0]
        t = w["A"] == 1
        c = ~t
        rows.append({
            "analysis": label,
            "n_total": len(w),
            "n_subtype": int(t.sum()),
            "n_continued": int(c.sum()),
            "weighted_risk_subtype": rt,
            "weighted_risk_continued": rc,
            "risk_difference": rd,
            "rd_ci95_low": float(rdrow["lower_95"]),
            "rd_ci95_high": float(rdrow["upper_95"]),
            "risk_ratio": rr,
            "rr_ci95_low": float(rrrow["lower_95"]),
            "rr_ci95_high": float(rrrow["upper_95"]),
            "bootstrap_successes": int(rdrow["n_success"]),
            "max_post_smd": float(bal["after"].max()),
            "ess_subtype": effective_sample_size(w.loc[t, "SW_A"]),
            "ess_continued": effective_sample_size(w.loc[c, "SW_A"]),
            "max_weight": float(w["SW_A"].max()),
            "interpretation": "Exploratory subtype sensitivity; not a replacement for the frozen binary primary estimand.",
        })
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "mimic_deescalation_subtype_adjusted_mortality.csv", index=False)
    return out


def diagnostics(w: pd.DataFrame, bal: pd.DataFrame) -> None:
    bal.to_csv(OUT / "mimic_primary_balance_before_after.csv", index=False)
    qrows = []
    for a, label in [(1, "deescalated/stopped"), (0, "continued broad")]:
        g = w.loc[w["A"] == a]
        for var in ["ps_den", "SW_A"]:
            x = pd.to_numeric(g[var], errors="coerce").dropna()
            qrows.append({
                "group": label,
                "variable": var,
                "n": len(x),
                "min": x.min(),
                "p01": x.quantile(.01),
                "p05": x.quantile(.05),
                "p25": x.quantile(.25),
                "median": x.quantile(.50),
                "p75": x.quantile(.75),
                "p95": x.quantile(.95),
                "p99": x.quantile(.99),
                "max": x.max(),
            })
    pd.DataFrame(qrows).to_csv(OUT / "mimic_ps_weight_quantiles.csv", index=False)

    top = bal.sort_values("before", ascending=False).head(35).sort_values("before")
    fig, ax = plt.subplots(figsize=(7.2, 8.2))
    y = np.arange(len(top))
    ax.scatter(top["before"], y, marker="o", label="Before weighting")
    ax.scatter(top["after"], y, marker="s", label="After weighting")
    ax.axvline(0.1, linestyle="--", linewidth=1)
    ax.set_yticks(y)
    ax.set_yticklabels([pretty_label(v) for v in top["variable"]], fontsize=7.5)
    ax.set_xlabel("Absolute standardized mean difference")
    ax.set_title("MIMIC-IV primary covariate balance")
    ax.legend(frameon=False)
    savefig(fig, "ESM_Fig1_mimic_balance_love")

    ps_edges = np.linspace(0, 1, 31)
    weight_edges = shared_histogram_edges(
        [w.loc[w["A"] == 1, "SW_A"], w.loc[w["A"] == 0, "SW_A"]], bins=30
    )
    fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.6))
    for a, label in [(1, "De-escalated/stopped"), (0, "Continued broad")]:
        axes[0].hist(
            w.loc[w["A"] == a, "ps_den"].to_numpy(float),
            bins=ps_edges,
            histtype="step",
            linewidth=1.5,
            label=label,
        )
        axes[1].hist(
            w.loc[w["A"] == a, "SW_A"].to_numpy(float),
            bins=weight_edges,
            histtype="step",
            linewidth=1.5,
            label=label,
        )
    axes[0].set_xlabel("Estimated propensity for de-escalation")
    axes[0].set_ylabel("Admissions")
    axes[0].set_title("Propensity-score overlap")
    axes[1].set_xlabel("Stabilized IPTW")
    axes[1].set_ylabel("Admissions")
    axes[1].set_title("Primary weight distribution")
    axes[1].set_xlim(left=0)
    for ax in axes:
        ax.legend(frameon=False, fontsize=8)
    savefig(fig, "ESM_Fig2_mimic_ps_weights")


def revised_main_figures() -> None:
    mort = pd.read_csv(HARM / "harmonized_mortality_results.csv")
    sec = apply_publication_secondary_overrides(
        pd.read_csv(HARM / "harmonized_secondary_outcomes.csv")
    )
    prog = pd.read_csv(HARM / "mimic_progressive_adjustment.csv")

    fig, ax = plt.subplots(figsize=(9.2, 2.6))
    ax.set_xlim(-4, 125); ax.set_ylim(-1.0, 1.0); ax.axis("off")
    ax.annotate("", xy=(120, 0), xytext=(0, 0), arrowprops=dict(arrowstyle="->", linewidth=2))
    for x, top, bottom in [
        (0, "t0", "First qualifying broad-spectrum antibiotic"),
        (72, "72 h", "Treatment decision"),
        (96, "96 h", "Landmark; follow-up starts"),
        (120, "30 d", "Post-landmark outcome horizon"),
    ]:
        ax.plot([x, x], [-0.12, 0.12], linewidth=1.7)
        ax.text(x, 0.35, top, ha="center", va="bottom", fontsize=10, fontweight="bold")
        ax.text(x, 0.58, bottom, ha="center", va="bottom", fontsize=8.5)
    ax.plot([0, 72], [-0.36, -0.36], linewidth=5, alpha=.45)
    ax.text(36, -0.58, "Pre-decision covariates", ha="center", fontsize=9)
    ax.plot([72, 96], [-0.36, -0.36], linewidth=5, alpha=.75)
    ax.text(84, -0.58, "Treatment classification", ha="center", fontsize=9)
    ax.plot([96, 120], [-0.36, -0.36], linewidth=5)
    ax.text(108, -0.58, "Outcome follow-up", ha="center", fontsize=9)
    savefig(fig, "Fig1_target_trial_timeline_revised")

    p = prepare_progressive_mortality(mort, prog)
    labels = ["M1  Demographics/comorbidity", "M2  + baseline severity", "M3  + day-3 clinical status", "M4  + trajectories/intensity"]
    est = 100 * p["risk_difference"].to_numpy(float)
    lo = 100 * p["rd_lower_95"].to_numpy(float)
    hi = 100 * p["rd_upper_95"].to_numpy(float)
    y = np.arange(len(p))[::-1]
    fig, ax = plt.subplots(figsize=(7.6, 4.2))
    ax.errorbar(est, y, xerr=[est-lo, hi-est], fmt="o", capsize=3, linewidth=1.5)
    ax.axvline(0, linestyle="--", linewidth=1)
    ax.set_yticks(y); ax.set_yticklabels(labels)
    ax.set_xlabel("30-day post-landmark mortality risk difference, percentage points")
    ax.set_title("Progressive adjustment in MIMIC-IV")
    for x, yy in zip(est, y):
        ax.text(x + 0.22, yy, f"{x:+.2f}", va="center", fontsize=8.5)
    savefig(fig, "Fig2_progressive_adjustment_revised")

    outcomes = [
        ("mortality", "30-day post-landmark mortality RD", "percentage points"),
        ("Antibiotic-free days", "Antibiotic-free days", "days"),
        ("Normalized systemic antibiotic exposure", "Normalized systemic exposure", "proportion"),
        ("Normalized broad-spectrum exposure", "Normalized broad-spectrum exposure", "proportion"),
    ]
    fig, axes = plt.subplots(1, 4, figsize=(13.0, 3.5))
    mm = mort.loc[mort["dataset_analysis"].str.startswith("MIMIC-IV primary")].iloc[0]
    pm = mort.loc[mort["dataset_analysis"].str.startswith("PSU modified")].iloc[0]
    for j, (key, title, unit) in enumerate(outcomes):
        ax = axes[j]
        if key == "mortality":
            vals = np.array([100*mm["mortality_rd"], 100*pm["mortality_rd"]], dtype=float)
            los = np.array([100*mm["rd_ci95_low"], 100*pm["rd_ci95_low"]], dtype=float)
            his = np.array([100*mm["rd_ci95_high"], 100*pm["rd_ci95_high"]], dtype=float)
        else:
            s = sec.loc[sec["outcome"] == key].set_index("dataset").loc[["MIMIC-IV", "PSU"]]
            vals = s["estimate"].to_numpy(float)
            los = s["ci95_low"].to_numpy(float)
            his = s["ci95_high"].to_numpy(float)
        yy = np.array([1, 0])
        ax.errorbar(vals, yy, xerr=[vals-los, his-vals], fmt="o", capsize=3, linewidth=1.4)
        ax.axvline(0, linestyle="--", linewidth=1)
        ax.set_yticks(yy)
        ax.set_yticklabels(["MIMIC-IV", "Penn State"] if j == 0 else [])
        ax.set_title(title, fontsize=9.5)
        ax.set_xlabel(unit, fontsize=8.5)
        ax.tick_params(axis="x", labelsize=8)
    fig.suptitle("Cross-dataset mortality and stewardship outcomes", fontsize=11, y=1.02)
    savefig(fig, "Fig3_cross_dataset_outcomes_revised")


def summary(flow: pd.DataFrame, desc: pd.DataFrame, sub: pd.DataFrame, bal: pd.DataFrame) -> None:
    narrowed = desc.loc[desc["deescalation_type"] == "narrowed_or_non_broad_only"].iloc[0]
    stopped = desc.loc[desc["deescalation_type"] == "stopped_all_observed_systemic_antibiotics"].iloc[0]
    text = f"""# Reviewer support package\n\nThis package preserves the frozen primary MIMIC-IV and PSU publication analyses. The only new inferential analyses are explicitly exploratory treatment-subtype mortality sensitivities using the final vital-corrected MIMIC cohort.\n\n## De-escalation composition\n\n- Narrowed/non-broad only: {int(narrowed['n']):,} ({narrowed['percent_within_deescalated']:.1f}% of de-escalated/stopped).\n- Stopped all observed systemic antibiotics: {int(stopped['n']):,} ({stopped['percent_within_deescalated']:.1f}% of de-escalated/stopped).\n\n## Exploratory subtype mortality sensitivities\n\n{sub.to_markdown(index=False)}\n\nThese subtype contrasts are exploratory and do not replace the frozen binary primary estimand.\n\n## Primary diagnostic context\n\n- Maximum post-weighting absolute SMD from the reproduced corrected primary fit: {bal['after'].max():.3f}.\n- Revised main figures use the already frozen manuscript-facing MIMIC/PSU estimates. Figure 2 presents the designated primary M4 confidence interval to avoid two manuscript-facing intervals for the same fully adjusted point estimate.\n\n## Cohort flow\n\n{flow.to_markdown(index=False)}\n"""
    (OUT / "reviewer_support_summary.md").write_text(text, encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for p in [COHORT, FLOW, HARM / "harmonized_mortality_results.csv", HARM / "harmonized_secondary_outcomes.csv", HARM / "mimic_progressive_adjustment.csv"]:
        if not p.exists():
            raise FileNotFoundError(p)
    d = pd.read_csv(COHORT, low_memory=False)
    w, bal = validate_primary(d)
    flow = build_flow()
    desc = subtype_descriptives(d)
    sub = subtype_sensitivity(d, reps=1000)
    diagnostics(w, bal)
    revised_main_figures()
    summary(flow, desc, sub, bal)
    metadata = {
        "corrected_cohort": str(COHORT),
        "n": len(d),
        "deescalated": int((d["A"] == 1).sum()),
        "continued": int((d["A"] == 0).sum()),
        "subtype_bootstrap_reps": 1000,
        "primary_science_changed": False,
    }
    (OUT / "reviewer_support_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
