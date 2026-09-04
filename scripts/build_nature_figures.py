#!/usr/bin/env python3
"""Build publication figures using the Nature-style visual system.

This builder is presentation-only. It consumes frozen manuscript-facing outputs,
uses the shared publication contracts, and does not alter any scientific estimand.
ESM Figure 2 necessarily refits the already-frozen MIMIC propensity-score model to
recover the plotting distributions because no frozen aggregate density artifact yet
exists; no inferential result is recomputed or exported.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde

import figstyle as fs
from figstyle import BLUE, FAINT, INK, MUTED, VERMILLION
from publication_figure_common import (
    apply_publication_secondary_overrides,
    prepare_progressive_mortality,
    pretty_label,
)

HARM = Path("outputs/publication_integration/harmonized")
BALANCE = Path("outputs/publication_integration/reviewer_support/mimic_primary_balance_before_after.csv")
FLOW = Path("outputs/publication_integration/reviewer_support/mimic_cohort_flow_reviewer.csv")
BASE_RUN = Path("outputs/mimic/mimic_iv_v5_7_final_20260820T003506Z")
COHORT = BASE_RUN / "audits/vital_repair/analysis_cohort_vital_corrected.csv"
OUT = Path("outputs/publication_integration/nature_figures")

SHORT_STAGE = {
    "Adult ICU admissions with hospital data": "Adult ICU admissions",
    "Early systemic IV broad-spectrum exposure; alive and hospitalized through 96 h":
        "Early broad-spectrum exposure,\nalive and hospitalized through 96 h",
    "Clinical microbiology sampled and no positive result available by 72 h":
        "Microbiology sampled, no positive\nresult available by 72 h",
    "No active vasopressor overlap during 66-72 h": "No vasopressor overlap, 66-72 h",
    "Systemic IV broad-spectrum coverage during 48-72 h":
        "Broad-spectrum coverage, 48-72 h",
}


def _draw_timeline(ax) -> None:
    """Draw the target-trial timing schematic as panel a of Figure 1."""
    ax.set_xlim(-5, 126)
    ax.set_ylim(-1.10, 1.18)
    ax.axis("off")
    ax.annotate("", xy=(121, 0), xytext=(0, 0),
                arrowprops=dict(arrowstyle="->", linewidth=1.1, color=INK))

    markers = [
        (0, "t0", "First qualifying broad-spectrum\nantibiotic", BLUE),
        (72, "72 h", "Treatment\ndecision", VERMILLION),
        (96, "96 h", "Landmark;\nfollow-up starts", "#2A9D45"),
        (121, "30 d", "Post-landmark\noutcome horizon", "#C62828"),
    ]
    for x, top, label, color in markers:
        ax.plot([x, x], [-0.16, 0.16], color=color, linewidth=1.5)
        ax.text(x, 0.40, top, ha="center", va="bottom", fontsize=8.5,
                fontweight="bold", color=INK)
        ax.text(x, 0.76, label, ha="center", va="bottom", fontsize=6.8,
                linespacing=1.05, color=INK)

    ax.plot([0, 72], [-0.43, -0.43], linewidth=5.0, alpha=0.45,
            color="#9273C5", solid_capstyle="butt")
    ax.text(36, -0.72, "Pre-decision covariates", ha="center", fontsize=7.0)
    ax.plot([72, 96], [-0.43, -0.43], linewidth=5.0, alpha=0.68,
            color="#A66C5B", solid_capstyle="butt")
    ax.text(84, -0.72, "Treatment classification", ha="center", fontsize=7.0)
    ax.plot([96, 121], [-0.43, -0.43], linewidth=5.0, alpha=0.68,
            color="#D454B6", solid_capstyle="butt")
    ax.text(108.5, -0.72, "Outcome follow-up", ha="center", fontsize=7.0)
    fs.panel_label(ax, "a", dx=-0.03, dy=0.98)


def _draw_attrition(ax, f: pd.DataFrame) -> None:
    """Draw sequential cohort attrition as panel b of Figure 1."""
    stages = f.loc[f["stage"].isin(SHORT_STAGE)].copy()
    n = stages["n"].to_numpy(float)
    prev = np.concatenate([[n[0]], n[:-1]])
    y = np.arange(len(stages))[::-1]

    col_n = n[0] * 1.04
    col_rm = n[0] * 1.34
    ax.barh(y, prev, color=FAINT, alpha=0.35, height=0.62, linewidth=0)
    ax.barh(y, n, color=MUTED, height=0.62, linewidth=0)
    for ni, pi, yi in zip(n, prev, y):
        ax.text(col_n, yi, f"{int(ni):,}", va="center", ha="left",
                fontsize=6.8, color=INK)
        if pi > ni:
            ax.text(col_rm, yi, f"−{int(pi - ni):,}", va="center", ha="right",
                    fontsize=6.6, color=MUTED)
    ax.text(col_n, len(stages) - 0.45, "Retained", fontsize=6.6,
            fontweight="bold", ha="left", va="center", color=INK)
    ax.text(col_rm, len(stages) - 0.45, "Excluded", fontsize=6.6,
            fontweight="bold", ha="right", va="center", color=MUTED)
    ax.set_yticks(y)
    ax.set_yticklabels([SHORT_STAGE[s] for s in stages["stage"]],
                       fontsize=6.6, linespacing=1.12)
    ax.set_xlim(0, col_rm * 1.02)
    ax.set_ylim(-0.70, len(stages) - 0.15)
    ax.set_xticks([])
    ax.spines["bottom"].set_visible(False)
    fs.strip_y_axis(ax)
    fs.panel_label(ax, "b", dx=-0.39, dy=0.99)


def _draw_treatment_split(ax, f: pd.DataFrame) -> None:
    """Draw the final analytic treatment split as panel c of Figure 1."""
    arms = f.loc[f["stage"].isin(["De-escalated/stopped", "Continued broad-spectrum"])]
    arm_n = arms.set_index("stage").loc[
        ["De-escalated/stopped", "Continued broad-spectrum"], "n"
    ].to_numpy(float)
    total = arm_n.sum()

    ax.barh([0], [arm_n[0]], color=BLUE, height=0.40, linewidth=0)
    ax.barh([0], [arm_n[1]], left=[arm_n[0]], color=VERMILLION,
            height=0.40, linewidth=0)
    for centre, val, color, label in [
        (arm_n[0] / 2, arm_n[0], BLUE, "De-escalated or stopped"),
        (arm_n[0] + arm_n[1] / 2, arm_n[1], VERMILLION,
         "Continued broad-spectrum"),
    ]:
        ax.text(centre, -0.40,
                f"{label}\n{int(val):,} ({100 * val / total:.1f}%)",
                ha="center", va="top", fontsize=6.8, color=color, linespacing=1.25)
    ax.set_xlim(0, total * 1.08)
    ax.set_ylim(-1.40, 0.42)
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ("bottom", "left"):
        ax.spines[s].set_visible(False)
    ax.text(-0.015, 0, f"Analytic cohort\nn = {int(total):,}",
            transform=ax.get_yaxis_transform(), fontsize=6.8,
            ha="right", va="center", color=INK, linespacing=1.15)
    fs.panel_label(ax, "c", dx=-0.39, dy=0.78)


def build_fig1() -> None:
    """Build one complete Figure 1: timeline, attrition, and treatment split."""
    f = pd.read_csv(FLOW)
    fig = plt.figure(figsize=(fs.DOUBLE, 5.15))
    gs = fig.add_gridspec(3, 1, height_ratios=[1.15, 2.55, 0.70], hspace=0.30)
    _draw_timeline(fig.add_subplot(gs[0]),)
    _draw_attrition(fig.add_subplot(gs[1]), f)
    _draw_treatment_split(fig.add_subplot(gs[2]), f)
    fig.subplots_adjust(left=0.30, right=0.97, top=0.97, bottom=0.06)
    fs.savefig(fig, OUT, "Fig1_target_trial_and_cohort")


def build_fig2() -> None:
    mort = pd.read_csv(HARM / "harmonized_mortality_results.csv")
    prog = pd.read_csv(HARM / "mimic_progressive_adjustment.csv")
    p = prepare_progressive_mortality(mort, prog)

    labels = [
        "M1  Demographics, comorbidity",
        "M2  + baseline severity",
        "M3  + day-3 clinical status",
        "M4  + trajectories, intensity",
    ]
    est = 100 * p["risk_difference"].to_numpy(float)
    lo = 100 * p["rd_lower_95"].to_numpy(float)
    hi = 100 * p["rd_upper_95"].to_numpy(float)
    y = np.arange(len(p))[::-1]

    fig, ax = plt.subplots(figsize=(fs.DOUBLE, 2.45))
    fig.subplots_adjust(left=0.22, right=0.72, top=0.90, bottom=0.34)
    fs.null_line(ax)
    ax.plot(est, y, color=FAINT, linewidth=0.8, zorder=1,
            solid_capstyle="round")

    colors = [MUTED] * (len(y) - 1) + [BLUE]
    sizes = [4.0] * (len(y) - 1) + [5.4]
    for xi, li, hi_i, yi, c, s in zip(est, lo, hi, y, colors, sizes):
        ax.plot([li, hi_i], [yi, yi], color=c, linewidth=1.1,
                solid_capstyle="butt", zorder=2)
        for cap in (li, hi_i):
            ax.plot([cap, cap], [yi - 0.14, yi + 0.14], color=c,
                    linewidth=1.1, zorder=2)
        ax.plot([xi], [yi], "o", color=c, markersize=s, zorder=3,
                markeredgecolor="white", markeredgewidth=0.45)

    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=7.6)
    ax.set_ylim(-0.70, len(y) - 0.30)
    ax.set_xlim(-7.5, 7.5)
    ax.set_xticks([-6, -4, -2, 0, 2, 4, 6])
    ax.tick_params(axis="x", labelsize=7.0)
    ax.set_xlabel("30-day mortality risk difference (percentage points)",
                  fontsize=7.8, labelpad=22)
    fs.strip_y_axis(ax)

    ax.text(1.06, 1.0, "RD (95% CI)", transform=ax.transAxes, fontsize=7.6,
            fontweight="bold", va="bottom", ha="left", color=INK)
    for xi, li, hi_i, yi, c in zip(est, lo, hi, y, colors):
        ax.text(1.06, yi, f"{xi:+.2f} ({li:+.2f}, {hi_i:+.2f})",
                transform=ax.get_yaxis_transform(), fontsize=7.4,
                va="center", ha="left", color=INK if c == BLUE else MUTED,
                clip_on=False)

    ax.text(0.0, -0.29, "← favours de-escalation", transform=ax.transAxes,
            fontsize=6.8, color=MUTED, ha="left", va="center")
    ax.text(1.0, -0.29, "favours continuation →", transform=ax.transAxes,
            fontsize=6.8, color=MUTED, ha="right", va="center")
    fs.savefig(fig, OUT, "Fig2_progressive_adjustment")


def build_fig3() -> None:
    """Build the accepted 2x2 cross-dataset outcome figure."""
    mort = pd.read_csv(HARM / "harmonized_mortality_results.csv")
    sec = apply_publication_secondary_overrides(
        pd.read_csv(HARM / "harmonized_secondary_outcomes.csv")
    )
    mm = mort.loc[mort["dataset_analysis"].str.startswith("MIMIC-IV primary")].iloc[0]
    pm = mort.loc[mort["dataset_analysis"].str.startswith("PSU modified")].iloc[0]

    panels = [
        ("mortality", "30-day mortality risk difference (percentage points)", 2),
        ("Antibiotic-free days", "Antibiotic-free days (difference in days)", 2),
        ("Normalized systemic antibiotic exposure",
         "Systemic antibiotic exposure (difference in proportion)", 3),
        ("Normalized broad-spectrum exposure",
         "Broad-spectrum exposure (difference in proportion)", 3),
    ]
    n_lab = [f"MIMIC-IV\nn = {int(mm['cohort_n']):,}",
             f"Penn State\nn = {int(pm['cohort_n']):,}"]

    fig, axes = plt.subplots(2, 2, figsize=(fs.DOUBLE, 3.15))
    axes = axes.ravel()
    for j, (key, unit, decimals) in enumerate(panels):
        ax = axes[j]
        if key == "mortality":
            vals = np.array([100 * mm["mortality_rd"], 100 * pm["mortality_rd"]], float)
            los = np.array([100 * mm["rd_ci95_low"], 100 * pm["rd_ci95_low"]], float)
            his = np.array([100 * mm["rd_ci95_high"], 100 * pm["rd_ci95_high"]], float)
        else:
            s = sec.loc[sec["outcome"] == key].set_index("dataset").loc[["MIMIC-IV", "PSU"]]
            vals = s["estimate"].to_numpy(float)
            los = s["ci95_low"].to_numpy(float)
            his = s["ci95_high"].to_numpy(float)

        fs.null_line(ax)
        yy = np.array([1, 0])
        colors = [BLUE, VERMILLION]
        for xi, li, hi_i, yi, c in zip(vals, los, his, yy, colors):
            ax.plot([li, hi_i], [yi, yi], color=c, linewidth=1.0,
                    solid_capstyle="butt", zorder=2)
            for cap in (li, hi_i):
                ax.plot([cap, cap], [yi - 0.1, yi + 0.1], color=c,
                        linewidth=1.0, zorder=2)
            ax.plot([xi], [yi], "o", color=c, markersize=4.2, zorder=3,
                    markeredgecolor="white", markeredgewidth=0.4)

        span = max(his.max(), 0) - min(los.min(), 0)
        pad = 0.18 * span
        ax.set_xlim(min(los.min(), 0) - pad, max(his.max(), 0) + pad)
        fmt = f"{{:+.{decimals}f}}"
        for xi, hi_i, yi, c in zip(vals, his, yy, colors):
            ax.text(hi_i + 0.04 * span, yi, fmt.format(xi), fontsize=6,
                    va="center", ha="left", color=c)

        ax.set_yticks(yy)
        ax.set_yticklabels(n_lab, fontsize=6)
        ax.set_ylim(-0.55, 1.55)
        ax.set_xlabel(unit, fontsize=6.5)
        fs.strip_y_axis(ax)
        fs.panel_label(ax, "abcd"[j], dx=-0.26)

    fig.subplots_adjust(left=0.11, right=0.98, top=0.92, bottom=0.13,
                        hspace=1.05, wspace=0.42)
    fs.savefig(fig, OUT, "Fig3_cross_dataset_outcomes")


def build_esm1() -> None:
    """Build absolute-SMD dumbbell plot with readable, unclipped annotations."""
    bal = pd.read_csv(BALANCE)
    top = bal.assign(
        before_abs=pd.to_numeric(bal["before"], errors="coerce").abs(),
        after_abs=pd.to_numeric(bal["after"], errors="coerce").abs(),
    ).sort_values("before_abs", ascending=False).head(35).sort_values("before_abs")

    before = top["before_abs"].to_numpy(float)
    after = top["after_abs"].to_numpy(float)
    y = np.arange(len(top))
    max_x = max(float(before.max()), float(after.max()), 0.1)
    right = max_x * 1.08
    left = -max_x * 0.025

    fig, ax = plt.subplots(figsize=(fs.ONE_HALF, 6.15))
    fig.subplots_adjust(left=0.46, right=0.96, top=0.955, bottom=0.08)
    ax.axvline(0.1, color=fs.RULE, linestyle=(0, (3, 2)), linewidth=0.6, zorder=0)

    for b, a, yi in zip(before, after, y):
        ax.plot([b, a], [yi, yi], color=FAINT, linewidth=0.8,
                zorder=1, solid_capstyle="round")
    ax.plot(before, y, "o", color=FAINT, markersize=3.2, zorder=2,
            markeredgecolor="white", markeredgewidth=0.3)
    ax.plot(after, y, "o", color=BLUE, markersize=3.6, zorder=3,
            markeredgecolor="white", markeredgewidth=0.3)

    ax.set_yticks(y)
    ax.set_yticklabels([pretty_label(v) for v in top["variable"]], fontsize=5.8)
    ax.set_ylim(-0.8, len(y) - 0.05)
    ax.set_xlim(left, right)
    ax.set_xlabel("Absolute standardized mean difference", fontsize=7.0)
    fs.strip_y_axis(ax)

    # Keep labels entirely inside the figure. Their horizontal locations also
    # reinforce which endpoint belongs to which series without needing a legend.
    after_anchor = min(max(float(after[-1]) / right, 0.03), 0.22)
    ax.text(after_anchor, 1.002, "after weighting", transform=ax.transAxes,
            fontsize=6.5, color=BLUE, ha="center", va="bottom", clip_on=False)
    ax.text(0.98, 1.002, "before weighting", transform=ax.transAxes,
            fontsize=6.5, color=MUTED, ha="right", va="bottom", clip_on=False)

    ax.text(0.1 + 0.012 * max_x, 0.02, "0.10 balance threshold",
            transform=ax.get_xaxis_transform(), rotation=90, fontsize=6.3,
            color=MUTED, va="bottom", ha="left")
    fs.savefig(fig, OUT, "ESM_Fig1_covariate_balance")


def build_esm2() -> None:
    """Build mirrored propensity-score and stabilized-weight densities."""
    from sepsis_deescalation.specification import CANDIDATE_PS_VARS
    from sepsis_deescalation.stats import fit_stabilized_iptw

    d = pd.read_csv(COHORT, low_memory=False)
    w, _, _ = fit_stabilized_iptw(d, CANDIDATE_PS_VARS)
    fig, axes = plt.subplots(1, 2, figsize=(fs.DOUBLE, 2.35))
    groups = [
        (1, "De-escalated/stopped", BLUE),
        (0, "Continued broad-spectrum", VERMILLION),
    ]

    for panel_index, (ax, (var, xlabel, lo, hi)) in enumerate(zip(
        axes,
        [("ps_den", "Estimated propensity for de-escalation", 0.0, 1.0),
         ("SW_A", "Stabilized IPTW", 0.0, None)],
    )):
        series = {
            a: pd.to_numeric(w.loc[w["A"] == a, var], errors="coerce")
            .dropna().to_numpy(float)
            for a, _, _ in groups
        }
        if hi is None:
            allw = np.concatenate(list(series.values()))
            hi = float(np.percentile(allw, 99.5))
            n_beyond = int((allw > hi).sum())
            ax.text(
                0.97, 0.52,
                f"{n_beyond} weights > {hi:.1f} not shown\n(max {allw.max():.1f})",
                transform=ax.transAxes, fontsize=7.0, color=MUTED,
                ha="right", va="center", linespacing=1.2,
            )
        grid = np.linspace(lo, hi, 512)
        for sign, (a, _, color) in zip((1, -1), groups):
            dens = gaussian_kde(series[a])(grid) * sign
            ax.fill_between(grid, 0, dens, color=color, alpha=0.30, linewidth=0)
            ax.plot(grid, dens, color=color, linewidth=1.0)

        ax.axhline(0, color=INK, linewidth=0.5)
        ax.set_xlim(lo, hi)
        ax.set_xlabel(xlabel, fontsize=7.6)
        ax.set_ylabel("Density", fontsize=7.4)
        ax.set_yticks([])
        ax.spines["left"].set_visible(False)
        fs.panel_label(ax, "ab"[panel_index], dx=-0.06)

        yl = ax.get_ylim()
        ax.text(0.97, 0.91, "de-escalated/stopped", transform=ax.transAxes,
                fontsize=7.2, color=BLUE, ha="right", va="top")
        ax.text(0.97, 0.09, "continued broad-spectrum", transform=ax.transAxes,
                fontsize=7.2, color=VERMILLION, ha="right", va="bottom")
        ax.set_ylim(yl)

    fig.subplots_adjust(left=0.06, right=0.98, top=0.88, bottom=0.22, wspace=0.16)
    fs.savefig(fig, OUT, "ESM_Fig2_propensity_overlap")


def main() -> None:
    fs.use_nature_style()
    OUT.mkdir(parents=True, exist_ok=True)

    required = [
        HARM / "harmonized_mortality_results.csv",
        HARM / "harmonized_secondary_outcomes.csv",
        HARM / "mimic_progressive_adjustment.csv",
        FLOW,
        BALANCE,
        COHORT,
    ]
    missing = [str(p) for p in required if not p.exists()]
    if missing:
        raise FileNotFoundError("Missing required inputs: " + ", ".join(missing))

    build_fig1()
    build_fig2()
    build_fig3()
    build_esm1()
    build_esm2()
    print("Built Figure 1-3 and ESM Figures 1-2")


if __name__ == "__main__":
    main()
