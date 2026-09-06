#!/usr/bin/env python3
"""Build manuscript-facing publication figures from frozen project outputs.

This builder is presentation-only. It consumes frozen manuscript-facing outputs,
uses shared publication contracts, and does not alter any scientific estimand.
ESM Figure 2 refits the already-frozen MIMIC propensity-score model only to recover
the plotting distributions; no inferential result is recomputed or exported.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde

import figstyle as fs
from figstyle import BLUE, FAINT, GREEN, INK, MUTED, RULE, VERMILLION
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

FLOW_STAGE_ORDER = [
    "Adult ICU admissions with hospital data",
    "Early systemic IV broad-spectrum exposure; alive and hospitalized through 96 h",
    "Clinical microbiology sampled and no positive result available by 72 h",
    "No active vasopressor overlap during 66-72 h",
    "Systemic IV broad-spectrum coverage during 48-72 h",
]

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
    """Draw a sparse target-trial timing schematic as panel a of Figure 1."""
    ax.set_xlim(-5, 126)
    ax.set_ylim(-1.00, 1.15)
    ax.axis("off")

    # Main decision timeline. The post-landmark period is schematic rather than
    # proportional to calendar time.
    ax.annotate(
        "",
        xy=(121, 0.05),
        xytext=(0, 0.05),
        arrowprops=dict(arrowstyle="->", linewidth=0.9, color=INK),
    )

    markers = [
        (0, "t0", "First qualifying broad-spectrum\nantibiotic exposure", MUTED, 1.0),
        (55, "72 h", "Treatment\ndecision", BLUE, 1.8),
        (86, "96 h", "Landmark\nfollow-up starts", INK, 1.5),
        (121, "30 d", "Mortality\nhorizon", MUTED, 1.0),
    ]
    for x, top, label, color, width in markers:
        ax.plot([x, x], [-0.10, 0.20], color=color, linewidth=width, zorder=3)
        ax.text(
            x,
            0.42,
            top,
            ha="center",
            va="bottom",
            fontsize=8.0,
            fontweight="bold",
            color=INK,
        )
        ax.text(
            x,
            0.72,
            label,
            ha="center",
            va="bottom",
            fontsize=6.5,
            linespacing=1.05,
            color=INK,
        )

    bands = [
        (0, 55, "Pre-decision covariates", "#E4E4E4"),
        (55, 86, "Treatment window", "#D9EBF5"),
        (86, 121, "Outcome follow-up", "#ECECEC"),
    ]
    for left, right, label, color in bands:
        ax.add_patch(
            Rectangle(
                (left, -0.50),
                right - left,
                0.16,
                facecolor=color,
                edgecolor="none",
                zorder=1,
            )
        )
        ax.text(
            (left + right) / 2,
            -0.69,
            label,
            ha="center",
            va="top",
            fontsize=6.3,
            color=INK,
        )
    fs.panel_label(ax, "a", dx=-0.02, dy=0.98)


def _draw_attrition(ax, f: pd.DataFrame) -> None:
    """Draw sequential MIMIC-IV cohort attrition as panel b of Figure 1."""
    by_stage = f.set_index("stage")
    missing = [s for s in FLOW_STAGE_ORDER if s not in by_stage.index]
    if missing:
        raise RuntimeError("Missing cohort-flow stages: " + ", ".join(missing))

    stages = by_stage.loc[FLOW_STAGE_ORDER].reset_index()
    n = stages["n"].to_numpy(float)
    if np.any(np.diff(n) > 0):
        raise RuntimeError("Cohort-flow retained counts must be non-increasing")

    prev = np.concatenate([[n[0]], n[:-1]])
    y = np.arange(len(stages))[::-1]
    max_n = float(n[0])
    retained_col = max_n * 1.035
    excluded_col = max_n * 1.31

    for idx, (ni, pi, yi) in enumerate(zip(n, prev, y)):
        ax.barh(
            yi,
            pi,
            color="#E9E9E9",
            height=0.58,
            linewidth=0,
            zorder=1,
        )
        ax.barh(
            yi,
            ni,
            color=INK if idx == len(n) - 1 else MUTED,
            height=0.58,
            linewidth=0,
            zorder=2,
        )
        ax.text(
            retained_col,
            yi,
            f"{int(ni):,}",
            va="center",
            ha="left",
            fontsize=6.7,
            color=INK,
            fontweight="bold" if idx == len(n) - 1 else "normal",
        )
        if pi > ni:
            ax.text(
                excluded_col,
                yi,
                f"-{int(pi - ni):,}",
                va="center",
                ha="right",
                fontsize=6.4,
                color=MUTED,
            )

    ax.text(
        retained_col,
        len(stages) - 0.42,
        "Retained",
        fontsize=6.4,
        fontweight="bold",
        ha="left",
        va="center",
        color=INK,
    )
    ax.text(
        excluded_col,
        len(stages) - 0.42,
        "Excluded",
        fontsize=6.4,
        fontweight="bold",
        ha="right",
        va="center",
        color=MUTED,
    )
    ax.set_yticks(y)
    tick_labels = [SHORT_STAGE[s] for s in stages["stage"]]
    ax.set_yticklabels(tick_labels, fontsize=6.4, linespacing=1.10)
    ax.get_yticklabels()[-1].set_fontweight("bold")
    ax.set_xlim(0, excluded_col * 1.025)
    ax.set_ylim(-0.65, len(stages) - 0.15)
    ax.set_xticks([])
    ax.spines["bottom"].set_visible(False)
    fs.strip_y_axis(ax)
    fs.panel_label(ax, "b", dx=-0.38, dy=0.99)


def _draw_treatment_split(ax, f: pd.DataFrame) -> None:
    """Draw the final treatment split with treatment-specific, local color semantics."""
    arms = f.loc[f["stage"].isin(["De-escalated/stopped", "Continued broad-spectrum"])]
    arm_n = arms.set_index("stage").loc[
        ["De-escalated/stopped", "Continued broad-spectrum"], "n"
    ].to_numpy(float)
    total = float(arm_n.sum())
    if total <= 0:
        raise RuntimeError("Treatment split must contain a positive analytic cohort")
    pct = 100.0 * arm_n / total

    ax.barh([0], [pct[0]], color=GREEN, height=0.34, linewidth=0)
    ax.barh([0], [pct[1]], left=[pct[0]], color=MUTED, height=0.34, linewidth=0)

    ax.text(
        pct[0] / 2,
        0,
        f"{pct[0]:.1f}%",
        ha="center",
        va="center",
        fontsize=6.4,
        fontweight="bold",
        color="white",
    )
    ax.text(
        pct[0] + pct[1] / 2,
        0,
        f"{pct[1]:.1f}%",
        ha="center",
        va="center",
        fontsize=6.4,
        fontweight="bold",
        color="white",
    )
    ax.text(
        0,
        -0.42,
        f"De-escalated or stopped\n{int(arm_n[0]):,}",
        ha="left",
        va="top",
        fontsize=6.3,
        color=INK,
        linespacing=1.15,
    )
    ax.text(
        100,
        -0.42,
        f"Continued broad-spectrum\n{int(arm_n[1]):,}",
        ha="right",
        va="top",
        fontsize=6.3,
        color=INK,
        linespacing=1.15,
    )
    ax.text(
        0,
        0.50,
        f"Analytic cohort, n = {int(total):,}",
        ha="left",
        va="bottom",
        fontsize=6.7,
        fontweight="bold",
        color=INK,
    )
    ax.set_xlim(0, 100)
    ax.set_ylim(-1.05, 0.78)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ("bottom", "left"):
        ax.spines[spine].set_visible(False)
    fs.panel_label(ax, "c", dx=-0.06, dy=0.98)


def build_fig1() -> None:
    """Build Figure 1: target-trial timing, cohort attrition, and treatment split."""
    f = pd.read_csv(FLOW)
    fig = plt.figure(figsize=(fs.DOUBLE, 4.35))
    gs = fig.add_gridspec(
        3,
        1,
        height_ratios=[1.00, 2.15, 0.78],
        hspace=0.34,
    )
    _draw_timeline(fig.add_subplot(gs[0]))
    _draw_attrition(fig.add_subplot(gs[1]), f)
    _draw_treatment_split(fig.add_subplot(gs[2]), f)
    fig.subplots_adjust(left=0.30, right=0.97, top=0.97, bottom=0.07)
    fs.savefig(fig, OUT, "Fig1_target_trial_and_cohort")


def build_fig2() -> None:
    """Build Figure 2: progressive adjustment of the MIMIC-IV mortality association."""
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

    fig, ax = plt.subplots(figsize=(fs.DOUBLE, 2.60))
    fig.subplots_adjust(left=0.235, right=0.72, top=0.90, bottom=0.34)
    fs.null_line(ax)
    ax.plot(est, y, color="#D0D0D0", linewidth=0.55, zorder=1, solid_capstyle="round")

    colors = [MUTED] * (len(y) - 1) + [BLUE]
    sizes = [3.8] * (len(y) - 1) + [5.2]
    for xi, li, hi_i, yi, color, size in zip(est, lo, hi, y, colors, sizes):
        ax.plot([li, hi_i], [yi, yi], color=color, linewidth=0.95, solid_capstyle="butt", zorder=2)
        for cap in (li, hi_i):
            ax.plot([cap, cap], [yi - 0.13, yi + 0.13], color=color, linewidth=0.95, zorder=2)
        ax.plot(
            [xi],
            [yi],
            "o",
            color=color,
            markersize=size,
            zorder=3,
            markeredgecolor="white",
            markeredgewidth=0.4,
        )

    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=7.4)
    ax.get_yticklabels()[-1].set_fontweight("bold")
    ax.set_ylim(-0.72, len(y) - 0.28)
    ax.set_xlim(-7.5, 7.5)
    ax.set_xticks([-6, -4, -2, 0, 2, 4, 6])
    ax.tick_params(axis="x", labelsize=6.8)
    ax.set_xlabel("30-day mortality risk difference (percentage points)", fontsize=7.5, labelpad=21)
    fs.strip_y_axis(ax)

    ax.text(
        1.06,
        1.0,
        "RD (95% CI), pp",
        transform=ax.transAxes,
        fontsize=7.2,
        fontweight="bold",
        va="bottom",
        ha="left",
        color=INK,
    )
    for xi, li, hi_i, yi, color in zip(est, lo, hi, y, colors):
        ax.text(
            1.06,
            yi,
            f"{xi:+.2f} ({li:+.2f}, {hi_i:+.2f})",
            transform=ax.get_yaxis_transform(),
            fontsize=7.1,
            va="center",
            ha="left",
            color=INK if color == BLUE else MUTED,
            clip_on=False,
        )
    ax.text(
        1.06,
        y[-1] - 0.36,
        "designated primary model",
        transform=ax.get_yaxis_transform(),
        fontsize=5.9,
        va="center",
        ha="left",
        color=MUTED,
        clip_on=False,
    )

    ax.text(
        0.0,
        -0.29,
        "favours de-escalation/stopping",
        transform=ax.transAxes,
        fontsize=6.4,
        color=MUTED,
        ha="left",
        va="center",
    )
    ax.text(
        1.0,
        -0.29,
        "favours continuation",
        transform=ax.transAxes,
        fontsize=6.4,
        color=MUTED,
        ha="right",
        va="center",
    )
    fs.savefig(fig, OUT, "Fig2_progressive_adjustment")


def _effect_limits(los: np.ndarray, his: np.ndarray) -> tuple[float, float]:
    """Return compact limits that include the null and leave room for CI caps."""
    lo = min(float(np.nanmin(los)), 0.0)
    hi = max(float(np.nanmax(his)), 0.0)
    span = hi - lo
    if span <= 0:
        span = max(abs(lo), abs(hi), 1.0)
    pad = 0.12 * span
    return lo - pad, hi + pad


def _draw_two_dataset_effect_row(
    ax,
    *,
    title: str,
    xlabel: str,
    vals: np.ndarray,
    los: np.ndarray,
    his: np.ndarray,
    decimals: int,
) -> None:
    """Draw one outcome row using the same MIMIC-IV/Penn-State grammar."""
    fs.null_line(ax)
    y = np.array([1.0, 0.0])
    colors = [BLUE, VERMILLION]

    for xi, li, hi_i, yi, color in zip(vals, los, his, y, colors):
        ax.plot([li, hi_i], [yi, yi], color=color, linewidth=0.95, solid_capstyle="butt", zorder=2)
        for cap in (li, hi_i):
            ax.plot([cap, cap], [yi - 0.11, yi + 0.11], color=color, linewidth=0.95, zorder=2)
        ax.plot(
            [xi],
            [yi],
            "o",
            color=color,
            markersize=4.4,
            markeredgecolor="white",
            markeredgewidth=0.4,
            zorder=3,
        )

    ax.set_xlim(*_effect_limits(los, his))
    ax.set_ylim(-0.45, 1.45)
    ax.set_yticks(y)
    ax.set_yticklabels(["MIMIC-IV", "Penn State"], fontsize=6.4)
    fs.strip_y_axis(ax)
    ax.tick_params(axis="x", labelsize=6.3)
    ax.set_xlabel(xlabel, fontsize=6.6, labelpad=3)
    ax.text(
        0,
        1.10,
        title,
        transform=ax.transAxes,
        fontsize=7.0,
        fontweight="bold",
        ha="left",
        va="bottom",
        color=INK,
    )

    fmt = f"{{:+.{decimals}f}} ({{:+.{decimals}f}}, {{:+.{decimals}f}})"
    for xi, li, hi_i, yi in zip(vals, los, his, y):
        ax.text(
            1.04,
            yi,
            fmt.format(xi, li, hi_i),
            transform=ax.get_yaxis_transform(),
            fontsize=6.5,
            color=INK,
            ha="left",
            va="center",
            clip_on=False,
        )


def build_fig3() -> None:
    """Build Figure 3 as an integrated cross-dataset comparison."""
    mort = pd.read_csv(HARM / "harmonized_mortality_results.csv")
    sec = apply_publication_secondary_overrides(
        pd.read_csv(HARM / "harmonized_secondary_outcomes.csv")
    )
    mm = mort.loc[mort["dataset_analysis"].str.startswith("MIMIC-IV primary")].iloc[0]
    pm = mort.loc[mort["dataset_analysis"].str.startswith("PSU modified")].iloc[0]

    rows = []
    rows.append(
        {
            "title": "30-day post-landmark mortality",
            "xlabel": "Risk difference (percentage points)",
            "decimals": 2,
            "vals": np.array([100 * mm["mortality_rd"], 100 * pm["mortality_rd"]], float),
            "los": np.array([100 * mm["rd_ci95_low"], 100 * pm["rd_ci95_low"]], float),
            "his": np.array([100 * mm["rd_ci95_high"], 100 * pm["rd_ci95_high"]], float),
        }
    )

    secondary_rows = [
        (
            "Antibiotic-free days",
            "Antibiotic-free days",
            "Difference in days",
            2,
        ),
        (
            "Normalized systemic antibiotic exposure",
            "Normalized systemic antibiotic exposure",
            "Difference in proportion",
            3,
        ),
        (
            "Normalized broad-spectrum exposure",
            "Normalized broad-spectrum exposure",
            "Difference in proportion",
            3,
        ),
    ]
    for outcome, title, xlabel, decimals in secondary_rows:
        s = sec.loc[sec["outcome"] == outcome].set_index("dataset").loc[["MIMIC-IV", "PSU"]]
        rows.append(
            {
                "title": title,
                "xlabel": xlabel,
                "decimals": decimals,
                "vals": s["estimate"].to_numpy(float),
                "los": s["ci95_low"].to_numpy(float),
                "his": s["ci95_high"].to_numpy(float),
            }
        )

    fig = plt.figure(figsize=(fs.DOUBLE, 5.75))
    gs = fig.add_gridspec(
        6,
        1,
        height_ratios=[0.18, 1.0, 0.22, 1.0, 1.0, 1.0],
        hspace=1.40,
    )

    head1 = fig.add_subplot(gs[0])
    head1.axis("off")
    head1.text(0, 0.55, "Clinical outcome", fontsize=7.2, fontweight="bold", ha="left", va="center")
    ax_mort = fig.add_subplot(gs[1])

    head2 = fig.add_subplot(gs[2])
    head2.axis("off")
    head2.text(0, 0.55, "Stewardship outcomes", fontsize=7.2, fontweight="bold", ha="left", va="center")
    axes = [ax_mort, fig.add_subplot(gs[3]), fig.add_subplot(gs[4]), fig.add_subplot(gs[5])]

    for panel_index, (ax, row) in enumerate(zip(axes, rows)):
        _draw_two_dataset_effect_row(ax, **row)
        fs.panel_label(ax, "abcd"[panel_index], dx=-0.17, dy=1.08)

    handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=BLUE, markeredgecolor="white",
               markeredgewidth=0.4, markersize=5, label=f"MIMIC-IV (n = {int(mm['cohort_n']):,})"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=VERMILLION, markeredgecolor="white",
               markeredgewidth=0.4, markersize=5, label=f"Penn State (n = {int(pm['cohort_n']):,})"),
    ]
    fig.legend(
        handles=handles,
        loc="upper right",
        bbox_to_anchor=(0.97, 0.987),
        ncol=2,
        frameon=False,
        handletextpad=0.35,
        columnspacing=1.2,
        fontsize=6.3,
    )
    fig.text(
        0.77,
        0.948,
        "Modified external replication; estimates are not pooled",
        ha="right",
        va="top",
        fontsize=5.8,
        color=MUTED,
    )
    fig.subplots_adjust(left=0.22, right=0.77, top=0.94, bottom=0.08)
    fs.savefig(fig, OUT, "Fig3_cross_dataset_outcomes")


def build_esm1() -> None:
    """Build the MIMIC-IV absolute-SMD balance diagnostic."""
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
    fig.subplots_adjust(left=0.46, right=0.96, top=0.94, bottom=0.08)
    ax.axvline(0.1, color=RULE, linestyle=(0, (3, 2)), linewidth=0.6, zorder=0)

    for b, a, yi in zip(before, after, y):
        ax.plot([b, a], [yi, yi], color=FAINT, linewidth=0.8, zorder=1, solid_capstyle="round")
    ax.plot(
        before, y, "o", color=FAINT, markersize=3.2, zorder=2,
        markeredgecolor="white", markeredgewidth=0.3,
    )
    ax.plot(after, y, "o", color=BLUE, markersize=3.6, zorder=3, markeredgecolor="white", markeredgewidth=0.3)

    ax.set_yticks(y)
    ax.set_yticklabels([pretty_label(v) for v in top["variable"]], fontsize=5.8)
    ax.set_ylim(-0.8, len(y) + 0.75)
    ax.set_xlim(left, right)
    ax.set_xlabel("Absolute standardized mean difference", fontsize=7.0)
    fs.strip_y_axis(ax)

    after_anchor = min(max(float(after[-1]) / right, 0.03), 0.22)
    ax.text(after_anchor, 1.002, "after weighting", transform=ax.transAxes,
            fontsize=6.5, color=INK, ha="center", va="bottom", clip_on=False)
    ax.text(0.98, 1.002, "before weighting", transform=ax.transAxes,
            fontsize=6.5, color=INK, ha="right", va="bottom", clip_on=False)
    ax.text(0.1 + 0.012 * max_x, 0.02, "0.10 threshold",
            transform=ax.get_xaxis_transform(), rotation=90, fontsize=6.0,
            color=MUTED, va="bottom", ha="left")
    ax.text(
        0.02,
        1.045,
        f"Maximum post-weight |SMD| = {float(np.nanmax(after)):.3f}",
        transform=ax.transAxes,
        fontsize=6.1,
        color=MUTED,
        ha="left",
        va="bottom",
        clip_on=False,
    )
    fs.savefig(fig, OUT, "ESM_Fig1_covariate_balance")


def _effective_sample_size(values: np.ndarray) -> float:
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return float("nan")
    denom = float(np.sum(x * x))
    if denom <= 0:
        return float("nan")
    return float(np.sum(x) ** 2 / denom)


def build_esm2() -> None:
    """Build mirrored propensity-score and stabilized-weight diagnostics."""
    from sepsis_deescalation.specification import CANDIDATE_PS_VARS
    from sepsis_deescalation.stats import fit_stabilized_iptw

    d = pd.read_csv(COHORT, low_memory=False)
    w, _, _ = fit_stabilized_iptw(d, CANDIDATE_PS_VARS)
    fig, axes = plt.subplots(1, 2, figsize=(fs.DOUBLE, 2.40))
    groups = [
        (1, "De-escalated/stopped", GREEN),
        (0, "Continued broad-spectrum", MUTED),
    ]

    for panel_index, (ax, (var, xlabel, lo, hi)) in enumerate(
        zip(
            axes,
            [
                ("ps_den", "Estimated propensity for de-escalation", 0.0, 1.0),
                ("SW_A", "Stabilized IPTW", 0.0, None),
            ],
        )
    ):
        series = {
            a: pd.to_numeric(w.loc[w["A"] == a, var], errors="coerce").dropna().to_numpy(float)
            for a, _, _ in groups
        }
        n_beyond = None
        max_weight = None
        if hi is None:
            allw = np.concatenate(list(series.values()))
            hi = float(np.percentile(allw, 99.5))
            n_beyond = int((allw > hi).sum())
            max_weight = float(np.nanmax(allw))

        grid = np.linspace(lo, hi, 512)
        for sign, (a, _, color) in zip((1, -1), groups):
            dens = gaussian_kde(series[a])(grid) * sign
            ax.fill_between(grid, 0, dens, color=color, alpha=0.24, linewidth=0)
            ax.plot(grid, dens, color=color, linewidth=0.95)

        ax.axhline(0, color=INK, linewidth=0.5)
        ax.set_xlim(lo, hi)
        ax.set_xlabel(xlabel, fontsize=7.4)
        ax.set_ylabel("Density", fontsize=7.2)
        ax.set_yticks([])
        ax.spines["left"].set_visible(False)
        fs.panel_label(ax, "ab"[panel_index], dx=-0.06)

        top_label = "de-escalated/stopped"
        bottom_label = "continued broad-spectrum"
        if panel_index == 1:
            ess_treated = _effective_sample_size(series[1])
            ess_cont = _effective_sample_size(series[0])
            top_label += f"  (ESS {int(round(ess_treated)):,})"
            bottom_label += f"  (ESS {int(round(ess_cont)):,})"
            ax.text(
                0.97,
                0.50,
                f"{n_beyond} weights > {hi:.1f} not shown\nmaximum = {max_weight:.2f}",
                transform=ax.transAxes,
                fontsize=6.3,
                color=MUTED,
                ha="right",
                va="center",
                linespacing=1.18,
            )

        yl = ax.get_ylim()
        ax.text(0.97, 0.91, top_label, transform=ax.transAxes,
                fontsize=6.6, color=INK, ha="right", va="top")
        ax.text(0.97, 0.09, bottom_label, transform=ax.transAxes,
                fontsize=6.6, color=INK, ha="right", va="bottom")
        ax.set_ylim(yl)

    fig.subplots_adjust(left=0.06, right=0.98, top=0.88, bottom=0.22, wspace=0.18)
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
