#!/usr/bin/env python3
"""Build the ICM submission figures from the preferred frozen-data figure system.

This wrapper reuses the publication-locked Figure 1, Figure 2, and ESM diagnostics
from ``build_nature_figures.py`` while giving Figure 3 a more clinically focused
editorial role. No inferential analysis is rerun and no scientific value is changed.
"""
from __future__ import annotations

from pathlib import Path
import shutil

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

import build_nature_figures as base
from figstyle import BLUE, INK, MUTED, VERMILLION
from publication_figure_common import apply_publication_secondary_overrides

OUT = Path("outputs/publication_integration/submission_figures")


def _copy_alias(stem_from: str, stem_to: str) -> None:
    for ext in ("png", "pdf"):
        src = OUT / f"{stem_from}.{ext}"
        dst = OUT / f"{stem_to}.{ext}"
        if not src.exists():
            raise FileNotFoundError(src)
        shutil.copyfile(src, dst)


def build_fig3() -> None:
    """Build a clinically prioritized cross-dataset comparison.

    The figure gives the primary mortality result and hospital-free days the
    clinical-outcome section, then shows antibiotic-free days as the single
    stewardship outcome. Normalized exposure measures remain available in the
    publication tables/ESM but are not given equal visual weight in the main figure.
    """
    mort = pd.read_csv(base.HARM / "harmonized_mortality_results.csv")
    sec = apply_publication_secondary_overrides(
        pd.read_csv(base.HARM / "harmonized_secondary_outcomes.csv")
    )
    mm = mort.loc[mort["dataset_analysis"].str.startswith("MIMIC-IV primary")].iloc[0]
    pm = mort.loc[mort["dataset_analysis"].str.startswith("PSU modified")].iloc[0]

    rows = [
        {
            "title": "30-day post-landmark mortality",
            "xlabel": "Risk difference (percentage points)",
            "decimals": 2,
            "vals": np.array([100 * mm["mortality_rd"], 100 * pm["mortality_rd"]], float),
            "los": np.array([100 * mm["rd_ci95_low"], 100 * pm["rd_ci95_low"]], float),
            "his": np.array([100 * mm["rd_ci95_high"], 100 * pm["rd_ci95_high"]], float),
        },
    ]

    for outcome, title, xlabel in [
        ("Hospital-free days", "Hospital-free days", "Difference in days"),
        ("Antibiotic-free days", "Antibiotic-free days", "Difference in days"),
    ]:
        s = sec.loc[sec["outcome"] == outcome].set_index("dataset").loc[["MIMIC-IV", "PSU"]]
        rows.append(
            {
                "title": title,
                "xlabel": xlabel,
                "decimals": 2,
                "vals": s["estimate"].to_numpy(float),
                "los": s["ci95_low"].to_numpy(float),
                "his": s["ci95_high"].to_numpy(float),
            }
        )

    fig = plt.figure(figsize=(base.fs.DOUBLE, 4.65))
    gs = fig.add_gridspec(
        5,
        1,
        height_ratios=[0.18, 1.0, 1.0, 0.22, 1.0],
        hspace=1.42,
    )

    head_clinical = fig.add_subplot(gs[0])
    head_clinical.axis("off")
    head_clinical.text(
        0, 0.55, "Clinical outcomes", fontsize=7.2, fontweight="bold",
        ha="left", va="center", color=INK,
    )
    axes = [fig.add_subplot(gs[1]), fig.add_subplot(gs[2])]

    head_stewardship = fig.add_subplot(gs[3])
    head_stewardship.axis("off")
    head_stewardship.text(
        0, 0.55, "Stewardship outcome", fontsize=7.2, fontweight="bold",
        ha="left", va="center", color=INK,
    )
    axes.append(fig.add_subplot(gs[4]))

    for panel_index, (ax, row) in enumerate(zip(axes, rows)):
        base._draw_two_dataset_effect_row(ax, **row)
        base.fs.panel_label(ax, "abc"[panel_index], dx=-0.17, dy=1.08)

    handles = [
        Line2D(
            [0], [0], marker="o", color="none", markerfacecolor=BLUE,
            markeredgecolor="white", markeredgewidth=0.4, markersize=5,
            label=f"MIMIC-IV (n = {int(mm['cohort_n']):,})",
        ),
        Line2D(
            [0], [0], marker="o", color="none", markerfacecolor=VERMILLION,
            markeredgecolor="white", markeredgewidth=0.4, markersize=5,
            label=f"Penn State (n = {int(pm['cohort_n']):,})",
        ),
    ]
    fig.legend(
        handles=handles,
        loc="upper right",
        bbox_to_anchor=(0.97, 0.988),
        ncol=2,
        frameon=False,
        handletextpad=0.35,
        columnspacing=1.2,
        fontsize=6.3,
    )
    fig.text(
        0.77,
        0.947,
        "Secondary modified replication; estimates are not pooled",
        ha="right",
        va="top",
        fontsize=5.8,
        color=MUTED,
    )
    fig.subplots_adjust(left=0.22, right=0.77, top=0.935, bottom=0.085)
    base.fs.savefig(fig, OUT, "Fig3_cross_dataset_outcomes")


def main() -> None:
    base.fs.use_nature_style()
    OUT.mkdir(parents=True, exist_ok=True)
    base.OUT = OUT

    required = [
        base.HARM / "harmonized_mortality_results.csv",
        base.HARM / "harmonized_secondary_outcomes.csv",
        base.HARM / "mimic_progressive_adjustment.csv",
        base.FLOW,
        base.BALANCE,
        base.COHORT,
    ]
    missing = [str(p) for p in required if not p.exists()]
    if missing:
        raise FileNotFoundError("Missing required inputs: " + ", ".join(missing))

    base.build_fig1()
    base.build_fig2()
    build_fig3()
    base.build_esm1()
    base.build_esm2()

    _copy_alias("Fig1_target_trial_and_cohort", "Fig1_target_trial_timeline")
    _copy_alias("ESM_Fig1_covariate_balance", "ESM_Fig1_mimic_balance_love")
    _copy_alias("ESM_Fig2_propensity_overlap", "ESM_Fig2_mimic_ps_weights")
    print("Built ICM submission Figure 1-3 and ESM Figure 1-2 from frozen inputs")


if __name__ == "__main__":
    main()
