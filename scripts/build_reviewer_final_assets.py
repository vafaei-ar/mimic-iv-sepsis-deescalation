#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HARM = Path("outputs/publication_integration/harmonized")
OUT = Path("outputs/publication_integration/reviewer_support")


def savefig(fig: plt.Figure, stem: str) -> None:
    fig.tight_layout()
    fig.savefig(OUT / f"{stem}.png", dpi=400, bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)


def build_fig1() -> None:
    fig, ax = plt.subplots(figsize=(9.2, 2.8))
    ax.set_xlim(-6, 126)
    ax.set_ylim(-1.15, 1.22)
    ax.axis("off")
    ax.annotate("", xy=(121, 0), xytext=(0, 0), arrowprops=dict(arrowstyle="->", linewidth=2.0))

    markers = [
        (0, "t0", "First qualifying broad-spectrum\nantibiotic", 0.66),
        (72, "72 h", "Treatment\ndecision", 0.66),
        (96, "96 h", "Landmark;\nfollow-up starts", 0.66),
        (121, "30 d", "Post-landmark\noutcome horizon", 0.66),
    ]
    for x, top, label, ylab in markers:
        ax.plot([x, x], [-0.13, 0.13], linewidth=1.8)
        ax.text(x, 0.32, top, ha="center", va="bottom", fontsize=11, fontweight="bold")
        ax.text(x, ylab, label, ha="center", va="bottom", fontsize=8.6, linespacing=1.05)

    ax.plot([0, 72], [-0.40, -0.40], linewidth=6, alpha=.45)
    ax.text(36, -0.63, "Pre-decision covariates", ha="center", fontsize=9.4)
    ax.plot([72, 96], [-0.40, -0.40], linewidth=6, alpha=.72)
    ax.text(84, -0.63, "Treatment classification", ha="center", fontsize=9.4)
    ax.plot([96, 121], [-0.40, -0.40], linewidth=6)
    ax.text(108.5, -0.63, "Outcome follow-up", ha="center", fontsize=9.4)
    savefig(fig, "Fig1_target_trial_timeline_final")


def build_fig2() -> None:
    mort = pd.read_csv(HARM / "harmonized_mortality_results.csv")
    prog = pd.read_csv(HARM / "mimic_progressive_adjustment.csv")
    primary = mort.loc[mort["dataset_analysis"].str.startswith("MIMIC-IV primary")].iloc[0]
    p = prog.copy()
    last = p.index[-1]
    p.loc[last, "rd_lower_95"] = primary["rd_ci95_low"]
    p.loc[last, "rd_upper_95"] = primary["rd_ci95_high"]
    labels = [
        "M1  Demographics/comorbidity",
        "M2  + baseline severity",
        "M3  + day-3 clinical status",
        "M4  + trajectories/intensity",
    ]
    est = 100 * p["risk_difference"].to_numpy(float)
    lo = 100 * p["rd_lower_95"].to_numpy(float)
    hi = 100 * p["rd_upper_95"].to_numpy(float)
    y = np.arange(len(p))[::-1]
    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    ax.errorbar(est, y, xerr=[est-lo, hi-est], fmt="o", capsize=3, linewidth=1.5)
    ax.axvline(0, linestyle="--", linewidth=1)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xlabel("30-day post-landmark mortality risk difference, percentage points")
    ax.set_title("Progressive adjustment in MIMIC-IV")
    for x, yy in zip(est, y):
        ax.text(x + 0.22, yy, f"{x:+.2f}", va="center", fontsize=8.5)
    savefig(fig, "Fig2_progressive_adjustment_final")


def build_fig3() -> None:
    mort = pd.read_csv(HARM / "harmonized_mortality_results.csv")
    sec = pd.read_csv(HARM / "harmonized_secondary_outcomes.csv")
    mm = mort.loc[mort["dataset_analysis"].str.startswith("MIMIC-IV primary")].iloc[0]
    pm = mort.loc[mort["dataset_analysis"].str.startswith("PSU modified")].iloc[0]

    panels = [
        ("mortality", "Mortality RD", "percentage points"),
        ("Antibiotic-free days", "Antibiotic-free days", "days"),
        ("Normalized systemic antibiotic exposure", "Systemic exposure", "proportion"),
        ("Normalized broad-spectrum exposure", "Broad-spectrum exposure", "proportion"),
    ]
    fig, axes = plt.subplots(1, 4, figsize=(12.8, 3.7))
    for j, (key, title, unit) in enumerate(panels):
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
        ax.set_ylim(-0.35, 1.35)
        ax.set_yticks(yy)
        ax.set_yticklabels(["MIMIC-IV", "Penn State"] if j == 0 else [])
        ax.set_title(title, fontsize=9.5)
        ax.set_xlabel(unit, fontsize=8.5)
        ax.tick_params(axis="x", labelsize=8)
        for x, y in zip(vals, yy):
            span = max(abs(ax.get_xlim()[1] - ax.get_xlim()[0]), 1e-6)
            ax.text(x + 0.025*span, y, f"{x:+.2f}" if key in {"mortality", "Antibiotic-free days"} else f"{x:+.3f}", va="center", fontsize=7.8)
    fig.suptitle("Cross-dataset mortality and stewardship outcomes", fontsize=11, y=1.02)
    savefig(fig, "Fig3_cross_dataset_outcomes_final")


def combine_subtypes() -> pd.DataFrame:
    files = [
        OUT / "mimic_deescalation_subtype_narrowed_adjusted_mortality.csv",
        OUT / "mimic_deescalation_subtype_stopped_adjusted_mortality.csv",
    ]
    missing = [str(p) for p in files if not p.exists()]
    if missing:
        raise FileNotFoundError("Missing subtype outputs: " + ", ".join(missing))
    out = pd.concat([pd.read_csv(p) for p in files], ignore_index=True)
    out.to_csv(OUT / "mimic_deescalation_subtype_adjusted_mortality_combined.csv", index=False)
    return out


def build_summary(sub: pd.DataFrame) -> None:
    flow = pd.read_csv(OUT / "mimic_cohort_flow_reviewer.csv")
    desc = pd.read_csv(OUT / "mimic_deescalation_subtype_descriptives.csv")
    narrowed = desc.loc[desc["deescalation_type"] == "narrowed_or_non_broad_only"].iloc[0]
    stopped = desc.loc[desc["deescalation_type"] == "stopped_all_observed_systemic_antibiotics"].iloc[0]
    text = f"""# Final reviewer-support package\n\nThe frozen primary MIMIC-IV and Penn State analyses are unchanged. New treatment-subtype analyses are exploratory only.\n\n## De-escalation composition\n- Narrowed/non-broad only: {int(narrowed['n']):,} ({narrowed['percent_within_deescalated']:.1f}% of de-escalated/stopped).\n- Stopped all observed systemic antibiotics: {int(stopped['n']):,} ({stopped['percent_within_deescalated']:.1f}% of de-escalated/stopped).\n\n## Exploratory subtype mortality sensitivities\n{sub.to_markdown(index=False)}\n\nThe complete-stopping contrast has materially poorer overlap and a very low treated effective sample size, so it should be interpreted cautiously and retained in supplementary material rather than used as a separate causal claim.\n\n## Cohort flow\n{flow.to_markdown(index=False)}\n"""
    (OUT / "reviewer_support_final_summary.md").write_text(text, encoding="utf-8")
    metadata = {
        "primary_science_changed": False,
        "subtype_analyses_exploratory": True,
        "fig2_m4_ci_source": "designated manuscript primary 1000-replicate mortality CI",
    }
    (OUT / "reviewer_support_final_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for p in [
        HARM / "harmonized_mortality_results.csv",
        HARM / "harmonized_secondary_outcomes.csv",
        HARM / "mimic_progressive_adjustment.csv",
        OUT / "mimic_cohort_flow_reviewer.csv",
        OUT / "mimic_deescalation_subtype_descriptives.csv",
    ]:
        if not p.exists():
            raise FileNotFoundError(p)
    sub = combine_subtypes()
    build_fig1()
    build_fig2()
    build_fig3()
    build_summary(sub)


if __name__ == "__main__":
    main()
