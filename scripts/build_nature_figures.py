#!/usr/bin/env python3
"""Nature-specification rebuild of the submission figures.

Each builder is sized from a Nature column width rather than arbitrary inches,
and each makes one deliberate form change over the current version:

  Fig 2  progressive adjustment -> the four models are nested, so the estimate
         drifting toward the null IS the finding. A connecting trajectory plus
         emphasis on the M4 primary estimate encodes that; four equal-weight
         dots do not. Numbers move to a right-hand column so they stop
         colliding with the interval whiskers.
  Fig 3  cross-dataset outcomes -> 2x2 at 183 mm instead of 1x4 at 325 mm, so
         the figure fits the page it is going on. Dataset identity becomes a
         validated colour pair instead of one blue for everything.
  ESM 1  covariate balance -> before/after per covariate is a dumbbell. Two
         disconnected scatter series make the reader pair up 70 dots by eye;
         a connecting segment shows the movement directly.
  ESM 2  propensity overlap -> mirrored densities. Two overlaid histograms
         with `bins=30` derive separate bin edges per group, so their heights
         are not comparable; sharing one axis removes the defect and the
         overlap region reads at a glance.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

import figstyle as fs
from figstyle import BLUE, FAINT, INK, MUTED, VERMILLION

import matplotlib.pyplot as plt

HARM = Path("outputs/publication_integration/harmonized")
BALANCE = Path("outputs/publication_integration/reviewer_support/mimic_primary_balance_before_after.csv")
BASE_RUN = Path("outputs/mimic/mimic_iv_v5_7_final_20260820T003506Z")
COHORT = BASE_RUN / "audits/vital_repair/analysis_cohort_vital_corrected.csv"
OUT = Path("outputs/publication_integration/nature_figures")

# Frozen publication value for the PSU antibiotic-free-days contrast, applied to
# every figure that shows it so the three packages cannot drift apart again.
SECONDARY_OVERRIDES = {
    ("PSU", "Antibiotic-free days"): {"estimate": 3.16, "ci95_low": 2.85, "ci95_high": 3.47},
}


def apply_overrides(sec: pd.DataFrame) -> pd.DataFrame:
    out = sec.copy()
    for (dataset, outcome), values in SECONDARY_OVERRIDES.items():
        mask = (out["dataset"] == dataset) & (out["outcome"] == outcome)
        if int(mask.sum()) != 1:
            raise RuntimeError(f"Expected exactly one row to override: {dataset}, {outcome}")
        for column, value in values.items():
            out.loc[mask, column] = value
    return out


# --------------------------------------------------------------------------- #
# Figure 2 - progressive adjustment
# --------------------------------------------------------------------------- #

def build_fig2() -> None:
    mort = pd.read_csv(HARM / "harmonized_mortality_results.csv")
    prog = pd.read_csv(HARM / "mimic_progressive_adjustment.csv")
    primary = mort.loc[mort["dataset_analysis"].str.startswith("MIMIC-IV primary")].iloc[0]

    p = prog.copy()
    last = p.index[-1]
    # The manuscript designates the primary 1000-replicate interval for M4. Assert
    # the point estimates agree before grafting one file's interval onto another
    # file's estimate, so a regenerated CSV cannot silently desynchronise them.
    if abs(float(p.loc[last, "risk_difference"]) - float(primary["mortality_rd"])) > 1e-9:
        raise RuntimeError("M4 point estimate disagrees with the harmonized primary result")
    p.loc[last, "rd_lower_95"] = primary["rd_ci95_low"]
    p.loc[last, "rd_upper_95"] = primary["rd_ci95_high"]

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

    # Double column. A forest plot spends its width on three columns - model
    # labels, the plot, the numeric column - so at 120 mm the plotting area is
    # only ~40 mm and the orientation cues below it overlap each other.
    fig, ax = plt.subplots(figsize=(fs.DOUBLE, 2.0))
    fig.subplots_adjust(left=0.21, right=0.72, top=0.90, bottom=0.38)
    fs.null_line(ax)

    # The attenuation trajectory: these models are nested, so the drift from M1
    # to M4 is itself the result. A hairline through the estimates encodes it.
    ax.plot(est, y, color=FAINT, linewidth=0.7, zorder=1, solid_capstyle="round")

    # Emphasis, not categorical colour: M4 is the designated primary estimate,
    # M1-M3 are the context that gets it there.
    colors = [MUTED] * (len(y) - 1) + [BLUE]
    sizes = [3.4] * (len(y) - 1) + [4.6]
    for xi, li, hi_i, yi, c, s in zip(est, lo, hi, y, colors, sizes):
        ax.plot([li, hi_i], [yi, yi], color=c, linewidth=0.9, solid_capstyle="butt", zorder=2)
        for cap in (li, hi_i):
            ax.plot([cap, cap], [yi - 0.13, yi + 0.13], color=c, linewidth=0.9, zorder=2)
        ax.plot([xi], [yi], "o", color=c, markersize=s, zorder=3,
                markeredgecolor="white", markeredgewidth=0.4)

    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_ylim(-0.7, len(y) - 0.3)
    ax.set_xlim(-7.5, 7.5)
    ax.set_xticks([-6, -4, -2, 0, 2, 4, 6])
    ax.set_xlabel("30-day mortality risk difference (percentage points)", labelpad=22)
    fs.strip_y_axis(ax)

    # Numeric column, clear of the whiskers instead of sitting on top of them.
    ax.text(1.06, 1.0, "RD (95% CI)", transform=ax.transAxes, fontsize=6,
            fontweight="bold", va="bottom", ha="left", color=INK)
    for xi, li, hi_i, yi, c in zip(est, lo, hi, y, colors):
        ax.text(1.06, yi, f"{xi:+.2f} ({li:+.2f}, {hi_i:+.2f})",
                transform=ax.get_yaxis_transform(), fontsize=6,
                va="center", ha="left", color=INK if c == BLUE else MUTED,
                clip_on=False)

    # Orientation cues get their own band below the tick labels; inside the axes
    # they collide with the M4 interval at this column width.
    ax.text(0.0, -0.30, "← favours de-escalation", transform=ax.transAxes,
            fontsize=5.5, color=MUTED, ha="left", va="center")
    ax.text(1.0, -0.30, "favours continuation →", transform=ax.transAxes,
            fontsize=5.5, color=MUTED, ha="right", va="center")

    fs.savefig(fig, OUT, "Fig2_progressive_adjustment")


# --------------------------------------------------------------------------- #
# Figure 3 - cross-dataset outcomes
# --------------------------------------------------------------------------- #

def build_fig3() -> None:
    mort = pd.read_csv(HARM / "harmonized_mortality_results.csv")
    sec = apply_overrides(pd.read_csv(HARM / "harmonized_secondary_outcomes.csv"))
    mm = mort.loc[mort["dataset_analysis"].str.startswith("MIMIC-IV primary")].iloc[0]
    pm = mort.loc[mort["dataset_analysis"].str.startswith("PSU modified")].iloc[0]

    # No panel titles: Nature carries the panel description in the caption, so the
    # axis label states the quantity in full and the panel letter does the rest.
    panels = [
        ("mortality", "30-day mortality risk difference (percentage points)", 2),
        ("Antibiotic-free days", "Antibiotic-free days (difference in days)", 2),
        ("Normalized systemic antibiotic exposure", "Systemic antibiotic exposure (difference in proportion)", 3),
        ("Normalized broad-spectrum exposure", "Broad-spectrum exposure (difference in proportion)", 3),
    ]
    n_lab = [f"MIMIC-IV\nn = {int(mm['cohort_n']):,}", f"Penn State\nn = {int(pm['cohort_n']):,}"]

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
            ax.plot([li, hi_i], [yi, yi], color=c, linewidth=1.0, solid_capstyle="butt", zorder=2)
            for cap in (li, hi_i):
                ax.plot([cap, cap], [yi - 0.1, yi + 0.1], color=c, linewidth=1.0, zorder=2)
            ax.plot([xi], [yi], "o", color=c, markersize=4.2, zorder=3,
                    markeredgecolor="white", markeredgewidth=0.4)

        # Pad the axis first, then place labels in data space with a margin that
        # is a fixed fraction of the FINAL span - so no label can land on a whisker.
        span = max(his.max(), 0) - min(los.min(), 0)
        pad = 0.18 * span
        ax.set_xlim(min(los.min(), 0) - pad, max(his.max(), 0) + pad)

        fmt = f"{{:+.{decimals}f}}"
        for xi, hi_i, yi, c in zip(vals, his, yy, colors):
            ax.text(hi_i + 0.04 * span, yi, fmt.format(xi), fontsize=6,
                    va="center", ha="left", color=c)

        ax.set_yticks(yy)
        # Cohort size rides with the dataset label rather than floating in-panel;
        # the two arms have different n, so one corner annotation would be wrong.
        ax.set_yticklabels(n_lab, fontsize=6)
        ax.set_ylim(-0.55, 1.55)
        ax.set_xlabel(unit, fontsize=6.5)
        fs.strip_y_axis(ax)
        fs.panel_label(ax, "abcd"[j], dx=-0.26)

    fig.subplots_adjust(left=0.11, right=0.98, top=0.92, bottom=0.13,
                        hspace=1.05, wspace=0.42)
    fs.savefig(fig, OUT, "Fig3_cross_dataset_outcomes")


# --------------------------------------------------------------------------- #
# ESM Figure 1 - covariate balance, as a dumbbell
# --------------------------------------------------------------------------- #

PRETTY = {
    "broad_abx_hours_pre72": "Broad-spectrum antibiotic hours, pre-72 h",
    "systemic_abx_hours_pre72": "Systemic antibiotic hours, pre-72 h",
    "antipseudomonal_pre72": "Antipseudomonal therapy, pre-72 h",
    "broad_abx_agents_pre72": "Broad-spectrum agents, pre-72 h",
    "micro_records_pre72": "Microbiology records, pre-72 h",
    "anaerobic_coverage_pre72": "Anaerobic coverage, pre-72 h",
    "strict_culture_records_pre72": "Culture records, pre-72 h",
    "carbapenem_pre72": "Carbapenem use, pre-72 h",
    "distinct_specimen_types_pre72": "Distinct specimen types, pre-72 h",
    "cardiac_icu": "Cardiac ICU",
    "respiratory_culture_pre72": "Respiratory culture, pre-72 h",
    "temperature_48_72h": "Temperature, 48-72 h",
    "vent_proc": "Mechanical ventilation",
    "sterile_fluid_culture_pre72": "Sterile-fluid culture, pre-72 h",
    "sofa_like_change_pre72": "SOFA-like change, pre-72 h",
    "sofa_like_48_72h": "SOFA-like score, 48-72 h",
    "fever_last12h_pre72": "Fever in prior 12 h, pre-72 h",
    "repeat_micro_48_72h": "Repeat microbiology, 48-72 h",
    "blood_culture_pre72": "Blood culture, pre-72 h",
    "hr_max_pre72": "Maximum heart rate, pre-72 h",
    "platelet_late_worst_48_72h": "Platelet count, worst 48-72 h",
    "sofa_like_improved_pre72": "SOFA-like improvement, pre-72 h",
    "lactate_last_pre72": "Last lactate, pre-72 h",
    "bilirubin_rising_pre72": "Rising bilirubin, pre-72 h",
    "wbc_late_last_48_72h": "Last WBC count, 48-72 h",
    "sicu": "Surgical ICU",
    "systemic_abx_agents_pre72": "Systemic antibiotic agents, pre-72 h",
    "lactate_rising_pre72": "Rising lactate, pre-72 h",
    "hours_admit_to_icu": "Admission-to-ICU interval, h",
    "micu": "Medical ICU",
    "white_blood_cells_last_pre72": "Last WBC count, pre-72 h",
    "wbc_rising_pre72": "Rising WBC count, pre-72 h",
    "severity_pre72": "Severity index, pre-72 h",
    "bilirubin_late_worst_48_72h": "Bilirubin, worst 48-72 h",
    "vasopressor_stopped_before_72h": "Vasopressor stopped before 72 h",
}

ABBREV = {"icu": "ICU", "wbc": "WBC", "sofa": "SOFA", "iv": "IV",
          "bmi": "BMI", "spo2": "SpO2", "map": "MAP", "rr": "RR", "hr": "HR"}


def pretty_label(name: str) -> str:
    if name in PRETTY:
        return PRETTY[name]
    words = [ABBREV.get(w.lower(), w) for w in str(name).replace("_", " ").split()]
    text = " ".join(words)
    return text[:1].upper() + text[1:]


def build_esm1() -> None:
    bal = pd.read_csv(BALANCE)
    top = bal.sort_values("before", ascending=False).head(35).sort_values("before")
    before = top["before"].to_numpy(float)
    after = top["after"].to_numpy(float)
    y = np.arange(len(top))

    # 35 covariate names need roughly half the page width; the depth stays under
    # the 247 mm Nature ceiling at 35 rows of ~4 mm.
    fig, ax = plt.subplots(figsize=(fs.ONE_HALF, 5.9))
    fig.subplots_adjust(left=0.46, right=0.98, top=0.97, bottom=0.08)

    # The 0.1 conventional balance threshold, as recessive non-data ink.
    ax.axvline(0.1, color=fs.RULE, linestyle=(0, (3, 2)), linewidth=0.5, zorder=0)

    # The connecting segment is the point of a dumbbell: it shows each covariate's
    # movement, which two independent scatter series leave the reader to infer.
    for b, a, yi in zip(before, after, y):
        ax.plot([b, a], [yi, yi], color=FAINT, linewidth=0.7, zorder=1,
                solid_capstyle="round")
    ax.plot(before, y, "o", color=FAINT, markersize=3.0, zorder=2,
            markeredgecolor="white", markeredgewidth=0.3, label="Before weighting")
    ax.plot(after, y, "o", color=BLUE, markersize=3.4, zorder=3,
            markeredgecolor="white", markeredgewidth=0.3, label="After weighting")

    ax.set_yticks(y)
    ax.set_yticklabels([pretty_label(v) for v in top["variable"]], fontsize=5.6)
    ax.set_ylim(-0.8, len(y) - 0.2)
    ax.set_xlim(left=0)
    ax.set_xlabel("Absolute standardized mean difference")
    fs.strip_y_axis(ax)
    # Threshold label sits at the foot of its own rule; at the top it collides with
    # the direct labels, since the best-balanced covariates cluster near 0.1.
    ax.text(0.1, -0.62, " 0.10 balance threshold", fontsize=5.5, color=MUTED,
            va="center", ha="left")
    # Direct labels on the top row's two markers instead of a legend box floating
    # in the empty lower-right region.
    ax.annotate("before weighting", xy=(before[-1], y[-1]), xytext=(0, 9),
                textcoords="offset points", fontsize=5.5, color=MUTED, ha="center")
    ax.annotate("after weighting", xy=(after[-1], y[-1]), xytext=(0, 9),
                textcoords="offset points", fontsize=5.5, color=BLUE, ha="center")
    fs.savefig(fig, OUT, "ESM_Fig1_covariate_balance")


# --------------------------------------------------------------------------- #
# ESM Figure 2 - propensity overlap and weights, as mirrored densities
# --------------------------------------------------------------------------- #

def build_esm2() -> None:
    from scipy.stats import gaussian_kde

    from sepsis_deescalation.specification import CANDIDATE_PS_VARS
    from sepsis_deescalation.stats import fit_stabilized_iptw

    d = pd.read_csv(COHORT, low_memory=False)
    w, _, _ = fit_stabilized_iptw(d, CANDIDATE_PS_VARS)

    fig, axes = plt.subplots(1, 2, figsize=(fs.DOUBLE, 2.05))
    groups = [(1, "De-escalated or stopped", BLUE), (0, "Continued broad-spectrum", VERMILLION)]

    for ax, (var, xlabel, lo, hi) in zip(
        axes,
        [("ps_den", "Estimated propensity for de-escalation", 0.0, 1.0),
         ("SW_A", "Stabilized IPTW", 0.0, None)],
    ):
        series = {a: pd.to_numeric(w.loc[w["A"] == a, var], errors="coerce").dropna().to_numpy()
                  for a, _, _ in groups}
        if hi is None:
            # A handful of extreme stabilized weights run to ~30. Scaling the axis
            # to the maximum compresses 99.5% of the mass into an unreadable spike,
            # so clip the view and state the tail in text instead of drawing it -
            # the extreme weights are a diagnostic fact, not a shape worth plotting.
            allw = np.concatenate(list(series.values()))
            hi = float(np.percentile(allw, 99.5))
            n_beyond = int((allw > hi).sum())
            ax.text(
                0.98, 0.5,
                f"{n_beyond} weights > {hi:.1f} not shown (max {allw.max():.1f})",
                transform=ax.transAxes, fontsize=5.2, color=MUTED,
                ha="right", va="center",
            )
        grid = np.linspace(lo, hi, 512)

        # One shared evaluation grid for both groups - the defect in the current
        # figure is that `bins=30` gives each group its own edges, so the two
        # curves are drawn against different x-quantisations and their heights
        # cannot be compared. Mirroring also removes the occlusion entirely.
        for sign, (a, label, color) in zip((1, -1), groups):
            dens = gaussian_kde(series[a])(grid) * sign
            ax.fill_between(grid, 0, dens, color=color, alpha=0.30, linewidth=0)
            ax.plot(grid, dens, color=color, linewidth=0.9)

        ax.axhline(0, color=INK, linewidth=0.5)
        ax.set_xlim(lo, hi)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Density")
        ax.set_yticks([])
        ax.spines["left"].set_visible(False)

    fs.panel_label(axes[0], "a", dx=-0.06)
    fs.panel_label(axes[1], "b", dx=-0.06)

    for ax in axes:
        yl = ax.get_ylim()
        ax.text(0.98, 0.92, "de-escalated", transform=ax.transAxes, fontsize=5.5,
                color=BLUE, ha="right", va="top")
        ax.text(0.98, 0.08, "continued", transform=ax.transAxes, fontsize=5.5,
                color=VERMILLION, ha="right", va="bottom")
        ax.set_ylim(yl)

    fig.subplots_adjust(left=0.06, right=0.98, top=0.88, bottom=0.22, wspace=0.16)
    fs.savefig(fig, OUT, "ESM_Fig2_propensity_overlap")


# --------------------------------------------------------------------------- #
# Figure 1b - cohort attrition cascade
# --------------------------------------------------------------------------- #

FLOW = Path("outputs/publication_integration/reviewer_support/mimic_cohort_flow_reviewer.csv")

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


def build_fig1b() -> None:
    """Two-part cascade: sequential exclusions, then the treatment split.

    The exclusions are strictly nested and never rejoin, so this is a waterfall,
    not an alluvial - a flow diagram would add ribbons that encode nothing. Each
    row is drawn to the PREVIOUS stage's width, so the pale remainder is the
    number excluded at that step and the loss is visible rather than annotated.
    """
    f = pd.read_csv(FLOW)
    stages = f.loc[f["stage"].isin(SHORT_STAGE)].copy()
    arms = f.loc[f["stage"].isin(["De-escalated/stopped", "Continued broad-spectrum"])]

    n = stages["n"].to_numpy(float)
    prev = np.concatenate([[n[0]], n[:-1]])
    y = np.arange(len(stages))[::-1]

    fig, (ax, ax2) = plt.subplots(
        2, 1, figsize=(fs.DOUBLE, 2.5),
        gridspec_kw={"height_ratios": [len(stages), 1.0]},
    )
    fig.subplots_adjust(left=0.30, right=0.97, top=0.95, bottom=0.10, hspace=0.30)

    # Counts live in two fixed text columns to the right of the longest bar. Placing
    # them at the bar ends instead collides as soon as consecutive stages are close
    # in size - which they are from the third exclusion onward.
    col_n = n[0] * 1.04
    col_rm = n[0] * 1.34
    ax.barh(y, prev, color=FAINT, alpha=0.35, height=0.62, linewidth=0)
    ax.barh(y, n, color=MUTED, height=0.62, linewidth=0)
    for ni, pi, yi in zip(n, prev, y):
        ax.text(col_n, yi, f"{int(ni):,}", va="center", ha="left", fontsize=6, color=INK)
        if pi > ni:
            ax.text(col_rm, yi, f"−{int(pi - ni):,}", va="center", ha="right",
                    fontsize=5.8, color=MUTED)
    ax.text(col_n, len(stages) - 0.45, "Retained", fontsize=5.8, fontweight="bold",
            ha="left", va="center", color=INK)
    ax.text(col_rm, len(stages) - 0.45, "Excluded", fontsize=5.8, fontweight="bold",
            ha="right", va="center", color=MUTED)
    ax.set_yticks(y)
    ax.set_yticklabels([SHORT_STAGE[s] for s in stages["stage"]], fontsize=6, linespacing=1.15)
    ax.set_xlim(0, col_rm * 1.02)
    ax.set_ylim(-0.7, len(stages) - 0.15)
    ax.set_xticks([])
    ax.spines["bottom"].set_visible(False)
    fs.strip_y_axis(ax)
    fs.panel_label(ax, "a", dx=-0.42, dy=0.99)

    arm_n = arms.set_index("stage").loc[
        ["De-escalated/stopped", "Continued broad-spectrum"], "n"].to_numpy(float)
    total = arm_n.sum()
    # Panel b gets its own scale. Sharing panel a's axis would render the analytic
    # cohort as a 10%-wide stub and waste the row.
    ax2.barh([0], [arm_n[0]], color=BLUE, height=0.42, linewidth=0)
    ax2.barh([0], [arm_n[1]], left=[arm_n[0]], color=VERMILLION, height=0.42, linewidth=0)
    for centre, val, color, label in [
        (arm_n[0] / 2, arm_n[0], BLUE, "De-escalated or stopped"),
        (arm_n[0] + arm_n[1] / 2, arm_n[1], VERMILLION, "Continued broad-spectrum"),
    ]:
        ax2.text(centre, -0.42, f"{label}\n{int(val):,} ({100*val/total:.1f}%)",
                 ha="center", va="top", fontsize=5.8, color=color, linespacing=1.3)
    ax2.set_xlim(0, total * 1.34 / 1.04)
    ax2.set_ylim(-1.45, 0.45)
    ax2.set_xticks([])
    ax2.set_yticks([])
    for s in ("bottom", "left"):
        ax2.spines[s].set_visible(False)
    ax2.text(-0.01, 0, f"Analytic cohort\nn = {int(total):,}", transform=ax2.get_yaxis_transform(),
             fontsize=6, ha="right", va="center", color=INK, linespacing=1.2)
    fs.panel_label(ax2, "b", dx=-0.42, dy=0.75)

    fs.savefig(fig, OUT, "Fig1b_cohort_attrition")


def main() -> None:
    fs.use_nature_style()
    OUT.mkdir(parents=True, exist_ok=True)

    build_fig2()
    build_fig3()
    print("Built Fig2, Fig3")

    if FLOW.exists():
        build_fig1b()
        print("Built Fig1b attrition cascade")
    else:
        print(f"skip Fig1b, missing {FLOW}")

    if BALANCE.exists():
        build_esm1()
        print("Built ESM Fig1")
    else:
        print(f"skip ESM Fig1, missing {BALANCE}")

    if COHORT.exists():
        build_esm2()
        print("Built ESM Fig2")
    else:
        print(f"skip ESM Fig2, missing {COHORT}")


if __name__ == "__main__":
    main()
