#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from sepsis_deescalation.specification import CANDIDATE_PS_VARS
from sepsis_deescalation.stats import fit_stabilized_iptw

BASE_RUN = Path("outputs/mimic/mimic_iv_v5_7_final_20260820T003506Z")
COHORT = BASE_RUN / "audits/vital_repair/analysis_cohort_vital_corrected.csv"
HARM = Path("outputs/publication_integration/harmonized")
OUT = Path("outputs/publication_integration/submission_figures")
BALANCE = Path("outputs/publication_integration/reviewer_support/mimic_primary_balance_before_after.csv")


def savefig(fig: plt.Figure, stem: str) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
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
        (0, "t0", "First qualifying broad-spectrum\nantibiotic"),
        (72, "72 h", "Treatment\ndecision"),
        (96, "96 h", "Landmark;\nfollow-up starts"),
        (121, "30 d", "Post-landmark\noutcome horizon"),
    ]
    for x, top, label in markers:
        ax.plot([x, x], [-0.13, 0.13], linewidth=1.8)
        ax.text(x, 0.32, top, ha="center", va="bottom", fontsize=11, fontweight="bold")
        ax.text(x, 0.66, label, ha="center", va="bottom", fontsize=8.6, linespacing=1.05)
    ax.plot([0, 72], [-0.40, -0.40], linewidth=6, alpha=0.45)
    ax.text(36, -0.63, "Pre-decision covariates", ha="center", fontsize=9.4)
    ax.plot([72, 96], [-0.40, -0.40], linewidth=6, alpha=0.72)
    ax.text(84, -0.63, "Treatment classification", ha="center", fontsize=9.4)
    ax.plot([96, 121], [-0.40, -0.40], linewidth=6)
    ax.text(108.5, -0.63, "Outcome follow-up", ha="center", fontsize=9.4)
    savefig(fig, "Fig1_target_trial_timeline")


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

    fig, ax = plt.subplots(figsize=(7.8, 4.6))
    ax.errorbar(est, y, xerr=[est - lo, hi - est], fmt="o", capsize=3, linewidth=1.6)
    ax.axvline(0, linestyle="--", linewidth=1)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_ylim(-0.45, 3.55)
    ax.set_xlabel("30-day post-landmark mortality risk difference, percentage points")
    ax.set_title("Progressive adjustment in MIMIC-IV")

    # Keep estimate and CI text clear of the confidence-interval line.
    for x, l, h, yy in zip(est, lo, hi, y):
        ax.annotate(
            f"{x:+.2f} ({l:+.2f}, {h:+.2f})",
            xy=(x, yy),
            xytext=(0, 9),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=8.4,
        )
    savefig(fig, "Fig2_progressive_adjustment")


def build_fig3() -> None:
    mort = pd.read_csv(HARM / "harmonized_mortality_results.csv")
    sec = pd.read_csv(HARM / "harmonized_secondary_outcomes.csv")
    mm = mort.loc[mort["dataset_analysis"].str.startswith("MIMIC-IV primary")].iloc[0]
    pm = mort.loc[mort["dataset_analysis"].str.startswith("PSU modified")].iloc[0]

    panels = [
        ("mortality", "Mortality risk difference", "percentage points"),
        ("Antibiotic-free days", "Antibiotic-free days", "days"),
        ("Normalized systemic antibiotic exposure", "Systemic antibiotic exposure", "proportion"),
        ("Normalized broad-spectrum exposure", "Broad-spectrum antibiotic exposure", "proportion"),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(9.2, 6.9))
    axes = axes.ravel()
    for j, (key, title, unit) in enumerate(panels):
        ax = axes[j]
        if key == "mortality":
            vals = np.array([100 * mm["mortality_rd"], 100 * pm["mortality_rd"]], dtype=float)
            los = np.array([100 * mm["rd_ci95_low"], 100 * pm["rd_ci95_low"]], dtype=float)
            his = np.array([100 * mm["rd_ci95_high"], 100 * pm["rd_ci95_high"]], dtype=float)
            decimals = 2
        else:
            s = sec.loc[sec["outcome"] == key].set_index("dataset").loc[["MIMIC-IV", "PSU"]]
            vals = s["estimate"].to_numpy(float)
            los = s["ci95_low"].to_numpy(float)
            his = s["ci95_high"].to_numpy(float)
            decimals = 2 if key == "Antibiotic-free days" else 3

        yy = np.array([1, 0])
        ax.errorbar(vals, yy, xerr=[vals - los, his - vals], fmt="o", capsize=3, linewidth=1.5)
        ax.axvline(0, linestyle="--", linewidth=1)
        ax.set_ylim(-0.45, 1.55)
        ax.set_yticks(yy)
        ax.set_yticklabels(["MIMIC-IV", "Penn State"])
        ax.set_title(title, fontsize=10.2)
        ax.set_xlabel(unit, fontsize=9.0)
        ax.tick_params(axis="both", labelsize=8.6)

        for x, l, h, y0 in zip(vals, los, his, yy):
            txt = (
                f"{x:+.2f} ({l:+.2f}, {h:+.2f})"
                if decimals == 2
                else f"{x:+.3f} ({l:+.3f}, {h:+.3f})"
            )
            ax.annotate(
                txt,
                xy=(x, y0),
                xytext=(0, 9),
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=8.1,
            )

    fig.suptitle("Cross-dataset mortality and stewardship outcomes", fontsize=11.5, y=1.01)
    savefig(fig, "Fig3_cross_dataset_outcomes")


PRETTY_LABELS = {
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
    "vent_proc": "Mechanical ventilation procedure",
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


def pretty_label(name: str) -> str:
    if name in PRETTY_LABELS:
        return PRETTY_LABELS[name]
    text = str(name).replace("_", " ")
    replacements = {
        "icu": "ICU",
        "wbc": "WBC",
        "sofa": "SOFA",
        "iv": "IV",
        "bmi": "BMI",
        "spo2": "SpO2",
        "map": "MAP",
        "rr": "RR",
        "hr": "HR",
    }
    words = [replacements.get(w.lower(), w) for w in text.split()]
    text = " ".join(words)
    return text[:1].upper() + text[1:]


def build_esm1() -> None:
    bal = pd.read_csv(BALANCE)
    top = bal.sort_values("before", ascending=False).head(35).sort_values("before")
    fig, ax = plt.subplots(figsize=(8.4, 9.4))
    y = np.arange(len(top))
    ax.scatter(top["before"], y, marker="o", label="Before weighting")
    ax.scatter(top["after"], y, marker="s", label="After weighting")
    ax.axvline(0.1, linestyle="--", linewidth=1)
    ax.set_yticks(y)
    ax.set_yticklabels([pretty_label(v) for v in top["variable"]], fontsize=8.3)
    ax.set_xlabel("Absolute standardized mean difference")
    ax.set_title("MIMIC-IV primary covariate balance")
    ax.legend(frameon=False, fontsize=10)
    savefig(fig, "ESM_Fig1_mimic_balance_love")


def build_esm2() -> None:
    d = pd.read_csv(COHORT, low_memory=False)
    w, _, _ = fit_stabilized_iptw(d, CANDIDATE_PS_VARS)
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.3))
    for a, label in [(1, "De-escalated/stopped"), (0, "Continued broad-spectrum")]:
        axes[0].hist(
            w.loc[w["A"] == a, "ps_den"],
            bins=np.linspace(0, 1, 31),
            histtype="step",
            linewidth=1.7,
            label=label,
        )
        axes[1].hist(
            w.loc[w["A"] == a, "SW_A"],
            bins=30,
            histtype="step",
            linewidth=1.7,
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
        ax.legend(frameon=False, fontsize=12)
        ax.tick_params(axis="both", labelsize=9.2)
    savefig(fig, "ESM_Fig2_mimic_ps_weights")


def main() -> None:
    required = [
        HARM / "harmonized_mortality_results.csv",
        HARM / "harmonized_secondary_outcomes.csv",
        HARM / "mimic_progressive_adjustment.csv",
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


if __name__ == "__main__":
    main()
