#!/usr/bin/env python3
"""Build a Nature Portfolio / Nature Communications-compatible publication package.

This script consumes only frozen aggregate publication outputs already present in the project.
It does not read patient-level data, refit models, redefine the estimand, or change any frozen
scientific result. It produces manuscript text, submission/checklist material, and publication
figures from the frozen harmonized tables.
"""
from __future__ import annotations

from pathlib import Path
import textwrap

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HARM = Path("outputs/publication_integration/harmonized")
MAN = Path("outputs/publication_integration/manuscript_package")
OUT = Path("outputs/publication_integration/nature_package")


def pct(x: float) -> float:
    return 100.0 * float(x)


def save_fig(fig: plt.Figure, stem: str) -> None:
    fig.savefig(OUT / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)


def figure_timeline() -> None:
    fig, ax = plt.subplots(figsize=(7.2, 2.4))
    ax.set_xlim(-4, 124)
    ax.set_ylim(-0.8, 1.2)
    ax.axis("off")
    ax.hlines(0, 0, 120, linewidth=1.4)
    ax.annotate("", xy=(120, 0), xytext=(116, 0), arrowprops=dict(arrowstyle="->", lw=1.4))

    marks = [
        (0, "First qualifying\nbroad-spectrum antibiotic"),
        (72, "Day-3 decision"),
        (96, "Landmark /\nfollow-up start"),
        (120, "Outcome horizon\ncontinues to day 30"),
    ]
    for x, label in marks:
        ax.vlines(x, -0.12, 0.12, linewidth=1.4)
        ax.text(x, 0.23, label, ha="center", va="bottom", fontsize=8)

    ax.annotate("Pre-decision covariates", xy=(35, -0.30), ha="center", fontsize=8)
    ax.annotate("Exposure classification", xy=(84, -0.30), ha="center", fontsize=8)
    ax.annotate("72-96 h", xy=(84, -0.48), ha="center", fontsize=8)
    ax.annotate("Follow-up", xy=(108, -0.30), ha="center", fontsize=8)
    ax.plot([0, 72], [-0.17, -0.17], lw=3, solid_capstyle="butt")
    ax.plot([72, 96], [-0.17, -0.17], lw=3, solid_capstyle="butt")
    ax.plot([96, 120], [-0.17, -0.17], lw=3, solid_capstyle="butt")
    ax.text(0.01, 0.97, "a", transform=ax.transAxes, fontweight="bold", fontsize=9, va="top")
    ax.set_title("Target-trial timing and landmark design", fontsize=10, pad=8)
    save_fig(fig, "figure1_target_trial_timeline")


def figure_progressive(progressive: pd.DataFrame) -> None:
    d = progressive.copy()
    y = np.arange(len(d))[::-1]
    est = 100 * d["risk_difference"].to_numpy()
    lo = 100 * d["rd_lower_95"].to_numpy()
    hi = 100 * d["rd_upper_95"].to_numpy()
    xerr = np.vstack([est - lo, hi - est])

    fig, ax = plt.subplots(figsize=(7.2, 3.8))
    ax.errorbar(est, y, xerr=xerr, fmt="o", capsize=3, lw=1.2)
    ax.axvline(0, lw=1, linestyle="--")
    ax.set_yticks(y)
    ax.set_yticklabels([
        "M1 demographics/comorbidity",
        "M2 + baseline severity/labs",
        "M3 + near-decision clinical status",
        "M4 + trajectories/treatment intensity",
    ], fontsize=8)
    ax.set_xlabel("30-day mortality risk difference (percentage points)", fontsize=9)
    ax.set_title("Mortality association attenuates with richer pre-decision adjustment", fontsize=10)
    ax.grid(axis="x", linewidth=0.4, alpha=0.4)
    for yi, e in zip(y, est):
        ax.text(e + 0.18, yi, f"{e:+.2f}", va="center", fontsize=8)
    ax.text(0.01, 0.98, "a", transform=ax.transAxes, fontweight="bold", fontsize=9, va="top")
    fig.tight_layout()
    save_fig(fig, "figure2_mimic_progressive_adjustment")


def figure_mortality(mortality: pd.DataFrame) -> None:
    rows = mortality.iloc[:4].copy()
    labels = [
        "MIMIC-IV primary",
        "PSU primary modified replication",
        "PSU MED_ADMIN sensitivity",
        "PSU lenient landmark sensitivity",
    ]
    est = 100 * rows["mortality_rd"].to_numpy()
    lo = 100 * rows["rd_ci95_low"].to_numpy()
    hi = 100 * rows["rd_ci95_high"].to_numpy()
    y = np.arange(len(rows))[::-1]
    xerr = np.vstack([est - lo, hi - est])

    fig, ax = plt.subplots(figsize=(7.2, 4.0))
    ax.errorbar(est, y, xerr=xerr, fmt="o", capsize=3, lw=1.2)
    ax.axvline(0, lw=1, linestyle="--")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("30-day mortality risk difference (percentage points)", fontsize=9)
    ax.set_title("Mortality estimates differ across the primary and modified replication analyses", fontsize=10)
    ax.grid(axis="x", linewidth=0.4, alpha=0.4)
    for yi, e, l, h in zip(y, est, lo, hi):
        ax.text(h + 0.15, yi, f"{e:+.2f} ({l:+.2f}, {h:+.2f})", va="center", fontsize=7.5)
    ax.text(0.01, 0.98, "a", transform=ax.transAxes, fontweight="bold", fontsize=9, va="top")
    fig.tight_layout()
    save_fig(fig, "figure3_cross_dataset_mortality")


def figure_stewardship(secondary: pd.DataFrame) -> None:
    outcomes = [
        ("Antibiotic-free days", "days"),
        ("Normalized systemic antibiotic exposure", "proportion"),
        ("Normalized broad-spectrum exposure", "proportion"),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(9.0, 3.4))
    panel = ["a", "b", "c"]
    for ax, (outcome, unit), p in zip(axes, outcomes, panel):
        d = secondary.loc[secondary["outcome"] == outcome].copy()
        labels = d["dataset"].tolist()
        est = d["estimate"].to_numpy(dtype=float)
        lo = d["ci95_low"].to_numpy(dtype=float)
        hi = d["ci95_high"].to_numpy(dtype=float)
        y = np.arange(len(d))[::-1]
        xerr = np.vstack([est - lo, hi - est])
        ax.errorbar(est, y, xerr=xerr, fmt="o", capsize=3, lw=1.2)
        ax.axvline(0, lw=1, linestyle="--")
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=8)
        ax.grid(axis="x", linewidth=0.4, alpha=0.4)
        ax.set_title(outcome, fontsize=8.5)
        ax.set_xlabel(unit, fontsize=8)
        ax.text(0.02, 0.98, p, transform=ax.transAxes, fontweight="bold", fontsize=9, va="top")
    fig.suptitle("Antibiotic burden is lower after de-escalation in both datasets", fontsize=10)
    fig.tight_layout()
    save_fig(fig, "figure4_stewardship_outcomes")


def main() -> None:
    required = [
        HARM / "harmonized_mortality_results.csv",
        HARM / "harmonized_secondary_outcomes.csv",
        HARM / "mimic_progressive_adjustment.csv",
        HARM / "weighting_diagnostics.csv",
        MAN / "manuscript_results_discussion_draft.md",
    ]
    missing = [str(p) for p in required if not p.exists()]
    if missing:
        raise SystemExit("Missing frozen publication inputs: " + ", ".join(missing))

    OUT.mkdir(parents=True, exist_ok=True)
    mortality = pd.read_csv(required[0])
    secondary = pd.read_csv(required[1])
    progressive = pd.read_csv(required[2])
    weights = pd.read_csv(required[3])

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 8,
        "axes.spines.top": False,
        "axes.spines.right": False,
    })

    m = mortality.iloc[0]
    p = mortality.iloc[1]
    mw = weights.loc[weights["dataset"] == "MIMIC-IV"].iloc[0]
    pw = weights.loc[weights["dataset"] == "PSU"].iloc[0]
    sec = {(r.dataset, r.outcome): r for r in secondary.itertuples(index=False)}

    title = "Day three antibiotic de-escalation in suspected intensive care infection across two health systems"

    abstract = f"""Day-three reassessment of empiric broad-spectrum antibiotics is a common stewardship decision, but observational comparisons are confounded by clinical recovery. We emulated a 72-hour treatment decision in MIMIC-IV and performed a modified external replication in Penn State clinical data. In MIMIC-IV, 9,589 eligible admissions were classified as de-escalated or stopped (n=1,863) versus continued broad-spectrum therapy (n=7,726). After stabilized inverse-probability weighting with pre-decision clinical trajectories and treatment intensity, the 30-day mortality risk difference was {pct(m.mortality_rd):+.2f} percentage points (95% CI {pct(m.rd_ci95_low):+.2f} to {pct(m.rd_ci95_high):+.2f}). The apparent protective association attenuated progressively with richer adjustment. In the modified Penn State replication (n=19,841), the mortality risk difference was {pct(p.mortality_rd):+.2f} percentage points (95% CI {pct(p.rd_ci95_low):+.2f} to {pct(p.rd_ci95_high):+.2f}), but measurement and eligibility differences preclude treating the estimates as exchangeable. Across both datasets, de-escalation was consistently associated with more antibiotic-free days and lower systemic and broad-spectrum antibiotic exposure. These findings support antibiotic-burden reduction as the most reproducible cross-dataset signal, while mortality effects remain sensitive to adjustment and data-model differences."""

    intro = """## Introduction\n\nEarly broad-spectrum antibiotic therapy is central to the treatment of suspected sepsis, but continued exposure after the first several days creates avoidable selection pressure, toxicity and treatment burden. Current stewardship practice therefore emphasizes reassessment and de-escalation when clinical and microbiological information becomes available.1,2 Culture-negative suspected intensive care infection is a particularly difficult setting because negative cultures reduce microbiological justification for broad therapy without proving that infection is absent.\n\nEvidence on de-escalation is difficult to interpret causally because the treatment decision is strongly linked to recovery. Patients whose physiology, organ dysfunction and diagnostic findings improve before day three are more likely to have antibiotics narrowed or stopped and are also more likely to survive. Observational estimates that do not represent this pre-decision trajectory can therefore attribute the prognosis of recovery to the antibiotic decision itself. Prior studies generally have not identified increased mortality after de-escalation,3-6 and recent cohorts have emphasized the importance of early clinical trajectories.12,13\n\nWe therefore emulated a day-three treatment decision in MIMIC-IV using a 96-hour landmark and adjustment that explicitly incorporated pre-decision clinical trajectories and treatment intensity. We evaluated transportability in a modified external replication using Penn State clinical data. The replication preserved the decision and landmark logic where feasible while retaining documented differences in culture-result availability, medication representation and timing semantics. We focused on 30-day mortality as the primary safety outcome and separated this from antibiotic-burden outcomes that are expected to lie directly downstream of the treatment strategy.\n"""

    m_af = sec[("MIMIC-IV", "Antibiotic-free days")]
    p_af = sec[("PSU", "Antibiotic-free days")]
    m_sys = sec[("MIMIC-IV", "Normalized systemic antibiotic exposure")]
    p_sys = sec[("PSU", "Normalized systemic antibiotic exposure")]
    m_broad = sec[("MIMIC-IV", "Normalized broad-spectrum exposure")]
    p_broad = sec[("PSU", "Normalized broad-spectrum exposure")]
    m_hfd = sec[("MIMIC-IV", "Hospital-free days")]
    p_hfd = sec[("PSU", "Hospital-free days")]
    m_late = sec[("MIMIC-IV", "Late recurrent/persistent antibiotic course")]
    p_late = sec[("PSU", "Late recurrent/persistent antibiotic course")]

    results = f"""## Results\n\n### Cohorts and treatment groups\n\nThe corrected MIMIC-IV analytic cohort included 9,589 admissions, of which 1,863 (19.4%) were de-escalated or stopped during 72-96 hours and 7,726 (80.6%) continued broad-spectrum therapy. The strict Penn State modified-replication cohort included 19,841 encounters, with 5,346 (26.9%) de-escalated and 14,495 continued. The prespecified lenient-landmark Penn State sensitivity included 23,937 encounters (Fig. 1; Extended Data Table 1).\n\n### Progressive adjustment changes the MIMIC-IV mortality association\n\nThe MIMIC-IV mortality association attenuated as progressively richer pre-decision information was added. The 30-day mortality risk difference was -4.14 percentage points (95% CI -5.97 to -2.29) with demographics and comorbidity alone, -3.19 points (95% CI -5.01 to -1.39) after baseline severity and laboratory values, -2.01 points (95% CI -3.95 to -0.10) after near-decision clinical status, and +0.84 points after adding improvement trajectories and treatment intensity (Fig. 2). This pattern is consistent with substantial confounding by indication and recovery trajectory.\n\nIn the corrected fully adjusted MIMIC-IV model, weighted 30-day mortality was 18.3% under de-escalation/stopping and 17.5% under continued broad-spectrum therapy. The risk difference was {pct(m.mortality_rd):+.2f} percentage points (95% bootstrap CI {pct(m.rd_ci95_low):+.2f} to {pct(m.rd_ci95_high):+.2f}) and the risk ratio was {m.mortality_rr:.3f}. The confidence interval was compatible with both benefit and harm and does not establish equivalence. Residual balance and positivity limitations remained in the primary ATE weighting (maximum post-weighting absolute SMD {mw.max_post_smd:.3f}; treated effective sample size {mw.ess_deescalated:.0f}; maximum weight {mw.max_weight:.2f}).\n\n### The modified external replication yields a different mortality estimate\n\nIn Penn State, the publication-locked primary modified-replication estimate was a 30-day mortality risk difference of {pct(p.mortality_rd):+.2f} percentage points (95% CI {pct(p.rd_ci95_low):+.2f} to {pct(p.rd_ci95_high):+.2f}) and a risk ratio of {p.mortality_rr:.3f} (95% CI {p.rr_ci95_low:.3f} to {p.rr_ci95_high:.3f}). Prespecified administration-based exposure and lenient-landmark sensitivities were similar (Fig. 3). Measured balance was strong (maximum post-weighting absolute SMD {pw.max_post_smd:.3f}), but this does not remove unmeasured confounding. Because Penn State lacked a faithful day-three positive-culture availability rule, exact ICU timing and verified intravenous route, we treat these results as a modified external replication rather than an exchangeable estimate of the MIMIC-IV effect.\n\n### Antibiotic burden is consistently lower after de-escalation\n\nThe most consistent cross-dataset signal was lower antibiotic burden (Fig. 4). Antibiotic-free days were higher by {m_af.estimate:.2f} days (95% CI {m_af.ci95_low:.2f} to {m_af.ci95_high:.2f}) in MIMIC-IV and {p_af.estimate:.2f} days (95% CI {p_af.ci95_low:.2f} to {p_af.ci95_high:.2f}) in Penn State. Normalized systemic antibiotic exposure was lower by {m_sys.estimate:.3f} (95% CI {m_sys.ci95_low:.3f} to {m_sys.ci95_high:.3f}) and {p_sys.estimate:.3f} (95% CI {p_sys.ci95_low:.3f} to {p_sys.ci95_high:.3f}), respectively. Normalized broad-spectrum exposure was lower by {m_broad.estimate:.3f} (95% CI {m_broad.ci95_low:.3f} to {m_broad.ci95_high:.3f}) in MIMIC-IV and {p_broad.estimate:.3f} (95% CI {p_broad.ci95_low:.3f} to {p_broad.ci95_high:.3f}) in Penn State. Hospital-free days were higher by {m_hfd.estimate:.2f} days (95% CI {m_hfd.ci95_low:.2f} to {m_hfd.ci95_high:.2f}) in MIMIC-IV and {p_hfd.estimate:.2f} days (95% CI {p_hfd.ci95_low:.2f} to {p_hfd.ci95_high:.2f}) in Penn State. The late recurrent or persistent antibiotic-course outcome was lower in both datasets (MIMIC-IV RD {pct(m_late.estimate):+.2f} percentage points; Penn State RD {pct(p_late.estimate):+.2f} percentage points) but remains exploratory because post-discharge observation differs.\n"""

    discussion = """## Discussion\n\nThe central finding is not that antibiotic de-escalation improves survival. In MIMIC-IV, an apparent mortality advantage under limited adjustment progressively disappeared after accounting for near-decision clinical status, improvement trajectories and treatment intensity. The fully adjusted estimate provided no clear evidence of mortality benefit or harm. This attenuation is consistent with substantial confounding by indication and recovery trajectory because patients who are improving may be both more likely to have antibiotics narrowed or stopped and more likely to survive.\n\nThe Penn State replication produced a more favorable mortality estimate, but the two data sources should not be forced into numerical agreement. The external analysis is modified rather than exact. Local data did not support a faithful day-three positive-culture availability rule, exact ICU timing was unavailable, the primary antibiotic exposure was a broad-spectrum prescribing proxy with usually unspecified route, and several outcomes relied on calendar-date rather than exact timestamp semantics. These differences change the eligible population and measurement process. The Penn State estimate therefore represents a robust observational association within that data representation, not proof that de-escalation causally reduces mortality.\n\nThe most reproducible finding across datasets was treatment separation in antibiotic burden. De-escalation was associated with more antibiotic-free days and lower systemic and broad-spectrum exposure in both datasets. These outcomes are expected to be mechanically downstream of the treatment strategy and should not be interpreted as independent causal evidence of clinical benefit. Instead, they quantify the stewardship consequence of the day-three strategy while the mortality analysis addresses safety under the assumptions of the observational design.\n\nThis study has several strengths. The target-trial and landmark framework separates the 72-hour decision from the 96-hour follow-up start, restricts confounder measurement to the pre-decision period and directly models clinical-improvement trajectories. Inference used 1,000 bootstrap replicates, weighting diagnostics were prespecified, and the external analysis was conducted after covariate and outcome definitions were frozen. The Penn State analysis also included prespecified exposure and landmark robustness checks rather than effect-directed tuning.\n\nLimitations remain substantial. MIMIC-IV exposure is prescription-based and may not equal bedside administration; culture-negative status depends on available microbiology documentation; and the analysis excludes patients who died or were discharged before the 96-hour landmark. Residual imbalance and low treated effective sample size in MIMIC-IV leave both residual measured and unmeasured confounding possible. Penn State lacked faithful culture-result availability semantics, exact ICU timing and verified intravenous route; prescribing data represent ordered treatment and date-level outcome and antibiotic-day calculations are approximations. Neither dataset captures all determinants of de-escalation, including source control, imaging, clinician judgment, goals of care or immune status.\n\nTaken together, the cross-dataset evidence most consistently supports day-three reassessment as a stewardship strategy that reduces antibiotic exposure. Stronger conclusions about mortality require prospective evaluation or more closely harmonized multi-site observational designs with richer measurement of the clinical decision process.\n"""

    methods = """## Methods\n\n### Study design and target-trial framework\n\nWe conducted a retrospective cohort study using a target-trial emulation and a 96-hour landmark design. In MIMIC-IV, the trial clock began at first qualifying prescribed systemic intravenous broad-spectrum antibiotic exposure. Eligibility was evaluated at a 72-hour decision time, treatment strategies were classified during 72-96 hours, and outcome follow-up began at the 96-hour landmark. The estimand therefore applies to patients who remained alive and hospitalized through 96 hours. The Penn State analysis retained the same conceptual decision and landmark structure where feasible but was prespecified as a modified external replication because exact ICU timing, verified intravenous route and faithful culture-result availability were not available in the local data model.\n\n### Data sources and eligibility\n\nWe used MIMIC-IV version 3.1, a de-identified electronic health record database from Beth Israel Deaconess Medical Center.7 MIMIC-IV data sources included admissions, ICU stays, prescriptions, microbiology events, laboratory events, input/output events, diagnoses, procedures and death records. The modified external replication used Penn State clinical data organized in a PCORnet-like data model containing encounter, diagnosis, laboratory, prescribing, medication-administration and death information. Raw or row-level Penn State data remained within the institutional environment.\n\nMIMIC-IV eligibility required adult ICU admission; early prescribed systemic intravenous broad-spectrum antibiotic exposure; at least one qualifying non-screening clinical microbiology record from 24 hours before first exposure through the 72-hour decision time; no positive qualifying clinical culture result available by 72 hours; broad-spectrum coverage during 48-72 hours; no active vasopressor overlap during the 6 hours before the decision; and survival and hospitalization through 96 hours. The Penn State replication used the same decision and landmark concepts but did not impose the day-three culture-negative rule because independent result-availability semantics were not available.\n\n### Treatment strategies\n\nIn MIMIC-IV, de-escalation/stopping required no prescribed systemic intravenous broad-spectrum overlap during 72-96 hours, whereas continuation required any broad-spectrum overlap. The de-escalated group could retain non-broad systemic therapy or stop all observed systemic antibiotics. In Penn State, the primary exposure used PRESCRIBING records and a frozen systemic broad-spectrum antibiotic proxy because route was usually unspecified. MED_ADMIN exposure classification was prespecified as an administration-based sensitivity. Exposure therefore represents an operational treatment phenotype rather than clinician intent.\n\n### Covariates\n\nAll propensity-score covariates were measured before the 72-hour decision. The final MIMIC-IV model included demographics, comorbidities, baseline and near-decision physiology and laboratory values, organ-dysfunction and vital-sign trajectories, vasopressor trajectories, diagnostic and microbiology intensity, antibiotic intensity, steroids and a BMI proxy where available. Invalid or non-informative direct GCS and FiO2 features identified during measurement audit were excluded from the corrected primary propensity specification. Penn State covariates were frozen before outcome estimation and harmonized to available constructs, including demographics, comorbidity, selected laboratory and vital-sign trajectories, derived mean arterial pressure and vasopressor history. Missing continuous values were handled with the prespecified median-imputation strategy without adding missingness indicators in the primary model.\n\n### Outcomes\n\nThe primary safety outcome was all-cause mortality within 30 days after the 96-hour landmark. Hospital-free days were death-to-zero and otherwise represented days alive and outside the index hospitalization through 30 days. Antibiotic-free days were death-to-zero and otherwise represented days without observed systemic antibiotic coverage. Normalized systemic and broad-spectrum exposure divided observed antibiotic days by days alive through the 30-day horizon. Late recurrent or persistent antibiotic-course use from day 7 onward was exploratory because observation is affected by discharge timing. Penn State outcomes used date-level encounter and prescribing semantics and should be interpreted as approximations of the corresponding MIMIC-IV constructs.\n\n### Statistical analysis\n\nWe estimated stabilized inverse-probability treatment weights for the average treatment effect. Continuous propensity-score covariates were standardized; extreme linear predictors were clipped for numerical stability; and treatment probabilities used for stabilized weights were clipped to [0.001, 0.999]. Balance was evaluated using absolute standardized mean differences. Effective sample size and maximum weight were reported. Percentile bootstrap confidence intervals used 1,000 resamples. MIMIC-IV sensitivity analyses included overlap weighting and stabilized-weight truncation at the 1st/99th and 2.5th/97.5th percentiles. Penn State prespecified robustness analyses used MED_ADMIN exposure classification and a lenient 96-hour calendar-date landmark; bootstrap analyses refit the propensity model within each resample.\n\nA causal interpretation requires consistency, positivity, correct exposure and outcome classification, correct model specification, and conditional exchangeability given measured pre-decision covariates. Conditional exchangeability remains uncertain because clinician-perceived improvement, source control, imaging and clinical judgment are incompletely captured.\n\n### Reporting framework\n\nThe manuscript is organized to support STROBE, RECORD and TARGET reporting for an observational target-trial emulation using routinely collected health data.8-10 Detailed phenotype definitions, code lists, timing specifications, balance diagnostics and sensitivity analyses are assigned to Extended Data and Supplementary Information.\n"""

    references = """## References\n\n1. Evans, L. et al. Surviving sepsis campaign: international guidelines for management of sepsis and septic shock 2021. *Intensive Care Med.* **47**, 1181-1247 (2021).\n2. Tabah, A. et al. Antimicrobial de-escalation in critically ill patients: a position statement from a task force of the ESICM and ESCMID Critically Ill Patients Study Group. *Intensive Care Med.* **46**, 245-265 (2020).\n3. Guo, Y., Gao, W., Yang, H., Ma, C. & Sui, S. De-escalation of empiric antibiotics in patients with severe sepsis or septic shock: a meta-analysis. *Heart Lung* **45**, 454-459 (2016).\n4. Roper, S., Wingler, M. J. B. & Cretella, D. A. Antibiotic de-escalation in critically ill patients with negative clinical cultures. *Pharmacy* **11**, 104 (2023).\n5. Kim, Y. C. et al. Discontinuation of glycopeptides in patients with culture negative severe sepsis or septic shock: a propensity-matched retrospective cohort study. *Antibiotics* **9**, 250 (2020).\n6. Patanwala, A. E. et al. Antibiotic de-escalation practices in the intensive care unit: a multicenter observational study. *Ann. Pharmacother.* **59**, 311-318 (2025).\n7. Johnson, A. E. W. et al. MIMIC-IV, a freely accessible electronic health record dataset. *Sci. Data* **10**, 1 (2023).\n8. Matthews, A. A., Danaei, G., Islam, N. & Kurth, T. Target trial emulation: applying principles of randomised trials to observational studies. *BMJ* **378**, e071108 (2022).\n9. Cashin, A. G. et al. Transparent Reporting of Observational Studies Emulating a Target Trial: The TARGET Statement. *JAMA* **334**, 1084-1093 (2025).\n10. Benchimol, E. I. et al. The REporting of studies Conducted using Observational Routinely collected health Data (RECORD) Statement. *PLoS Med.* **12**, e1001885 (2015).\n11. Burrows, P., Brown, R.-A., Samuelsen, A. & Bonavia, A. S. Association between in-hospital antibiotic use and long-term outcomes in critically ill patients. *Antimicrob. Steward. Healthc. Epidemiol.* (2025). doi:10.1017/ash.2025.10054.\n12. Ohnuma, T. et al. Antibiotic de-escalation and 30-day mortality in patients with suspected bacterial culture-negative sepsis. *J. Crit. Care* **95**, 155646 (2026).\n13. Makino, J. et al. Clinical impact of antimicrobial de-escalation in critically ill patients: a single-center cohort study. *J. Infect. Chemother.* **32**, 102946 (2026).\n"""

    declarations = """## Data availability\n\nMIMIC-IV is available to credentialed users through PhysioNet under its data use agreement. The Penn State source data used for the modified external replication are not publicly available because of institutional data-use restrictions. Aggregate publication outputs and phenotype specifications can be shared through the project repository subject to institutional policy.\n\n## Code availability\n\nAnalysis and publication-generation code will be available in the public project repository: https://github.com/vafaei-ar/mimic-iv-sepsis-deescalation. No restricted patient-level data are included in the repository.\n\n## Acknowledgements\n\nThe authors acknowledge the developers and maintainers of MIMIC-IV and PhysioNet. [Add collaborator and institutional acknowledgements before submission.]\n\n## Funding\n\n[Insert complete funding statement before submission.]\n\n## Author contributions\n\n[Insert CRediT-style author contribution statement before submission.]\n\n## Competing interests\n\n[Insert complete competing-interests declaration before submission.]\n\n## Ethics\n\nMIMIC-IV contains de-identified data. This secondary analysis used de-identified MIMIC-IV data. The Penn State modified external replication used institutional clinical data under local governance. [Insert the applicable Penn State IRB protocol number or determination before submission.]\n\n## Correspondence\n\nCorrespondence and requests for materials should be addressed to Alireza Vafaei Sadr. [Insert submission email.]\n"""

    legends = """## Figure legends\n\n**Figure 1 | Target-trial timing and landmark design.** The trial clock begins at first qualifying broad-spectrum antibiotic exposure. Covariates are measured before the 72-hour decision, treatment is classified during 72-96 hours, and outcome follow-up begins at the 96-hour landmark.\n\n**Figure 2 | Mortality association attenuates with richer pre-decision adjustment in MIMIC-IV.** Points show weighted 30-day mortality risk differences for the progressive M1-M4 propensity specifications; bars show 95% bootstrap confidence intervals. The fully adjusted M4 model includes improvement trajectories and treatment intensity.\n\n**Figure 3 | Thirty-day mortality estimates in MIMIC-IV and the modified Penn State replication.** Points show risk differences for the corrected MIMIC-IV primary stabilized-IPTW analysis, the Penn State primary PRESCRIBING analysis, and prespecified Penn State MED_ADMIN and lenient-landmark sensitivities. Bars show 95% bootstrap confidence intervals. The estimates are displayed side-by-side and are not pooled because the external replication differs in eligibility and measurement semantics.\n\n**Figure 4 | Antibiotic burden after de-escalation in MIMIC-IV and Penn State.** Points show de-escalation-minus-continuation contrasts for antibiotic-free days, normalized systemic antibiotic exposure and normalized broad-spectrum antibiotic exposure. Bars show 95% bootstrap confidence intervals. These outcomes quantify treatment separation and should not be interpreted as independent proof of clinical benefit.\n"""

    manuscript = "\n".join([
        f"# {title}",
        "\n**Authors**: Alireza Vafaei Sadr, [student collaborators], Anthony S. Bonavia, Vida Abedi, [additional collaborators]",
        "\n**Affiliations**: Department of Public Health Sciences, College of Medicine, Pennsylvania State University, Hershey, Pennsylvania, USA; [additional affiliations]",
        "\n## Abstract\n\n" + abstract,
        intro,
        results,
        discussion,
        methods,
        references,
        legends,
        declarations,
    ])
    (OUT / "nature_manuscript.md").write_text(manuscript)

    cover = f"""# Draft cover letter for a Nature Portfolio journal\n\nDear Editors,\n\nPlease consider our manuscript, **\"{title}\"**, for publication as an Article.\n\nDay-three antibiotic reassessment is a common intensive-care stewardship decision, yet observational estimates are particularly vulnerable to confounding because clinicians preferentially de-escalate patients who are already recovering. We address this problem using a target-trial emulation in MIMIC-IV with explicit pre-decision clinical-trajectory and treatment-intensity adjustment, followed by a prespecified modified external replication in Penn State clinical data.\n\nThe principal result is methodological and clinical. In MIMIC-IV, an apparent protective mortality association attenuates from -4.14 percentage points under limited adjustment to {pct(m.mortality_rd):+.2f} percentage points after richer pre-decision adjustment, with a 95% confidence interval spanning both benefit and harm. The Penn State replication yields a lower-mortality association but differs materially in culture-result, timing and medication semantics, making transportability limitations visible rather than concealing them. Across both datasets, the most reproducible signal is reduced antibiotic burden after de-escalation.\n\nWe believe the work will interest readers in critical care, antimicrobial stewardship, causal inference and real-world evidence because it shows how clinical-trajectory adjustment can alter the apparent effect of a clinician-driven treatment decision and how the same target-trial framework behaves across distinct health-system data models. The analysis is reported with STROBE, RECORD and TARGET principles, and the project repository contains code and aggregate reproducibility outputs without restricted patient-level data.\n\nThis manuscript is not under consideration elsewhere. [Confirm before submission.] All authors have approved the manuscript and its submission. [Confirm before submission.]\n\nSincerely,\nAlireza Vafaei Sadr\n[Institution]\n[Email]\n"""
    (OUT / "nature_cover_letter.md").write_text(cover)

    checklist = f"""# Nature Portfolio publication assembly checklist\n\n## Format selected for this package\n- Working target: Nature Communications-compatible Article structure within broader Nature Portfolio style.\n- Title: {len(title.split())} words; Nature Communications guidance is <=15 words.\n- Abstract: {len(abstract.split())} words; guidance is <=200 words and no references.\n- Main structure: Introduction, Results, Discussion, Methods, Data Availability, Code Availability, References, Acknowledgements, Author Contributions, Competing Interests.\n- Figures use sans-serif lettering and simple layouts; panel labels are bold upright lowercase letters.\n\n## Scientific freeze\n- [x] MIMIC-IV corrected primary result preserved.\n- [x] PSU publication-locked primary result preserved.\n- [x] No pooling of MIMIC-IV and PSU mortality estimates.\n- [x] PSU described as a modified external replication.\n- [x] Late recurrent/persistent antibiotic-course outcome remains exploratory.\n- [x] No model refitting or new patient-level analysis in publication assembly.\n\n## Reporting standards\n- [x] Target-trial decision time, exposure window and landmark stated.\n- [x] Pre-decision confounder timing stated.\n- [x] Weighting diagnostics stated.\n- [x] Residual confounding and positivity limitations stated.\n- [x] Data-model differences in external replication stated.\n- [ ] Complete STROBE checklist for submission supplement.\n- [ ] Complete RECORD checklist for submission supplement.\n- [ ] Complete TARGET checklist for submission supplement.\n\n## Required author/institution items before submission\n- [ ] Final author list and order.\n- [ ] Complete affiliations and present addresses.\n- [ ] Corresponding-author email.\n- [ ] Penn State IRB protocol number or formal determination.\n- [ ] Funding statement.\n- [ ] CRediT author-contribution statement.\n- [ ] Competing-interest declaration.\n- [ ] Confirm all authors approve submission.\n- [ ] Confirm manuscript is not under consideration elsewhere.\n\n## Data and code\n- [x] MIMIC-IV credentialed-access statement.\n- [x] Penn State restricted-data statement.\n- [x] Public repository URL included for code.\n- [ ] Confirm exact release/tag/commit to cite at submission.\n- [ ] Add repository DOI if archived through Zenodo or equivalent.\n\n## Display items\n- Main Fig. 1: target-trial timing and landmark design.\n- Main Fig. 2: MIMIC-IV progressive adjustment.\n- Main Fig. 3: cross-dataset mortality estimates and Penn State sensitivities.\n- Main Fig. 4: stewardship outcomes in both datasets.\n- Main Table 1: concise cohort and primary/secondary outcomes.\n- Extended Data: cohort flow, balance/overlap diagnostics, phenotype definitions, sensitivity analyses.\n\n## Nature-style editorial checks\n- [x] Main message stated early and directly.\n- [x] Mortality not framed as a proven benefit.\n- [x] Technical audit detail moved out of the main narrative.\n- [x] Figure legends define the estimands and confidence intervals.\n- [ ] Final reference formatting and DOI verification.\n- [ ] Final line-numbered Word or PDF submission file.\n- [ ] Journal-specific reporting and editorial-policy forms.\n"""
    (OUT / "nature_submission_checklist.md").write_text(checklist)

    extended = """# Extended Data and Supplementary Information plan\n\n## Extended Data display items\n1. **Extended Data Fig. 1: cohort construction and eligibility flow.** Include MIMIC-IV and Penn State flow counts, with explicit separation of the modified external-replication eligibility rules.\n2. **Extended Data Fig. 2: propensity-score overlap and weight distributions.** Show MIMIC-IV primary stabilized IPTW and Penn State primary weighting.\n3. **Extended Data Fig. 3: Love plot / covariate balance.** Report pre- and post-weighting absolute SMDs; call out the MIMIC-IV maximum post-weighting SMD of 0.133.\n4. **Extended Data Table 1: target-trial protocol and cross-dataset harmonization.** Eligibility, decision time, treatment classification, landmark, outcome horizon and estimand.\n5. **Extended Data Table 2: progressive MIMIC-IV adjustment models.** M1-M4 covariate domains, estimates, confidence intervals and balance.\n6. **Extended Data Table 3: mortality sensitivity analyses.** MIMIC-IV overlap/truncation plus Penn State MED_ADMIN and lenient-landmark analyses.\n7. **Extended Data Table 4: antibiotic and route definitions.** Broad-spectrum and systemic definitions, exclusions and route rules.\n8. **Extended Data Table 5: microbiology phenotype and culture-result availability rules.** Primary and strict sensitivity definitions.\n9. **Extended Data Table 6: clinical-improvement covariate timing and missingness.** Include measurement windows and imputation.\n10. **Extended Data Table 7: reporting checklist mapping.** STROBE, RECORD and TARGET items and manuscript locations.\n\n## Supplementary Information\n- Detailed code lists and medication mappings.\n- Microbiology term audit and strict culture sensitivity.\n- Additional MIMIC-IV weighting sensitivities and diagnostic tables.\n- Full reporting checklists.\n- Reproducibility manifest identifying frozen publication source files and exact project commit.\n"""
    (OUT / "nature_extended_data_plan.md").write_text(extended)

    wordcounts = [
        ("title_words", len(title.split())),
        ("abstract_words", len(abstract.split())),
        ("introduction_words", len(intro.split())),
        ("results_words", len(results.split())),
        ("discussion_words", len(discussion.split())),
        ("methods_words", len(methods.split())),
        ("main_text_words_intro_results_discussion", len((intro + results + discussion).split())),
    ]
    pd.DataFrame(wordcounts, columns=["section", "words"]).to_csv(OUT / "nature_wordcounts.csv", index=False)

    figure_timeline()
    figure_progressive(progressive)
    figure_mortality(mortality)
    figure_stewardship(secondary)


if __name__ == "__main__":
    main()
