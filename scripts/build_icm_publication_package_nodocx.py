#!/usr/bin/env python3
"""Build an Intensive Care Medicine submission package without DOCX dependencies.

Consumes only frozen aggregate publication outputs. No patient-level data are read, no models
are refit, and no estimands are redefined. DOCX assembly is intentionally deferred to the
separate document-production step after this aggregate package is validated.
"""
from __future__ import annotations

from pathlib import Path
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HARM = Path("outputs/publication_integration/harmonized")
OUT = Path("outputs/publication_integration/icm_package")


def pct(x: float) -> float:
    return 100.0 * float(x)


def wc(text: str) -> int:
    text = re.sub(r"[#*`|]", " ", text)
    return len(re.findall(r"\b[\w’'-]+\b", text))


def save_fig(fig, stem: str) -> None:
    fig.savefig(OUT / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)


def fig1() -> None:
    fig, ax = plt.subplots(figsize=(7.2, 2.1))
    ax.set_xlim(-4, 124); ax.set_ylim(-0.7, 1.0); ax.axis("off")
    ax.hlines(0, 0, 120, linewidth=1.5, color="#1f4e79")
    ax.annotate("", xy=(120, 0), xytext=(116, 0), arrowprops=dict(arrowstyle="->", lw=1.5, color="#1f4e79"))
    for x, label in [(0, "First qualifying\nbroad-spectrum antibiotic"), (72, "Day-3\ndecision"), (96, "Landmark /\nfollow-up start"), (120, "30-day outcome\nhorizon continues")]:
        ax.vlines(x, -0.12, 0.12, linewidth=1.5, color="#1f4e79")
        ax.text(x, 0.18, label, ha="center", va="bottom", fontsize=8)
    ax.plot([0,72],[-0.18,-0.18],lw=4,color="#9ecae1",solid_capstyle="butt")
    ax.plot([72,96],[-0.18,-0.18],lw=4,color="#4f81bd",solid_capstyle="butt")
    ax.plot([96,120],[-0.18,-0.18],lw=4,color="#1f4e79",solid_capstyle="butt")
    ax.text(36,-0.36,"Pre-decision covariates",ha="center",fontsize=8)
    ax.text(84,-0.36,"Exposure classification (72–96 h)",ha="center",fontsize=8)
    ax.text(108,-0.36,"Outcome follow-up",ha="center",fontsize=8)
    fig.tight_layout(); save_fig(fig, "Fig1_target_trial_timeline")


def fig2(progressive: pd.DataFrame) -> None:
    y=np.arange(len(progressive))[::-1]
    est=100*progressive["risk_difference"].to_numpy(float)
    lo=100*progressive["rd_lower_95"].to_numpy(float)
    hi=100*progressive["rd_upper_95"].to_numpy(float)
    fig,ax=plt.subplots(figsize=(7.2,3.4))
    ax.errorbar(est,y,xerr=np.vstack([est-lo,hi-est]),fmt="o",capsize=3,lw=1.3,color="#1f4e79",ecolor="#4f81bd")
    ax.axvline(0,lw=1,linestyle="--",color="#666666")
    ax.set_yticks(y); ax.set_yticklabels(["M1 demographics/comorbidity","M2 + baseline severity/labs","M3 + near-decision status","M4 + trajectories/intensity"],fontsize=8)
    ax.set_xlabel("30-day mortality risk difference, percentage points",fontsize=9)
    ax.grid(axis="x",linewidth=0.4,alpha=0.3)
    for yi,e in zip(y,est): ax.text(e+0.18,yi,f"{e:+.2f}",va="center",fontsize=8)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    fig.tight_layout(); save_fig(fig, "Fig2_progressive_adjustment")


def fig3(mortality: pd.DataFrame, secondary: pd.DataFrame) -> None:
    fig,axes=plt.subplots(1,2,figsize=(8.5,3.6))
    d=mortality.iloc[:2]
    est=100*d["mortality_rd"].to_numpy(float); lo=100*d["rd_ci95_low"].to_numpy(float); hi=100*d["rd_ci95_high"].to_numpy(float)
    axes[0].errorbar(est,[1,0],xerr=np.vstack([est-lo,hi-est]),fmt="o",capsize=3,lw=1.3,color="#1f4e79",ecolor="#4f81bd")
    axes[0].axvline(0,lw=1,linestyle="--",color="#666666")
    axes[0].set_yticks([1,0]); axes[0].set_yticklabels(["MIMIC-IV","Penn State"],fontsize=8)
    axes[0].set_xlabel("Mortality RD, percentage points",fontsize=8); axes[0].set_title("a  30-day mortality",loc="left",fontsize=9,fontweight="bold")
    axes[0].grid(axis="x",linewidth=0.4,alpha=0.3)
    ax=axes[1]; ax.axis("off"); ax.set_title("b  Stewardship outcomes",loc="left",fontsize=9,fontweight="bold")
    y0=0.88
    for outcome,label in [("Antibiotic-free days","Antibiotic-free days"),("Normalized systemic antibiotic exposure","Systemic exposure"),("Normalized broad-spectrum exposure","Broad-spectrum exposure")]:
        ax.text(0.02,y0,label,transform=ax.transAxes,fontsize=8.5,fontweight="bold",va="top"); y0-=0.10
        for _,r in secondary[secondary["outcome"]==outcome].iterrows():
            unit=" days" if r.unit=="days" else ""
            ax.text(0.05,y0,f"{r.dataset}: {r.estimate:+.2f} ({r.ci95_low:+.2f}, {r.ci95_high:+.2f}){unit}",transform=ax.transAxes,fontsize=8,va="top"); y0-=0.085
        y0-=0.04
    fig.tight_layout(); save_fig(fig, "Fig3_cross_dataset_results")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    mortality=pd.read_csv(HARM/"harmonized_mortality_results.csv")
    secondary=pd.read_csv(HARM/"harmonized_secondary_outcomes.csv")
    progressive=pd.read_csv(HARM/"mimic_progressive_adjustment.csv")
    weights=pd.read_csv(HARM/"weighting_diagnostics.csv")
    m,p=mortality.iloc[0],mortality.iloc[1]
    mw=weights[weights.dataset=="MIMIC-IV"].iloc[0]; pw=weights[weights.dataset=="PSU"].iloc[0]
    sec={(r.dataset,r.outcome):r for r in secondary.itertuples(index=False)}
    maf=sec[("MIMIC-IV","Antibiotic-free days")]; paf=sec[("PSU","Antibiotic-free days")]
    ms=sec[("MIMIC-IV","Normalized systemic antibiotic exposure")]; ps=sec[("PSU","Normalized systemic antibiotic exposure")]
    mb=sec[("MIMIC-IV","Normalized broad-spectrum exposure")]; pb=sec[("PSU","Normalized broad-spectrum exposure")]
    mh=sec[("MIMIC-IV","Hospital-free days")]; ph=sec[("PSU","Hospital-free days")]

    title="Day-3 antibiotic de-escalation in suspected ICU infection: a target trial emulation"
    abstract=(
        "**Purpose:** To estimate mortality and antibiotic-burden outcomes after a day-3 broad-spectrum antibiotic de-escalation strategy while accounting for pre-decision clinical improvement.\n\n"
        "**Methods:** We emulated a 72-hour treatment decision with a 96-hour landmark in MIMIC-IV and performed a modified external replication in Penn State clinical data. Stabilized inverse probability treatment weighting adjusted for demographics, comorbidity, physiology, laboratory values, clinical trajectories and treatment intensity measured before 72 hours. The primary outcome was 30-day mortality.\n\n"
        f"**Results:** MIMIC-IV included 9,589 admissions (1,863 de-escalated/stopped; 7,726 continued). Fully adjusted mortality was 18.3% versus 17.5% (risk difference {pct(m.mortality_rd):+.2f} percentage points, 95% CI {pct(m.rd_ci95_low):+.2f} to {pct(m.rd_ci95_high):+.2f}; risk ratio {m.mortality_rr:.3f}). The apparent mortality advantage attenuated from -4.14 to +0.84 percentage points with richer adjustment. In Penn State (n=19,841), the modified-replication mortality risk difference was {pct(p.mortality_rd):+.2f} percentage points (95% CI {pct(p.rd_ci95_low):+.2f} to {pct(p.rd_ci95_high):+.2f}). Antibiotic burden was lower in both datasets.\n\n"
        "**Conclusions:** In MIMIC-IV, day-3 de-escalation was not clearly associated with higher or lower 30-day mortality after clinical-trajectory and treatment-intensity adjustment. The modified Penn State replication yielded a lower-mortality association, but differing eligibility and data semantics preclude treating the estimates as exchangeable; lower antibiotic burden was the most consistent cross-dataset finding."
    )
    intro=(
        "Early broad-spectrum antibiotics are central to suspected sepsis care, but prolonged exposure increases treatment burden and antimicrobial selection pressure. Stewardship guidance therefore emphasizes reassessment and de-escalation when clinical and microbiological information becomes available [1,2]. Culture-negative suspected ICU infection is a particularly difficult setting because negative cultures reduce microbiological justification for broad therapy without proving that infection is absent.\n\n"
        "Observational comparisons of de-escalation are vulnerable to confounding because clinicians preferentially narrow or stop therapy in patients already improving. Prior studies generally have not identified increased mortality after de-escalation [3-6], and recent cohorts have emphasized the importance of early clinical trajectories [12,13]. We therefore emulated a day-3 treatment decision in MIMIC-IV with explicit pre-decision trajectory adjustment and then evaluated transportability in a modified Penn State external replication."
    )
    methods=(
        "### Study design and treatment strategies\n\nWe used MIMIC-IV 3.1 and a retrospective target-trial emulation with a 96-hour landmark. The trial clock began at first qualifying prescribed systemic intravenous broad-spectrum antibiotic exposure. MIMIC-IV eligibility required adult ICU care, early broad-spectrum therapy, clinical microbiology sampling, no positive clinical culture result available by 72 hours, no active vasopressor overlap in the prior 6 hours, broad coverage during 48-72 hours, and survival and hospitalization through 96 hours. The Penn State analysis preserved the conceptual 72-hour decision, 72-96-hour exposure window and 96-hour landmark, but was a modified external replication because faithful culture-result availability, exact ICU timing and verified intravenous route were unavailable.\n\n"
        "De-escalation/stopping was defined as no broad-spectrum overlap during 72-96 hours; continuation required any broad-spectrum overlap. Penn State used a systemic broad-spectrum PRESCRIBING proxy, with MED_ADMIN as a prespecified sensitivity.\n\n"
        "### Outcomes and statistical analysis\n\nThe primary outcome was 30-day all-cause mortality from the 96-hour landmark. Secondary outcomes were hospital-free days, antibiotic-free days, normalized systemic antibiotic exposure and normalized broad-spectrum exposure. We estimated stabilized inverse probability treatment weights for the average treatment effect using only pre-72-hour covariates. The corrected MIMIC model included demographics, comorbidity, physiology, laboratory values, near-decision status, improvement trajectories and treatment intensity. Percentile bootstrap confidence intervals used 1,000 resamples. A causal interpretation requires conditional exchangeability, positivity, consistency and correct measurement; clinician-perceived improvement, source control, imaging and judgement remain incompletely captured."
    )
    results=(
        f"The corrected MIMIC-IV cohort included 9,589 admissions: 1,863 de-escalated/stopped and 7,726 continued. The strict Penn State modified-replication cohort included 19,841 encounters: 5,346 de-escalated and 14,495 continued.\n\n"
        f"In MIMIC-IV, the mortality risk difference moved from -4.14 percentage points with demographics/comorbidity adjustment to -3.19 after baseline severity/laboratory values, -2.01 after near-decision clinical status and {pct(m.mortality_rd):+.2f} after trajectories and treatment intensity (Fig. 2). Fully adjusted weighted mortality was 18.3% versus 17.5%, corresponding to a risk difference of {pct(m.mortality_rd):+.2f} percentage points (95% CI {pct(m.rd_ci95_low):+.2f} to {pct(m.rd_ci95_high):+.2f}) and risk ratio {m.mortality_rr:.3f}. The confidence interval was compatible with benefit and harm and does not establish equivalence. Residual weighting limitations remained (maximum post-weighting absolute SMD {mw.max_post_smd:.3f}; treated ESS {mw.ess_deescalated:.0f}; maximum weight {mw.max_weight:.2f}).\n\n"
        f"In Penn State, the publication-locked primary mortality risk difference was {pct(p.mortality_rd):+.2f} percentage points (95% CI {pct(p.rd_ci95_low):+.2f} to {pct(p.rd_ci95_high):+.2f}) and risk ratio {p.mortality_rr:.3f} (95% CI {p.rr_ci95_low:.3f} to {p.rr_ci95_high:.3f}). Measured balance was stronger (maximum post-weighting absolute SMD {pw.max_post_smd:.3f}), but this does not eliminate unmeasured confounding.\n\n"
        f"Antibiotic-free days were higher by {maf.estimate:.2f} days (95% CI {maf.ci95_low:.2f} to {maf.ci95_high:.2f}) in MIMIC-IV and {paf.estimate:.2f} days (95% CI {paf.ci95_low:.2f} to {paf.ci95_high:.2f}) in Penn State. Normalized systemic antibiotic exposure was lower by {ms.estimate:.3f} (95% CI {ms.ci95_low:.3f} to {ms.ci95_high:.3f}) and {ps.estimate:.3f} (95% CI {ps.ci95_low:.3f} to {ps.ci95_high:.3f}), respectively; normalized broad-spectrum exposure was lower by {mb.estimate:.3f} (95% CI {mb.ci95_low:.3f} to {mb.ci95_high:.3f}) and {pb.estimate:.3f} (95% CI {pb.ci95_low:.3f} to {pb.ci95_high:.3f}). Hospital-free days differed by {mh.estimate:.2f} days in MIMIC-IV and {ph.estimate:.2f} days in Penn State."
    )
    discussion=(
        "The central finding is not that antibiotic de-escalation improves survival. In MIMIC-IV, an apparent mortality advantage under limited adjustment progressively disappeared after accounting for near-decision clinical status, improvement trajectories and treatment intensity. This attenuation is consistent with substantial confounding by indication and recovery trajectory. The final MIMIC estimate therefore supports no clear evidence of mortality benefit or harm after measured adjustment, not proof of equivalence or safety.\n\n"
        "The Penn State analysis yielded a more favorable mortality estimate, but the two datasets should not be forced into numerical agreement. Penn State lacked faithful day-3 culture-result availability, exact ICU timing and verified intravenous route, and several outcomes used date-level semantics. These differences alter both eligibility and measurement. Excellent measured balance in Penn State does not eliminate unmeasured confounding, so the persistent lower-mortality association remains an observational finding rather than evidence that de-escalation causally reduces mortality.\n\n"
        "The most reproducible cross-dataset signal was reduced antibiotic burden. These outcomes are mechanically downstream of treatment strategy and should not be interpreted as independent causal proof of clinical benefit; they quantify the stewardship consequence of day-3 de-escalation.\n\n"
        "Strengths include the explicit target-trial and landmark framework, separation of treatment classification from follow-up, restriction of confounder measurement to the pre-decision period, progressive trajectory adjustment, 1,000-replicate bootstrap inference and a frozen external replication. Limitations include prescription-based MIMIC exposure, residual imbalance and low treated ESS, possible microbiology misclassification, exclusion of patients who died or were discharged before 96 hours, and incomplete capture of source control, imaging, clinician judgement, goals of care and immune status.\n\n"
        "In selected patients eligible for a day-3 decision, de-escalation reduced subsequent antibiotic exposure without a clear mortality penalty in MIMIC-IV. The modified Penn State replication did not show excess mortality but cannot be interpreted as confirmation of a causal survival benefit."
    )
    refs=[
        "Evans L, Rhodes A, Alhazzani W, et al. Surviving sepsis campaign: international guidelines for management of sepsis and septic shock 2021. Intensive Care Med. 2021;47:1181-1247.",
        "Tabah A, Bassetti M, Kollef MH, et al. Antimicrobial de-escalation in critically ill patients: a position statement. Intensive Care Med. 2020;46:245-265.",
        "Guo Y, Gao W, Yang H, Ma C, Sui S. De-escalation of empiric antibiotics in severe sepsis or septic shock: a meta-analysis. Heart Lung. 2016;45:454-459.",
        "Roper S, Wingler MJB, Cretella DA. Antibiotic de-escalation in critically ill patients with negative clinical cultures. Pharmacy. 2023;11:104.",
        "Kim YC, Kim JH, Ahn JY, et al. Discontinuation of glycopeptides in culture-negative severe sepsis or septic shock. Antibiotics. 2020;9:250.",
        "Patanwala AE, Abu Sardaneh A, Alffenaar J-WC, et al. Antibiotic de-escalation practices in the intensive care unit. Ann Pharmacother. 2025;59:311-318.",
        "Johnson AEW, Bulgarelli L, Shen L, et al. MIMIC-IV, a freely accessible electronic health record dataset. Sci Data. 2023;10:1.",
        "Matthews AA, Danaei G, Islam N, Kurth T. Target trial emulation: applying principles of randomised trials to observational studies. BMJ. 2022;378:e071108.",
        "Cashin AG, Hansford HJ, Hernán MA, et al. Transparent Reporting of Observational Studies Emulating a Target Trial: the TARGET Statement. JAMA. 2025;334:1084-1093.",
        "Benchimol EI, Smeeth L, Guttmann A, et al. The RECORD Statement. PLoS Med. 2015;12:e1001885.",
        "Burrows P, Brown R-A, Samuelsen A, Bonavia AS. Association between in-hospital antibiotic use and long-term outcomes in critically ill patients. Antimicrob Steward Healthc Epidemiol. 2025.",
        "Ohnuma T, Fuller M, Balamurugan P, et al. Antibiotic de-escalation and 30-day mortality in suspected bacterial culture-negative sepsis. J Crit Care. 2026;95:155646.",
        "Makino J, Honda M, Sato F, et al. Clinical impact of antimicrobial de-escalation in critically ill patients. J Infect Chemother. 2026;32:102946."
    ]
    take_home="Day-3 de-escalation reduced antibiotic burden in two datasets. In MIMIC-IV, an apparent mortality advantage disappeared after clinical-trajectory adjustment; Penn State showed a lower-mortality association but was a modified, nonexchangeable replication."
    tweet="Day-3 ICU antibiotic de-escalation lowered antibiotic burden; MIMIC mortality benefit disappeared after trajectory adjustment."
    table1=("| Protocol element | MIMIC-IV | Penn State modified replication |\n|---|---|---|\n| Decision time | 72 h after first qualifying broad-spectrum exposure | Same conceptual clock |\n| Exposure window | 72–96 h | 72–96 h |\n| De-escalation | No broad-spectrum overlap | PRESCRIBING proxy; MED_ADMIN sensitivity |\n| Follow-up | Starts at 96 h | Starts at 96 h |\n| Culture criterion | No positive clinical culture result available by 72 h | Not faithfully implementable |\n| Primary outcome | 30-day mortality | 30-day mortality |")
    table2=(f"| Outcome | MIMIC-IV | Penn State |\n|---|---:|---:|\n| 30-day mortality RD | {pct(m.mortality_rd):+.2f} pp ({pct(m.rd_ci95_low):+.2f}, {pct(m.rd_ci95_high):+.2f}) | {pct(p.mortality_rd):+.2f} pp ({pct(p.rd_ci95_low):+.2f}, {pct(p.rd_ci95_high):+.2f}) |\n| 30-day mortality RR | {m.mortality_rr:.3f} | {p.mortality_rr:.3f} ({p.rr_ci95_low:.3f}, {p.rr_ci95_high:.3f}) |\n| Hospital-free days | {mh.estimate:+.2f} ({mh.ci95_low:+.2f}, {mh.ci95_high:+.2f}) | {ph.estimate:+.2f} ({ph.ci95_low:+.2f}, {ph.ci95_high:+.2f}) |\n| Antibiotic-free days | {maf.estimate:+.2f} ({maf.ci95_low:+.2f}, {maf.ci95_high:+.2f}) | {paf.estimate:+.2f} ({paf.ci95_low:+.2f}, {paf.ci95_high:+.2f}) |\n| Normalized systemic exposure | {ms.estimate:+.3f} ({ms.ci95_low:+.3f}, {ms.ci95_high:+.3f}) | {ps.estimate:+.3f} ({ps.ci95_low:+.3f}, {ps.ci95_high:+.3f}) |\n| Normalized broad-spectrum exposure | {mb.estimate:+.3f} ({mb.ci95_low:+.3f}, {mb.ci95_high:+.3f}) | {pb.estimate:+.3f} ({pb.ci95_low:+.3f}, {pb.ci95_high:+.3f}) |")

    manuscript=f"""# {title}\n\nAuthor names: Alireza Vafaei Sadr, [student collaborators], Anthony S. Bonavia, Vida Abedi, and [additional collaborators]\n\nAffiliations: Department of Public Health Sciences, College of Medicine, Pennsylvania State University, Hershey, Pennsylvania, USA; [additional affiliations]\n\nCorresponding author: Alireza Vafaei Sadr, [email], [telephone]\n\nRunning title: Day-3 antibiotic de-escalation\n\n## Abstract\n\n{abstract}\n\n**Keywords:** antimicrobial stewardship; sepsis; intensive care; target trial emulation; de-escalation; MIMIC-IV\n\n## Take-home message\n\n{take_home}\n\n## 140-character social-media summary\n\n{tweet}\n\n## Introduction\n\n{intro}\n\n## Methods\n\n{methods}\n\n### Table 1. Target trial emulation and modified external replication\n\n{table1}\n\n## Results\n\n{results}\n\n## Discussion\n\n{discussion}\n\n### Table 2. Primary and secondary outcomes\n\n{table2}\n\n## Figure legends\n\n**Fig. 1** Target-trial timing and landmark design. Covariates were measured before the 72-hour decision; treatment was classified during 72–96 hours; follow-up began at 96 hours.\n\n**Fig. 2** Progressive adjustment of the MIMIC-IV mortality association. Risk differences are de-escalation/stopping minus continued broad-spectrum therapy; error bars show 95% confidence intervals.\n\n**Fig. 3** Cross-dataset results. Panel a shows primary 30-day mortality risk differences. Panel b summarizes stewardship outcomes. Penn State is a modified external replication and estimates are not pooled.\n\n## Declarations\n\n**Ethics approval:** MIMIC-IV contains de-identified data. The Penn State replication used institutional clinical data under local governance. [Insert applicable Penn State IRB protocol number or determination.]\n\n**Funding:** [Insert final funding statement.]\n\n**Conflicts of interest:** [Insert final statement.]\n\n**Data availability:** MIMIC-IV is available to credentialed users through PhysioNet. Penn State source data are not publicly available because of institutional data-use restrictions.\n\n**Code availability:** Analysis code, phenotype definitions, aggregate tables and figure-generation code will be available in the public project repository. [Insert repository URL before submission.]\n\n## References\n\n""" + "\n".join(f"{i}. {r}" for i,r in enumerate(refs,1)) + "\n"

    (OUT/"ICM_manuscript_draft.md").write_text(manuscript)
    cover=("Dear Editors,\n\nWe submit the Original Paper entitled \""+title+"\" for consideration in Intensive Care Medicine. The study addresses a common ICU antimicrobial-stewardship decision using an explicit day-3 target-trial emulation in MIMIC-IV and a prespecified modified external replication in Penn State clinical data. The central result is that an apparent mortality advantage in MIMIC-IV attenuated to the null after adjustment for clinical recovery and treatment intensity, whereas antibiotic exposure was consistently lower in both datasets. The Penn State analysis yielded a lower-mortality association, which we present as informative heterogeneity rather than as confirmation of a causal survival benefit.\n\nThe manuscript is designed to meet ICM Original Paper limits, with three figures and two tables in the main text and detailed audits/robustness analyses assigned to electronic supplementary material.\n\nSincerely,\nAlireza Vafaei Sadr\n")
    (OUT/"ICM_cover_letter.md").write_text(cover)
    checklist=("# ICM submission checklist\n\n- Original Paper main text target: <=3,000 words.\n- Structured abstract: Purpose, Methods, Results, Conclusions.\n- Keywords: 6.\n- Main display items: 5 total (2 tables, 3 figures).\n- References: 13.\n- Two-sentence take-home message included.\n- <=140-character social-media summary included.\n- Numbered square-bracket citations used in main text.\n- Detailed cohort, sensitivity, microbiology, medication, balance and TARGET/RECORD material assigned to ESM.\n- Complete authors, affiliations, corresponding-author contact, funding, conflicts and Penn State IRB details before submission.\n- DOCX will be assembled and visually QA'd after this frozen aggregate package passes.\n")
    (OUT/"ICM_submission_checklist.md").write_text(checklist)
    main_text="\n\n".join([intro,methods,results,discussion])
    pd.DataFrame([{"section":"Abstract","words":wc(abstract)},{"section":"Introduction","words":wc(intro)},{"section":"Methods","words":wc(methods)},{"section":"Results","words":wc(results)},{"section":"Discussion","words":wc(discussion)},{"section":"Main text total","words":wc(main_text)}]).to_csv(OUT/"ICM_wordcounts.csv",index=False)
    fig1(); fig2(progressive); fig3(mortality,secondary)
    if wc(main_text)>3000: raise SystemExit(f"ICM main text exceeds 3000 words: {wc(main_text)}")
    if not 150 <= wc(abstract) <= 250: raise SystemExit(f"ICM abstract outside 150-250 words: {wc(abstract)}")
    if len(tweet)>140: raise SystemExit(f"ICM social-media summary exceeds 140 characters: {len(tweet)}")


if __name__=="__main__":
    main()
