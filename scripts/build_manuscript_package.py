#!/usr/bin/env python3
"""Build manuscript-facing text and figure data from frozen harmonized outputs.

No models are refit and no patient-level data are read. This script consumes only the
aggregate harmonized publication tables created by build_publication_integration.py.
"""
from __future__ import annotations

from pathlib import Path
import pandas as pd

SRC = Path("outputs/publication_integration/harmonized")
OUT = Path("outputs/publication_integration/manuscript_package")


def pct(x: float) -> float:
    return 100.0 * float(x)


def main() -> None:
    mortality = pd.read_csv(SRC / "harmonized_mortality_results.csv")
    secondary = pd.read_csv(SRC / "harmonized_secondary_outcomes.csv")
    progressive = pd.read_csv(SRC / "mimic_progressive_adjustment.csv")
    weights = pd.read_csv(SRC / "weighting_diagnostics.csv")

    OUT.mkdir(parents=True, exist_ok=True)

    # Figure 2 data: MIMIC progressive adjustment on RD scale, percentage points.
    fig2 = progressive[["model", "risk_difference", "rd_lower_95", "rd_upper_95", "max_post_smd"]].copy()
    for col in ["risk_difference", "rd_lower_95", "rd_upper_95"]:
        fig2[col + "_pp"] = 100.0 * fig2[col]
    fig2.to_csv(OUT / "figure2_mimic_progressive_adjustment.csv", index=False)

    # Figure 4 data: MIMIC vs PSU primary mortality estimates. These are shown side-by-side,
    # not pooled, because the PSU analysis is a modified external replication.
    fig4 = mortality.iloc[:2][[
        "dataset_analysis", "cohort_n", "deescalated_n", "continued_n",
        "mortality_rd", "rd_ci95_low", "rd_ci95_high", "mortality_rr",
        "rr_ci95_low", "rr_ci95_high", "note"
    ]].copy()
    for col in ["mortality_rd", "rd_ci95_low", "rd_ci95_high"]:
        fig4[col + "_pp"] = 100.0 * fig4[col]
    fig4.to_csv(OUT / "figure4_mimic_psu_mortality.csv", index=False)

    # Main publication table: primary mortality plus stewardship/recovery outcomes.
    main_table = secondary.copy()
    main_table.to_csv(OUT / "table_main_secondary_outcomes.csv", index=False)

    m = mortality.iloc[0]
    p = mortality.iloc[1]
    prog_m1 = progressive.iloc[0]
    prog_m4 = progressive.iloc[-1]
    mw = weights.loc[weights["dataset"] == "MIMIC-IV"].iloc[0]
    pw = weights.loc[weights["dataset"] == "PSU"].iloc[0]

    sec = {(r.dataset, r.outcome): r for r in secondary.itertuples(index=False)}
    m_af = sec[("MIMIC-IV", "Antibiotic-free days")]
    p_af = sec[("PSU", "Antibiotic-free days")]
    m_sys = sec[("MIMIC-IV", "Normalized systemic antibiotic exposure")]
    p_sys = sec[("PSU", "Normalized systemic antibiotic exposure")]
    m_broad = sec[("MIMIC-IV", "Normalized broad-spectrum exposure")]
    p_broad = sec[("PSU", "Normalized broad-spectrum exposure")]

    text = f"""# Manuscript-facing frozen results package\n\nThis package is generated only from frozen aggregate harmonized outputs. No patient-level data are read and no models are refit.\n\n## Results draft\n\nIn MIMIC-IV, 9,589 eligible admissions were included, of whom 1,863 were classified as de-escalated/stopped and 7,726 as continued broad-spectrum therapy. After stabilized inverse probability weighting with the fully adjusted propensity model, 30-day mortality was 18.3% in the de-escalated/stopped group and 17.5% in the continued group, corresponding to a risk difference of {pct(m.mortality_rd):+.2f} percentage points (95% CI {pct(m.rd_ci95_low):+.2f} to {pct(m.rd_ci95_high):+.2f}) and a risk ratio of {m.mortality_rr:.3f}. The apparent mortality association attenuated progressively with richer adjustment: the risk difference moved from {pct(prog_m1.risk_difference):+.2f} percentage points in the demographics/comorbidity model to {pct(prog_m4.risk_difference):+.2f} percentage points in the model including near-decision trajectories and treatment intensity.\n\nThe PSU modified external replication included 19,841 encounters, with 5,346 de-escalated and 14,495 continued. Its publication-locked primary estimate was a 30-day mortality risk difference of {pct(p.mortality_rd):+.2f} percentage points (95% CI {pct(p.rd_ci95_low):+.2f} to {pct(p.rd_ci95_high):+.2f}) and a risk ratio of {p.mortality_rr:.3f}. Because the PSU implementation differs materially in data model, timing, route, and culture-result semantics, these estimates should be presented as a modified external replication rather than pooled as if they were exchangeable estimates of a single common effect.\n\nAcross both datasets, de-escalation was consistently associated with lower antibiotic burden. Antibiotic-free days were higher by {m_af.estimate:.2f} days (95% CI {m_af.ci95_low:.2f} to {m_af.ci95_high:.2f}) in MIMIC-IV and {p_af.estimate:.2f} days (95% CI {p_af.ci95_low:.2f} to {p_af.ci95_high:.2f}) in PSU. Normalized systemic antibiotic exposure was lower by {m_sys.estimate:.3f} (95% CI {m_sys.ci95_low:.3f} to {m_sys.ci95_high:.3f}) in MIMIC-IV and {p_sys.estimate:.3f} (95% CI {p_sys.ci95_low:.3f} to {p_sys.ci95_high:.3f}) in PSU. Normalized broad-spectrum exposure was lower by {m_broad.estimate:.3f} (95% CI {m_broad.ci95_low:.3f} to {m_broad.ci95_high:.3f}) in MIMIC-IV and {p_broad.estimate:.3f} (95% CI {p_broad.ci95_low:.3f} to {p_broad.ci95_high:.3f}) in PSU.\n\n## Discussion draft\n\nThe central finding is not that antibiotic de-escalation improves survival. In MIMIC-IV, an apparent mortality advantage under limited adjustment progressively disappeared after accounting for clinical improvement, near-decision status, and treatment intensity, leaving no clear evidence of mortality benefit or harm in the fully adjusted analysis. This attenuation is consistent with strong confounding by indication and by recovery trajectory, because patients who are improving are both more likely to have antibiotics narrowed or stopped and more likely to survive.\n\nThe PSU replication produced a more favorable mortality estimate, but the two data sources should not be forced into numerical agreement. The PSU analysis is a modified external replication rather than an exact transport of the MIMIC estimand and measurement process. Differences in culture-result availability, medication source and route semantics, target-clock implementation, and the underlying clinical data model plausibly contribute to heterogeneity in the mortality estimate. The most reproducible cross-dataset signal is therefore stewardship-related: de-escalation was associated with materially lower systemic and broad-spectrum antibiotic exposure and more antibiotic-free days in both datasets.\n\nThese findings support a cautious clinical interpretation. Among patients who remain clinically stable enough to be eligible for a de-escalation decision after empiric broad-spectrum therapy, de-escalation appears to reduce subsequent antibiotic exposure without a clear mortality penalty in MIMIC-IV. The more favorable PSU mortality estimate strengthens the case that the strategy is not obviously harmful, but it should not be interpreted as proof of a mortality benefit because residual confounding, positivity, and cross-dataset measurement differences remain.\n\n## Weighting diagnostics to report\n\n- MIMIC-IV maximum post-weighting absolute SMD: {mw.max_post_smd:.3f}; de-escalated ESS: {mw.ess_deescalated:.0f}; maximum weight: {mw.max_weight:.2f}.\n- PSU maximum post-weighting absolute SMD: {pw.max_post_smd:.3f}; de-escalated ESS: {pw.ess_deescalated:.0f}; maximum weight: {pw.max_weight:.2f}.\n\nThe late recurrent/persistent antibiotic-course outcome remains exploratory because observation is affected by discharge timing.\n"""
    (OUT / "manuscript_results_discussion_draft.md").write_text(text)


if __name__ == "__main__":
    main()
