# Manuscript-ready integration of MIMIC-IV and PSU results

Status: analysis-complete draft for manuscript integration. This document does not change any frozen analysis definition or estimate.

## Recommended framing

The manuscript should present MIMIC-IV as the primary target-trial emulation and Penn State as a modified external replication. The PSU analysis should not be described as an exact replication because several data elements could not be harmonized exactly: PSU uses a hospital-clock implementation rather than validated ICU timing, systemic broad-spectrum prescribing as the primary exposure proxy rather than verified IV treatment, no validated day-3 culture-result-availability restriction, and date-level discharge/death and medication interval semantics.

The central result is therefore not that two databases produced the same mortality effect. The central result is that day-3 de-escalation consistently reduced subsequent antibiotic burden, whereas the adjusted mortality association differed across data environments: essentially null in MIMIC-IV and lower mortality in the modified PSU replication.

## Results text

### MIMIC-IV primary analysis

In MIMIC-IV, the primary stabilized inverse-probability-weighted analysis did not show evidence of a difference in 30-day mortality between patients who de-escalated systemic broad-spectrum therapy at day 3 and those who continued broad-spectrum therapy. Weighted mortality risks were 18.3% and 17.5%, respectively (risk difference [RD], +0.84 percentage points; risk ratio [RR], 1.05; 95% bootstrap CI for the RD, -4.39 to +5.67 percentage points). Prespecified alternative weighting approaches were also compatible with no mortality difference. Overlap weighting produced an RD of -0.93 percentage points (RR, 0.94; 95% bootstrap CI for the RD, -4.36 to +2.25), while 1st/99th and 2.5th/97.5th percentile weight truncation produced RDs of -1.06 and -1.40 percentage points, respectively, with confidence intervals spanning the null.

Progressive adjustment materially changed the MIMIC-IV mortality association. The RD shifted from -4.14 percentage points with demographics/comorbidity adjustment to -3.19 after adding baseline severity, -2.01 after near-decision clinical status, and +0.84 after the full trajectory/treatment-intensity adjustment. This pattern indicates substantial confounding by clinical improvement and treatment intensity around the de-escalation decision.

De-escalation was nevertheless associated with lower subsequent treatment burden. In the primary weighted analysis, hospital-free days increased by 1.10 days (95% CI, -0.35 to +2.87), antibiotic-free days increased by 1.75 days (95% CI, +0.33 to +3.34), normalized systemic antibiotic exposure decreased by 0.117 (95% CI, -0.137 to -0.071), and normalized broad-spectrum exposure decreased by 0.169 (95% CI, -0.185 to -0.125). The lower late recurrent/persistent antibiotic-course estimate is retained as exploratory because post-discharge observation differs across patients.

### Modified PSU external replication

The modified PSU external replication included 19,841 eligible encounters, of which 5,346 (26.9%) were classified as de-escalated and 14,495 (73.1%) as continued broad-spectrum therapy. Covariate balance after stabilized IPTW was strong (maximum absolute post-weighting standardized mean difference, 0.023), with effective sample sizes of approximately 4,513 de-escalated and 14,222 continued encounters and a maximum stabilized weight of 7.20.

In contrast to MIMIC-IV, the PSU analysis showed a lower adjusted 30-day mortality risk among de-escalated encounters. Weighted mortality was approximately 10.1% with de-escalation and 12.6% with continued broad-spectrum therapy (RD, -2.56 percentage points; RR, 0.80). The 1,000-replicate bootstrap 95% CI was -3.61 to -1.55 percentage points for the RD and 0.72 to 0.88 for the RR. All 1,000 bootstrap replicates completed successfully.

This association was stable across prespecified weighting sensitivities. Overlap weighting produced an RD of approximately -2.37 percentage points, while 1st/99th and 2.5th/97.5th percentile truncation produced RDs of approximately -2.74 to -2.75 percentage points; all corresponding bootstrap confidence intervals excluded the null.

The PSU secondary outcomes were directionally concordant with the MIMIC-IV treatment-burden findings. In the primary IPTW analysis, de-escalation was associated with approximately 1.45 additional hospital-free days and 3.16 additional antibiotic-free days, together with substantially lower normalized systemic and broad-spectrum antibiotic exposure. The exploratory late recurrent/persistent antibiotic-course outcome was also lower with de-escalation but remains vulnerable to differential observation after discharge.

### Prespecified PSU robustness analyses

The PSU mortality association was not materially changed by either of the two prespecified implementation sensitivities identified before outcome-effect inspection. Reclassifying day-3 exposure using MED_ADMIN rather than PRESCRIBING yielded an RD of -2.54 percentage points (95% bootstrap CI, -3.44 to -1.55) and an RR of 0.80 (95% CI, 0.73 to 0.87). Expanding eligibility using the lenient 96-hour landmark rule increased the cohort to 23,937 encounters and yielded an RD of -2.39 percentage points (95% CI, -3.30 to -1.48) and an RR of 0.80 (95% CI, 0.73 to 0.87). Each sensitivity completed all 1,000 bootstrap replicates without failure.

## Cross-dataset results table

| Analysis | Cohort | De-escalated | Continued | Mortality RD | Mortality RR | 95% CI for RD | Max post-weighting |SMD| |
|---|---:|---:|---:|---:|---:|---:|---:|
| MIMIC-IV primary stabilized IPTW | 9,589 | 1,863 | 7,726 | +0.84 pp | 1.05 | -4.39 to +5.67 pp | 0.133 |
| MIMIC-IV overlap weighting | 9,589 | 1,863 | 7,726 | -0.93 pp | 0.94 | -4.36 to +2.25 pp | 0.054 |
| PSU modified replication, primary IPTW | 19,841 | 5,346 | 14,495 | -2.56 pp | 0.80 | -3.61 to -1.55 pp | 0.023 |
| PSU MED_ADMIN exposure sensitivity | 19,841 | 5,347 | 14,494 | -2.54 pp | 0.80 | -3.44 to -1.55 pp | 0.023 |
| PSU lenient-landmark sensitivity | 23,937 | 7,009 | 16,928 | -2.39 pp | 0.80 | -3.30 to -1.48 pp | 0.021 |

## Discussion text

In this target-trial emulation, day-3 de-escalation was consistently associated with lower subsequent antibiotic exposure, but the adjusted mortality association differed across datasets. In MIMIC-IV, an apparent mortality advantage under less complete adjustment attenuated toward the null after incorporating near-decision physiology, clinical trajectories, and treatment intensity. This attenuation supports the interpretation that clinical improvement strongly influences the decision to de-escalate and can create substantial confounding in observational comparisons. By contrast, the modified PSU external replication showed an approximately 20% lower relative 30-day mortality risk with de-escalation after weighting, and the estimate was stable across alternative weighting, MED_ADMIN exposure classification, and landmark definitions.

The difference between MIMIC-IV and PSU should not be interpreted as evidence that one analysis is correct and the other is incorrect. The two implementations differ in clinically important measurement features. MIMIC-IV permits a stricter ICU-based clock and day-3 culture-result-availability restriction, whereas the available PSU PCORnet data required a hospital-clock implementation, did not permit faithful reconstruction of microbiology result availability by day 3, used systemic broad-spectrum prescribing as the primary treatment proxy, and relied on date-level timing for several constructs. These differences can change both the eligible population and residual confounding structure. The PSU mortality result therefore provides evidence of external robustness of the association within a modified implementation, rather than a direct replication of the MIMIC-IV causal contrast.

Across both datasets, the most reproducible finding was reduced antibiotic burden after de-escalation. This treatment separation is expected in part because antibiotic exposure is downstream of the treatment strategy itself, and it should not be interpreted as independent proof of clinical benefit. However, the preservation of similar directions for hospital-free and antibiotic-free days, together with the absence of a detectable mortality penalty in MIMIC-IV and the lower adjusted mortality association in PSU, supports the clinical plausibility of day-3 de-escalation among appropriately selected patients while retaining uncertainty about the causal mortality effect.

## Limitations to add or revise

The PSU analysis is a modified external replication rather than an exact transport of the MIMIC-IV target trial. The local PCORnet extract did not support validated ICU-entry timing or faithful day-3 culture-result-availability reconstruction. PRESCRIBING therefore served as the primary systemic broad-spectrum treatment proxy, with MED_ADMIN used as a sensitivity analysis; route information was insufficient to label the PSU treatment phenotype as verified intravenous therapy. Discharge, death, and several medication interval fields were available only at calendar-date resolution. These limitations can introduce exposure, eligibility, and outcome misclassification and may partly explain the difference in mortality estimates across datasets.

Residual confounding remains possible in both datasets. In MIMIC-IV, the substantial attenuation across progressive adjustment models demonstrates the importance of clinical trajectories around the de-escalation decision. Although measured covariate balance was excellent in PSU, the modified data representation omits several MIMIC-IV adjustment constructs, including validated ICU timing, microbiology intensity, direct GCS/FiO2 terms, ventilation procedure indicators, urine output, and exact pre-decision broad-spectrum antibiotic hours. The PSU mortality association should therefore be described as adjusted and robust to the prespecified analyses, not as proof that de-escalation causally reduces mortality.

Late recurrent/persistent antibiotic-course outcomes are exploratory in both datasets because follow-up observation differs after hospital discharge. Antibiotic-free days and normalized antibiotic exposure are also mechanically influenced by the treatment definition and should be interpreted primarily as treatment-burden measures.

## Manuscript lock note

Before copying exact MIMIC-IV confidence intervals into the final submitted manuscript, reconcile the known discrepancy between the primary corrected bootstrap CI and the separately generated progressive-model M4 bootstrap CI, which share the same corrected point estimate but differ slightly in the saved interval. This is a reporting audit only and should not trigger model retuning or additional outcome-driven analyses.
