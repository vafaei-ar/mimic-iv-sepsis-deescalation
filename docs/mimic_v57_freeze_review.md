# MIMIC-IV v5.7 final corrected freeze review

This document records the final scientific/code-review status of the MIMIC-IV v5.7 analysis used in the manuscript. It replaces the earlier provisional freeze note that predated the vital-sign audit and corrected inference rerun.

## 1. Publication status

The MIMIC scientific definitions are frozen for publication. The final corrected inference uses the established 96-hour landmark target-trial emulation, the day-3 result-availability microbiology phenotype, the frozen exposure definition, and the final audited propensity-score implementation.

The original v5.7 run directory is `mimic_iv_v5_7_final_20260820T003506Z`. Corrected inference was generated after the vital-sign audit without changing the target trial, exposure, outcomes, microbiology eligibility, or weighting estimands.

Final analytic cohort:

- total: 9,589 admissions;
- de-escalated/stopped: 1,863;
- continued broad-spectrum: 7,726;
- 30-day deaths after the 96-hour landmark: 1,657.

## 2. Why a corrected inference rerun was required

The initial v5.7 review found three measurement problems in candidate vital-sign covariates:

1. baseline temperature aggregation could mix Fahrenheit and Celsius values before patient-level aggregation;
2. direct GCS extraction was non-informative in the saved cohort;
3. direct FiO2 extraction was nearly constant and did not have sufficiently trustworthy source/unit semantics for the frozen PS.

These were treated as measurement-audit problems, not as opportunities to tune the treatment effect. Temperature handling was corrected at the reading level, and the demonstrably invalid/non-informative direct GCS and FiO2 covariates were excluded from the final corrected PS implementation. This exclusion was identified after the original model specification, so it must be described transparently as an audit-driven correction rather than as a prespecified modeling choice.

No outcome-driven covariate tuning was performed.

## 3. Primary MIMIC estimand and mortality result

The primary estimand is the ATE estimated using stabilized inverse probability treatment weighting.

Final corrected weighted 30-day mortality:

- de-escalated/stopped risk: 0.182995;
- continued broad-spectrum risk: 0.174619;
- risk difference (de-escalated minus continued): +0.008376, or +0.84 percentage points;
- bootstrap 95% CI for the risk difference: -0.043903 to +0.056714;
- risk ratio: 1.047966;
- successful bootstrap replicates: 1,000.

The interval crosses the null widely. The MIMIC analysis therefore does not support a claim that day-3 de-escalation either lowers or increases 30-day mortality in the selected landmark population.

## 4. Progressive adjustment

The progressive sequence is retained because it demonstrates how strongly the crude/apparently protective association depends on clinical-improvement and treatment-intensity adjustment.

| Model | Mortality RD | Interpretation |
|---|---:|---|
| M1 demographics/comorbidity | -0.041432 | Strong apparent protective association with limited adjustment |
| M2 + baseline severity/labs | -0.031895 | Association attenuates |
| M3 + near-decision clinical status | -0.020128 | Further attenuation |
| M4 + trajectories/intensity | +0.008376 | Final fully adjusted point estimate crosses to near-null/slightly positive |

The M4 point estimate is identical to the primary corrected mortality point estimate. A later audit found that the separately generated progressive-M4 bootstrap CI differs slightly from the primary bootstrap CI despite the identical point estimate. The primary outcome file has explicit 1,000-replicate provenance, so manuscript reporting uses the primary mortality CI from `primary_secondary_outcomes.csv`. The progressive table is used to show attenuation of point estimates rather than to substitute a different M4 confidence interval.

Do not choose between those CIs based on favorability.

## 5. Primary weighting diagnostics

Final stabilized IPTW diagnostics:

- maximum post-weighting absolute SMD: 0.132610;
- worst-balanced covariate: `temp_48_72h`;
- treated/de-escalated ESS: 613.96 of 1,863;
- continued ESS: 7,119.58 of 7,726;
- maximum stabilized weight: 30.86.

These diagnostics indicate residual measured imbalance/limited effective sample size in the treated group. They must be reported as limitations. They do not invalidate the analysis, but they prevent describing the primary MIMIC weighting as having excellent balance.

## 6. Weighting sensitivities

These analyses were used as robustness checks, not as post hoc replacements for the primary ATE.

| Weighting analysis | Mortality RD | RR | Max post-SMD | Key interpretation |
|---|---:|---:|---:|---|
| Stabilized IPTW, primary ATE | +0.008376 | 1.04797 | 0.13261 | Primary estimand |
| Overlap weighting, ATO | -0.009341 | 0.93912 | 0.05374 | Better measured balance but different target population |
| IPTW truncated 1/99 | -0.010567 | 0.93955 | 0.17725 | Truncation changes weights and worsens measured balance |
| IPTW truncated 2.5/97.5 | -0.014043 | 0.91987 | 0.25067 | More aggressive truncation further worsens measured balance |

Bootstrap 95% mortality intervals for these sensitivity estimands all cross the null:

- overlap RD: -0.043568 to +0.022547; RR: 0.740636 to 1.175857;
- truncation 1/99 RD: -0.044087 to +0.023857; RR: 0.740636 to 1.181479;
- truncation 2.5/97.5 RD: -0.044980 to +0.022547; RR: 0.740636 to 1.175857.

Overlap weighting should not silently replace the primary result because it estimates the overlap-population estimand rather than the ATE.

## 7. Secondary outcomes

Final corrected primary IPTW estimates:

| Outcome | Estimate | 95% bootstrap CI |
|---|---:|---:|
| Hospital-free days | +1.10181 days | -0.35003 to +2.86773 |
| Antibiotic-free days | +1.75492 days | +0.32863 to +3.34471 |
| Normalized systemic antibiotic exposure | -0.117022 | -0.137390 to -0.070591 |
| Normalized broad-spectrum exposure | -0.168978 | -0.184525 to -0.125014 |
| Late recurrent/persistent antibiotic-course risk | RD -0.070916 | -0.103739 to -0.019983 |

The antibiotic-burden outcomes are mechanically downstream of the treatment strategy and should be interpreted primarily as treatment separation/burden outcomes, not as independent proof of clinical benefit. The late recurrent/persistent course outcome is exploratory because discharge changes how long inpatient treatment can be observed.

## 8. Microbiology and time semantics

The MIMIC primary culture rule uses the information available by the 72-hour treatment decision: qualifying clinical microbiology must have been sampled and no positive clinical culture result may be available by day 3. Eventual specimen-based negativity is a sensitivity analysis rather than the primary rule.

Treatment is classified during 72-96 hours after the first qualifying broad-spectrum exposure, and outcome follow-up starts at the 96-hour landmark. This is therefore a 96-hour landmark estimand, not an estimand for every patient who is alive at 72 hours.

## 9. Propensity-score implementation guardrails

The final PS implementation follows these frozen rules:

- in PS preparation, continuous covariates are filled with the analysis-sample mean, standardized using that mean and population SD, and clipped to +/-8;
- binary covariates are filled with 0;
- the PS-preparation function does not create new missingness indicators; any explicitly modeled missingness flags come from upstream feature construction and the frozen variable list;
- duplicate PS variables identified during v5.7 review are removed;
- the denominator model uses binomial GLM with the documented small-ridge regularized fallback only when needed;
- treatment probabilities are clipped to [0.001, 0.999];
- stabilized IPTW is the primary weighting scheme;
- bootstrap inference refits the PS within each resample.

Some upstream feature-construction functions use their own frozen source-layer imputation before PS preparation. The statement above describes the final PS-preparation behavior in `src/sepsis_deescalation/stats.py` and should not be replaced by a blanket claim of median imputation.

Do not add rank-pruning, new missingness indicators, new trimming thresholds, alternate ML propensity models, or additional covariates after seeing treatment effects unless explicitly labeled as a new post-publication analysis.

## 10. Final interpretation

The final MIMIC result is best summarized as follows: the apparent protective mortality association seen with limited adjustment disappears after accounting for near-decision clinical status, trajectories, and treatment/diagnostic intensity. Day-3 de-escalation is associated with lower antibiotic burden, but the fully adjusted mortality estimate is compatible with both modest benefit and modest harm.

This is the frozen MIMIC result that should be compared with the PSU modified external replication.
