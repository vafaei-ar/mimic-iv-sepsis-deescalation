# MIMIC-IV v5.7 publication freeze review

This document records the accepted MIMIC-IV publication analysis after the final vital-sign audit and correction. Historical pre-correction values are not manuscript-facing results.

## Accepted publication source

Base source-dependent run:

- `outputs/mimic/mimic_iv_v5_7_final_20260820T003506Z`

Accepted corrected inference rerun:

- `outputs/mimic/mimic_iv_v5_7_final_20260820T003506Z/inference_reruns/final_vital_corrected_final_20260825T010041Z`

The corrected rerun is the source for manuscript primary/secondary outcomes, progressive adjustment, and final weighting diagnostics. Source-dependent microbiology and missing-stop-time sensitivities still originate from the complete base run.

## Final cohort and estimand

- Cohort: 9,589 admissions
- De-escalated or stopped: 1,863
- Continued broad-spectrum: 7,726
- Primary estimand: stabilized-IPTW ATE
- Treatment decision: 72 h after first qualifying broad-spectrum exposure
- Treatment classification: 72-96 h
- Follow-up start: 96-h landmark
- Primary culture-negative phenotype: qualifying microbiology sampled with no positive clinical culture result available by 72 h

The eventual culture-negative phenotype is retained only as a sensitivity analysis because using eventual results for day-3 eligibility would use future information.

## Final corrected primary mortality result

The publication-facing corrected primary result is:

- Weighted mortality, de-escalated/stopped: 0.1829947002
- Weighted mortality, continued broad-spectrum: 0.1746188734
- Risk difference: +0.0083758268 (+0.84 percentage points)
- Risk ratio: 1.0479663318
- 95% bootstrap CI for RD: -0.0439033352 to +0.0567136492 (-4.39 to +5.67 percentage points)

Interpretation: the fully adjusted MIMIC analysis did not show clear evidence of increased 30-day mortality, but clinically meaningful harm could not be excluded. These data do not establish noninferiority or safety.

## Progressive adjustment

The M1-M4 sequence demonstrates attenuation as clinical improvement and treatment intensity are incorporated:

| Model | RD, percentage points | RR | max post-SMD |
|---|---:|---:|---:|
| M1 demographics/comorbidity | -4.14 | 0.7725 | 0.0133 |
| M2 + baseline severity | -3.19 | 0.8232 | 0.0266 |
| M3 + near-decision clinical status | -2.01 | 0.8872 | 0.0272 |
| M4 + trajectories/intensity | +0.84 | 1.0480 | 0.1326 |

The progressive-model M4 point estimate matches the primary analysis exactly. Its separately generated bootstrap interval is not the manuscript primary interval.

### Confidence-interval reporting rule

For the manuscript primary mortality result, report the 1,000-replicate interval from `primary_secondary_outcomes.csv`: **-4.39 to +5.67 percentage points**. When presenting the M1-M4 sequence, the M4 point estimate may be shown with this designated primary interval only after asserting point-estimate equality.

A separately generated progressive-model M4 interval is a different saved bootstrap realization and should not replace the designated primary interval. Do not choose between those CIs based on favorability.

## Final corrected secondary outcomes

- Hospital-free days: +1.1018 days, 95% CI -0.3500 to +2.8677
- Antibiotic-free days: +1.7549 days, 95% CI +0.3286 to +3.3447
- Normalized systemic antibiotic exposure: -0.1170, 95% CI -0.13739 to -0.07059
- Normalized broad-spectrum exposure: -0.16898, 95% CI -0.18452 to -0.12501
- Late recurrent/persistent antibiotic-course RD: -0.07092, 95% CI -0.10374 to -0.01998

The late-course outcome is exploratory because discharge changes observation opportunity. Antibiotic-burden outcomes are mechanically downstream of treatment classification and are treatment-separation/burden measures rather than independent proof of clinical benefit.

## Final weighting diagnostics

Primary stabilized IPTW:

- maximum post-weighting absolute SMD: 0.13261
- worst-balanced covariate: `temperature_48_72h`
- de-escalated/stopped ESS: 613.96
- continued ESS: 7,119.58
- maximum stabilized weight: 30.864

Prespecified weighting sensitivities:

| Analysis | RD, percentage points | RR | 95% CI for RD, percentage points | max post-SMD |
|---|---:|---:|---:|---:|
| Overlap weighting | -0.93 | 0.939 | -4.36 to +2.25 | 0.0537 |
| IPTW truncated 1/99 | -1.06 | 0.940 | -4.41 to +2.39 | 0.177 |
| IPTW truncated 2.5/97.5 | -1.40 | 0.920 | -4.50 to +2.25 | 0.251 |

Overlap weighting targets the overlap population rather than the primary ATE and must not replace the primary estimand post hoc.

## Vital-sign audit and final correction

The final audit identified three measurement issues:

- baseline temperature mixed Fahrenheit and Celsius values before aggregation;
- the direct GCS-total feature was not a trustworthy routine measurement;
- direct FiO2 mapping was unreliable/sparse for propensity-score inclusion.

The accepted correction:

1. converts temperature at the individual-reading level before aggregation;
2. reconstructs GCS only from same-timestamp eye + verbal + motor component triplets;
3. normalizes audited FiO2 readings to fractions;
4. recomputes the SOFA-like trajectory after the GCS correction;
5. excludes direct GCS and FiO2 terms from the primary propensity score while retaining corrected values for derived scores/descriptives;
6. reruns inference from the corrected cohort.

The primary PS preparation retains the frozen historical MIMIC rules: binary missing values are filled with 0; continuous missing values use mean imputation in PS preparation; continuous variables are standardized and clipped to [-8, +8]; and estimated treatment probabilities are clipped to [0.001, 0.999].

These changes were audit-driven measurement repairs, not outcome-driven model tuning.

## Exploratory treatment subtypes

Among 1,863 de-escalated/stopped admissions:

- narrowed/non-broad only: 1,077 (57.8%)
- stopped all observed systemic antibiotics: 786 (42.2%)

Exploratory adjusted mortality sensitivities:

- narrowing vs continued: RD -1.97 percentage points, 95% CI -6.86 to +3.33
- stopping vs continued: RD +2.55 percentage points, 95% CI -5.14 to +9.97

The complete-stopping contrast has poor overlap (treated ESS about 106, max post-SMD about 0.303, maximum weight about 47.9) and remains supplementary only. Neither subtype analysis replaces the frozen binary primary estimand.

## Publication lock

Status: **frozen for manuscript use**.

Do not alter the cohort definition, treatment window, microbiology rule, propensity-score specification, weighting estimand, trimming thresholds, outcome definitions, or bootstrap strategy in response to treatment-effect results. Any new scientific specification should be explicitly labeled as a new post hoc analysis rather than folded into the publication freeze.
