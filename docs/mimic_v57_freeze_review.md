# MIMIC-IV v5.7 freeze review

This document records the review of the completed MIMIC-IV v5.7 package before manuscript lock and PSU external replication.

## Reviewed run

- Run directory: `mimic_iv_v5_7_final_20260820T003506Z`
- MIMIC-IV version: 3.1
- Cohort: 9,589 admissions
- De-escalated/stopped: 1,863
- Continued broad-spectrum: 7,726
- Deaths by 30-day horizon: 1,657
- Run manifest commit: `54dbf4a487f980c95666058d5cad986f79a26982`

The run package is internally consistent across cohort flow, primary outcomes, progressive adjustment, mortality sensitivities, weight diagnostics, and final weighting sensitivities.

## Primary estimand and main result

The primary estimand remains the ATE estimated with stabilized IPTW.

- Weighted 30-day mortality, de-escalated/stopped: 0.176359
- Weighted 30-day mortality, continued broad-spectrum: 0.175395
- Risk difference: +0.000964 (+0.10 percentage points)
- Risk ratio: 1.005499
- Bootstrap 95% CI for RD: -0.042897 to +0.043607

The mortality analysis therefore does not support a claim of mortality benefit or harm.

## Progressive adjustment

The apparent mortality advantage attenuates as clinical status, trajectories, and treatment intensity are added:

| Model | RD | RR | max post-SMD |
|---|---:|---:|---:|
| M1 demographics/comorbidity | -0.041432 | 0.772489 | 0.013280 |
| M2 + baseline severity/labs | -0.031895 | 0.823249 | 0.026551 |
| M3 + near-decision clinical status | -0.024092 | 0.865404 | 0.027100 |
| M4 + trajectories/intensity | +0.000964 | 1.005499 | 0.122087 |

This attenuation is an important substantive result and should be presented as evidence of strong confounding by clinical improvement and treatment intensity rather than as evidence that de-escalation reduces mortality.

## Stewardship and recovery outcomes

The weighted analyses show lower antibiotic exposure and more hospital-/antibiotic-free days among de-escalated/stopped patients:

- Hospital-free days: +1.319 days, 95% CI 0.141 to 2.444
- Antibiotic-free days: +1.955 days, 95% CI 0.469 to 3.145
- Normalized systemic antibiotic exposure: -0.1183, 95% CI -0.1388 to -0.0899
- Normalized broad-spectrum exposure: -0.1692, 95% CI -0.1863 to -0.1413

These outcomes are downstream of the treatment strategy and should be interpreted as treatment-burden/recovery measures, not as independent evidence of causal clinical benefit.

The late recurrent/persistent antibiotic-course outcome remains exploratory because observation is affected by discharge timing.

## Positivity and weighting sensitivity

Primary stabilized IPTW has residual positivity/balance limitations:

- max post-weighting SMD: 0.122087 (`temperature_48_72h`)
- treated ESS: 672.3 / 1,863
- continued ESS: 7,123.5 / 7,726
- maximum weight: 25.38

Sensitivity analyses:

| Analysis | RD | RR | max post-SMD |
|---|---:|---:|---:|
| Overlap weighting | -0.013405 | 0.914626 | 0.053776 |
| IPTW truncated 1/99 | -0.013174 | 0.924971 | 0.170405 |
| IPTW truncated 2.5/97.5 | -0.016301 | 0.907428 | 0.250367 |

Overlap weighting improves balance substantially but targets the overlap population rather than the ATE. Truncation worsens measured balance. Neither should replace the primary ATE post hoc.

## Mortality sensitivity analyses

- Narrowed/non-broad only vs continued: RD -0.025432, 95% CI -0.066774 to +0.018567
- Stopped all observed systemic antibiotics vs continued: RD +0.019560, 95% CI -0.044193 to +0.092321
- Excluding discharge within 24 h after landmark: RD -0.027874, 95% CI -0.127730 to +0.069928
- Strict test-name culture-only: RD +0.000839, 95% CI -0.034769 to +0.043418
- Eventual culture-negative: RD -0.021098, 95% CI -0.052926 to +0.013520

The eventual-culture-negative sensitivity reproduces the earlier approximately -2.1 percentage-point estimate, confirming that the old cohort definition was effectively an eventual-negativity phenotype.

## Freeze blocker found during final review

The package should **not yet be declared fully frozen** because the baseline temperature covariate `temp_max_pre72` mixes Fahrenheit and Celsius values before aggregation.

Evidence from the saved analytic cohort:

- `temp_max_pre72` range: 32.0 to 107.4
- median: 99.6
- 73 cohort rows have `temp_max_pre72 < 60`, indicating Celsius-scale maxima are present in a feature otherwise dominated by Fahrenheit-scale values
- `temperature_48_72h` is already normalized to Celsius and has a clinically plausible range (approximately 35.1 to 41.1 C)

This is a covariate-measurement problem in the propensity model, not an outcome/exposure-definition problem. A quick diagnostic transformation of values <60 from C to F moved the full-model mortality RD only modestly (about +0.10 pp to +0.22 pp), but the correct fix must operate on individual temperature readings before patient-level maxima are calculated because some stays may contain both Fahrenheit and Celsius observations.

Two additional vital features require verification before final lock:

- `gcs_total_48_72h` is constant at 0 in the saved cohort and is therefore dropped from the PS fit; this likely reflects extraction/mapping rather than true clinical measurement.
- `fio2_48_72h` is almost constant (median 1.0; only two distinct values in the saved cohort), so its source-item and unit semantics need confirmation.

A targeted audit script, `scripts/audit_mimic_vital_units.py`, was added to inspect the underlying MIMIC D_ITEMS/CHARTEVENTS labels and values without rerunning the full pipeline.

## Freeze decision

Current status: **provisionally locked for exposure, outcome, microbiology, landmark, and weighting definitions; not yet locked for the final covariate implementation.**

Before manuscript lock and PSU harmonization, complete the targeted vital-unit audit and either:

1. correct temperature/GCS/FiO2 extraction and rerun only the affected feature/inference stages; or
2. remove a demonstrably invalid/non-informative feature from the harmonized PS specification with a documented rationale, then rerun inference.

Do not tune any correction based on whether it makes the mortality estimate more favorable.
