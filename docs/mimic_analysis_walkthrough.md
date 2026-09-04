# MIMIC-IV publication analysis: reviewer walkthrough

This document is the code-review map for the MIMIC-IV analysis used in the manuscript. It explains both **what to run** and **why several non-obvious implementation decisions exist**.

## 1. Scientific clock and estimand

The MIMIC analysis is a 96-hour landmark target-trial emulation.

- `t0`: first qualifying systemic IV broad-spectrum antibiotic exposure from ICU `intime - 6 h` through `intime + 24 h`.
- treatment decision: `t0 + 72 h`.
- treatment classification: observed broad-spectrum overlap during `72-96 h` after `t0`.
- `A=1`: no qualifying broad-spectrum overlap during 72-96 h; this includes narrowing to non-broad therapy and stopping all observed systemic antibiotics.
- `A=0`: any qualifying broad-spectrum overlap during 72-96 h.
- follow-up begins at `t0 + 96 h`.
- primary outcome: 30-day all-cause mortality from the 96-hour landmark.

Patients who die or leave the hospital through the 96-hour classification window are not part of this landmark estimand. This avoids assigning a treatment strategy before the full classification window has been observed.

## 2. Primary microbiology rule

The primary MIMIC phenotype is **result-availability aware**.

A qualifying clinical microbiology specimen must be collected by the 72-hour decision time, and no qualifying positive clinical culture result may be available by that decision time.

This is intentionally different from asking whether the specimen eventually remains negative. The eventual specimen-based culture-negative phenotype is retained only as a sensitivity analysis because using future results to define day-3 eligibility would use information unavailable to the clinician at the treatment decision.

The relevant code is in:

- `src/sepsis_deescalation/microbiology.py`
- `src/sepsis_deescalation/mimic_pipeline.py`

## 3. Medication phenotype

The MIMIC primary exposure uses structured prescription/order records and inferred coverage intervals. It is therefore an operational medication-coverage phenotype, not proof of bedside administration and not a measure of clinician intention.

The broad-spectrum and systemic-antibiotic definitions are implemented in:

- `src/sepsis_deescalation/antibiotics.py`

The missing-stop-time fill rule is frozen in configuration and is separately examined with deterministic stop-time sensitivity analyses.

## 4. Why clinical-improvement trajectories are in the PS

Confounding by clinical improvement is the central causal problem. Patients who are recovering are more likely to have therapy narrowed or stopped and are also more likely to survive.

The final propensity model therefore includes information measured **before** the 72-hour treatment decision, including baseline illness, near-decision physiology, changes from early to late windows, vasopressor trajectories, diagnostic intensity, antibiotic intensity, steroid exposure, and BMI proxy where validly measured.

The M1-M4 sequence is retained to show how the mortality association changes as those clinical-improvement and treatment-intensity variables are added. It is not a model-selection exercise.

## 5. Final v5.7 vital-sign correction

The original v5.7 source-dependent run exposed a measurement issue during the final audit:

- baseline temperature mixed Fahrenheit and Celsius before aggregation;
- the direct GCS-total feature was not a trustworthy routine measurement;
- the direct FiO2 feature had unreliable/sparse semantics for direct PS inclusion.

The correction was handled as an **audit-driven measurement repair**, not as outcome-driven model tuning.

The accepted publication workflow therefore has two stages:

1. complete the source-dependent v5.7 run, which establishes cohort membership, exposure, microbiology, outcomes, and source-dependent sensitivities;
2. run `scripts/repair_v57_vital_covariates.py` on that completed analytic cohort and then rerun statistical inference from the corrected cohort.

The repair script:

- converts temperature at the individual-reading level before calculating the baseline maximum;
- reconstructs GCS only from same-timestamp eye + verbal + motor component triplets;
- normalizes audited FiO2 readings to fractions;
- recomputes the SOFA-like score after the GCS correction;
- refits the propensity score after replacing the audited covariates.

Direct GCS and FiO2 are excluded from the primary PS in `src/sepsis_deescalation/specification.py`. Reconstructed GCS still affects the SOFA-like trajectory, which is why the measurement repair matters even though GCS is not a direct PS term.

This two-stage provenance is kept visible rather than silently rewriting the historical base run.

## 6. Propensity-score preparation

The final MIMIC PS implementation is in `src/sepsis_deescalation/stats.py` and the frozen variable list is in `src/sepsis_deescalation/specification.py`.

Important implementation choices:

- binary 0/1 variables: missing values are filled with 0;
- continuous variables: missing values are filled with the variable mean in PS preparation, then standardized by the sample mean/SD and clipped to `[-8, +8]` for numerical stability;
- no new missingness indicators are created by the PS-preparation function;
- treatment probabilities are clipped to `[0.001, 0.999]` to prevent numerical explosion at machine-level probabilities;
- ordinary binomial GLM is used first;
- a very small ridge-regularized binomial fallback (`alpha=0.001`, L2) is used only if the ordinary fit fails or does not converge;
- stabilized IPTW estimates the primary ATE.

The mean-imputation rule here is a historical/frozen MIMIC implementation detail and is intentionally not changed to match PSU's source-layer median imputation after treatment effects were known.

## 7. Duplicate and invalid direct PS terms

Two exact duplicate PS variables found during the v5.7 audit are excluded so they do not create redundant columns:

- `clinical_micro_records_pre72` duplicates `micro_records_pre72`;
- `white_blood_cells_missing_pre72` duplicates `creatinine_missing_pre72` in the audited cohort.

Direct `gcs_total_48_72h` and `fio2_48_72h` are also excluded after the targeted vital measurement audit/ablation. The reasons are recorded next to the exclusions in `src/sepsis_deescalation/specification.py`.

Do not add rank pruning, alternate PS algorithms, new missingness indicators, or additional post hoc covariates to the publication model.

## 8. Bootstrap inference

The publication analysis uses 1,000 nonparametric encounter/admission-level bootstrap replicates.

Within each bootstrap replicate, the propensity model is refit. The bootstrap does **not** hold the original weights fixed.

The optimized implementation in `src/sepsis_deescalation/fast_bootstrap.py` uses a numeric design matrix and caps the ordinary GLM attempt at 150 iterations. The source comments explain why: a matched-sample benchmark showed parity with the historical 200-iteration rule at 150, while lower caps did not. The ridge fallback remains unchanged.

This is a runtime optimization, not a change in estimand.

## 9. Weighting sensitivities

The primary result remains stabilized IPTW ATE.

Prespecified/frozen sensitivity approaches include:

- overlap weighting, which targets the overlap population (ATO), not the same ATE;
- stabilized-weight truncation at 1/99 and 2.5/97.5 percentiles.

The overlap result has better measured balance but must not replace the primary result post hoc because it changes the target population.

## 10. Outcomes

Primary:

- 30-day all-cause mortality from the 96-hour landmark.

Secondary:

- hospital-free days, death assigned 0;
- antibiotic-free days, death assigned 0;
- normalized systemic antibiotic exposure;
- normalized broad-spectrum exposure.

Exploratory:

- late recurrent/persistent antibiotic-course use.

Antibiotic-burden outcomes are mechanically downstream of treatment classification. They demonstrate treatment separation/burden and should not be interpreted as independent causal clinical benefits. The late-course outcome is especially vulnerable to differential observation after discharge.

## 11. Exact publication reproduction sequence

After installing dependencies and validating MIMIC-IV 3.1:

```bash
python scripts/run_mimic.py \
  --config config/mimic.yaml \
  --mode final \
  --jobs auto
```

Record the timestamped run directory printed by the command, for example:

```bash
export RUN_DIR=outputs/mimic/mimic_iv_v5_7_final_YYYYMMDDTHHMMSSZ
```

Apply the audited vital correction:

```bash
python scripts/repair_v57_vital_covariates.py "$RUN_DIR" \
  --config config/mimic.yaml
```

Then rerun final statistical inference from the corrected cohort:

```bash
python scripts/rerun_inference.py "$RUN_DIR" \
  --config config/mimic.yaml \
  --mode final \
  --jobs auto \
  --cohort-path "$RUN_DIR/audits/vital_repair/analysis_cohort_vital_corrected.csv" \
  --label vital_corrected_final
```

The manuscript primary/secondary, progressive-adjustment, and final weighting numbers come from the corrected `final_vital_corrected_final_*` inference rerun, not from the pre-repair point estimates in the base run directory.

Source-dependent microbiology and missing-stop-time sensitivities still require the complete base run because an inference-only rerun cannot reconstruct source-level eligibility changes.

## 12. Frozen publication values for parity review

A reviewer should expect the corrected primary mortality result to reproduce approximately:

- de-escalated/stopped weighted risk: `0.1829947`;
- continued weighted risk: `0.1746189`;
- RD: `+0.0083758`;
- RR: `1.0479663`;
- RD bootstrap 95% CI: `[-0.0439033, +0.0567136]`;
- max post-weighting absolute SMD: approximately `0.13261`;
- treated ESS: approximately `613.96`;
- continued ESS: approximately `7119.58`;
- maximum weight: approximately `30.86`.

The exact reporting rules and the audited M4-versus-primary bootstrap-CI provenance issue are documented in `docs/mimic_v57_freeze_review.md`.

## 13. What not to change during manuscript preparation

Do not modify the frozen publication analysis because a different specification produces a more favorable treatment effect. In particular, do not add new subgroups, interactions, ML propensity models, trimming thresholds, missingness indicators, or covariates after inspecting treatment effects unless the work is explicitly labeled as a new post hoc analysis.
