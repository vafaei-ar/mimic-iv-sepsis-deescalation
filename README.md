# MIMIC-IV Sepsis Day-3 De-escalation

Reproducible analysis code for the MIMIC-IV day-3 broad-spectrum antibiotic de-escalation target-trial emulation and the Penn State modified external replication.

## Publication status

The scientific definitions and publication analyses are frozen. `main` is the canonical branch for manuscript-facing code and documentation. New scientific specifications should be treated as post hoc analyses rather than folded into the frozen primary analysis.

## Start here

For scientific or code review, read:

1. `docs/target_trial_spec.md` - frozen target-trial contract.
2. `docs/mimic_analysis_walkthrough.md` - final MIMIC workflow, corrected-vital publication path, and parity targets.
3. `docs/psu_crosswalk.md` - MIMIC-to-Penn-State data-source crosswalk.
4. `docs/psu_analysis_walkthrough.md` - final Penn State modified-replication workflow and parity targets.
5. `docs/mimic_v57_freeze_review.md` - final MIMIC publication-freeze record.

## Scientific design

### MIMIC-IV primary analysis

The primary MIMIC-IV analysis is anchored at the first qualifying systemic intravenous broad-spectrum antibiotic exposure. The treatment decision occurs 72 hours later, treatment is classified during hours 72-96, and follow-up begins at the 96-hour landmark. The primary culture-negative definition requires qualifying microbiology sampling and no positive clinical culture result available by the 72-hour decision.

The primary estimand is the average treatment effect estimated with stabilized inverse probability treatment weighting. Covariates are restricted to information available before the 72-hour decision. The final model includes near-decision clinical status, recovery trajectories, diagnostic intensity, and treatment intensity.

### Penn State modified external replication

Penn State preserves the conceptual first-broad-spectrum anchor, 72-hour decision, 72-96-hour classification window, and 96-hour landmark, but it is not an exact MIMIC replication. Source semantics differ for hospital/ICU representation, medication exposure, microbiology, route information, and event-time resolution.

The primary Penn State medication phenotype uses order-based prescribing records as the closest defensible analogue to the MIMIC prescription/order construct. Medication-administration records are used as a prespecified sensitivity analysis. The frozen Penn State analytic cohort is drawn from the upstream sepsis encounter source and uses the adapted local sepsis definition documented in the Penn State walkthrough.

## Repository layout

```text
config/                  analysis configuration and site-mapping templates
src/sepsis_deescalation/ reusable MIMIC analysis package
scripts/                 final entry points plus retained audit/provenance scripts
tests/                   unit, smoke, and publication-contract tests
docs/                    scientific contracts, freeze records, and reviewer guides
outputs/                 generated local results and caches; ignored by git
.runrelay/                approved-machine execution manifest
```

## Final analysis entry points

### MIMIC-IV

```bash
python scripts/run_mimic.py --config config/mimic.yaml --mode final --jobs auto

python scripts/repair_v57_vital_covariates.py "$RUN_DIR" --config config/mimic.yaml

python scripts/rerun_inference.py "$RUN_DIR" \
  --config config/mimic.yaml \
  --mode final \
  --jobs auto \
  --cohort-path "$RUN_DIR/audits/vital_repair/analysis_cohort_vital_corrected.csv" \
  --label vital_corrected_final
```

The manuscript primary and secondary outcomes, progressive-adjustment sequence, and final weighting diagnostics come from the corrected vital-sign inference rerun. Source-dependent microbiology and missing-stop-time sensitivities remain tied to the complete base run.

### Penn State

The frozen Penn State publication sequence is:

```bash
python scripts/audit_psu_final_covariate_freeze.py "$PSU_DATA_ROOT" \
  --output-dir outputs/psu_final_covariate_freeze/latest

python scripts/audit_psu_ps_balance.py "$PSU_DATA_ROOT" \
  --output-dir outputs/psu_ps_balance/latest

python scripts/audit_psu_final_outcome_freeze.py "$PSU_DATA_ROOT" \
  --output-dir outputs/psu_final_outcome_freeze/latest

python scripts/run_psu_point_estimates.py "$PSU_DATA_ROOT" \
  --output-dir outputs/psu_point_estimates/latest

python scripts/run_psu_bootstrap_inference.py "$PSU_DATA_ROOT" \
  --output-dir outputs/psu_bootstrap_inference/latest

python scripts/run_psu_prespecified_robustness.py "$PSU_DATA_ROOT" \
  --output-dir outputs/psu_prespecified_robustness/latest

python scripts/run_psu_prespecified_robustness_bootstrap.py "$PSU_DATA_ROOT" \
  --output-dir outputs/psu_prespecified_robustness_bootstrap/latest
```

Expected parity targets and the rationale for each stage are documented in `docs/psu_analysis_walkthrough.md`.

## Publication-support builders

Manuscript-facing aggregate outputs are built from frozen analyses. Important entry points include:

- `scripts/build_nature_figures.py` - final manuscript and ESM figures.
- `scripts/build_publication_baseline_characteristics.py` - combined MIMIC-IV and Penn State baseline-characteristics aggregates.
- `scripts/build_publication_integration.py` - harmonized manuscript-facing aggregate tables.
- `scripts/build_manuscript_package.py` - publication-support text and figure/table data.

These builders export aggregate or sanitized publication outputs only.

## Historical audit scripts

The `audit_psu_*` scripts are retained intentionally as provenance. They document source-semantic checks performed before the Penn State phenotype, covariates, and outcomes were frozen. They are not all required for routine reproduction, but deleting them would remove the audit trail for non-obvious design decisions. The reviewer walkthrough distinguishes final entry points from historical audits.

## Data and privacy policy

Never commit patient-level MIMIC-IV, Penn State, PCORnet, or derived analytic data. Raw and patient-level data remain local. `outputs/` is ignored except for its placeholder file. RunRelay artifacts are limited to explicitly declared safe aggregate outputs.

The public repository contains code, configuration templates, tests, documentation, and non-sensitive aggregate publication support only.

## Installation and validation

Use a repository-local environment:

```bash
/usr/bin/python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e '.[dev]'
pytest -q
ruff check .
```

## Reproducibility rules

1. Work from a named git commit and a clean repository state.
2. Do not edit generated effect tables by hand.
3. Keep bootstrap seeds and publication replicate counts fixed.
4. Freeze scientific definitions before examining treatment-effect changes.
5. Document site-specific differences rather than silently forcing Penn State into MIMIC semantics.
6. Keep patient-level checkpoints and source data local and outside git or artifact transport.
7. Require parity on cohort counts, balance, weights, outcomes, point estimates, and bootstrap intervals before replacing frozen publication code.
8. Treat overlap weighting as a different estimand rather than a replacement for the primary ATE.
9. Describe Penn State as a modified external replication rather than an exact validation cohort.
