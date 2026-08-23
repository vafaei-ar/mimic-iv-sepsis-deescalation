# MIMIC-IV Sepsis Day-3 De-escalation

Reproducible analysis code for the MIMIC-IV day-3 broad-spectrum antibiotic de-escalation study and its harmonized Penn State/PCORnet external validation.

## Design

The primary MIMIC analysis is anchored at the first qualifying systemic IV broad-spectrum antibiotic exposure. The treatment decision time is 72 h later, treatment is classified over 72-96 h, and follow-up starts at the 96-h landmark. The primary culture-negative definition requires qualifying microbiology sampling by 72 h and no positive organism result available to clinicians by that decision time. The prior eventual-culture-negative definition is retained as a sensitivity analysis.

The external validation code uses the same conceptual clock, eligibility, treatment strategies, outcomes, and weighting framework. Site-specific source mappings live in configuration, not inside the statistical code.

## Repository layout

```text
config/                 analysis configuration and site mapping templates
src/sepsis_deescalation reusable analysis package
scripts/                command-line entry points
tests/                  unit/smoke tests
docs/                   target-trial and MIMIC-to-PCORnet crosswalks
outputs/                 generated results and local caches; ignored by git
```

## Data policy

**Never commit patient-level MIMIC-IV, PSU, PCORnet, or derived analytic data.** Patient-level CSV/Parquet files and `outputs/` are ignored by git. The repository should contain code, configuration templates, tests, documentation, and non-sensitive aggregate outputs only.

## Installation

Use a repository-local virtual environment:

```bash
/usr/bin/python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e '.[dev]'
pytest -q
```

## MIMIC-IV run

Set the MIMIC-IV v3.1 root if the configured local path is not already available:

```bash
export MIMIC_SOURCE=/home/asadr/datasets/MIMIC/physionet.org/files/mimiciv/3.1
```

Run a validation/preflight first:

```bash
python scripts/validate_mimic.py --config config/mimic.yaml
```

### Fast development mode

Use this after code changes. It keeps the same cohort and estimands but limits each bootstrap analysis to the configured development count (currently 100) and parallelizes replicates:

```bash
python scripts/run_mimic.py \
  --config config/mimic.yaml \
  --mode fast \
  --jobs auto
```

### Final publication mode

The final mode uses the configured publication bootstrap counts (currently 1,000):

```bash
python scripts/run_mimic.py \
  --config config/mimic.yaml \
  --mode final \
  --jobs auto
```

`--jobs auto` uses up to eight worker processes. BLAS/OpenMP libraries are restricted to one thread inside each bootstrap worker to avoid CPU oversubscription.

The optimized bootstrap engine uses a numeric design matrix instead of rebuilding Patsy formulas for every replicate. The six primary/secondary outcomes share one propensity-score fit per bootstrap replicate, and overlap/truncation weighting sensitivities also share one fit per replicate. Point-estimate definitions and target estimands are unchanged.

### Inference-only resume mode

Each complete run writes a local patient-level Parquet checkpoint under `outputs/cache/mimic/<analysis_version>/`. For later statistical development, primary/secondary, progressive, and weighting inference can be rerun without rereading MIMIC or rebuilding features:

```bash
python scripts/rerun_inference.py \
  outputs/mimic/mimic_iv_v5_7_final_YYYYMMDDTHHMMSSZ \
  --config config/mimic.yaml \
  --mode fast \
  --jobs auto
```

This is intentionally an inference-only checkpoint. It does not replace a complete final MIMIC run for source-dependent microbiology and missing-stop-time sensitivity analyses.

Each run creates a timestamped directory under `outputs/mimic/`, including a manifest, logs, cohort-flow counts, estimates, bootstrap CIs, balance diagnostics, sensitivity analyses, figure-ready CSVs, runtime timings, and a ZIP archive for review.

## PSU / PCORnet external validation

Start by copying the mapping template and adapting only the source paths/fields that differ at Penn State:

```bash
cp config/pcornet_psu.example.yaml config/pcornet_psu.local.yaml
python scripts/validate_pcornet.py --config config/pcornet_psu.local.yaml
python scripts/run_pcornet.py --config config/pcornet_psu.local.yaml
```

Local configuration files are ignored by git. The PSU code is designed to support both `PRESCRIBING` and `MED_ADMIN`. The primary harmonized replication can use the source closest to the frozen MIMIC exposure definition, while a second administration-based analysis can quantify treatment reclassification.

## Reproducibility rules

1. Run from a clean `.venv` environment/process.
2. Do not edit generated CSV estimates manually.
3. Keep random seeds in configuration.
4. Store a machine-readable run manifest with package versions, git commit, configuration hash, and timestamps.
5. Freeze MIMIC definitions before changing the PSU implementation.
6. Site-specific differences must be documented in `docs/psu_crosswalk.md`.
7. Use `--mode fast` for development and `--mode final` only for publication-quality inference.
