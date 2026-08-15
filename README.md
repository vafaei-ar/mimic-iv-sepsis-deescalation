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
outputs/                 generated aggregate results only; ignored by git
```

## Data policy

**Never commit patient-level MIMIC-IV, PSU, PCORnet, or derived analytic data.** This repository contains code, configuration templates, tests, documentation, and non-sensitive aggregate outputs only.

## Installation

```bash
conda create -n sepsis-deescalation python=3.11 -y
conda activate sepsis-deescalation
pip install -e .
```

For development/tests:

```bash
pip install -e '.[dev]'
pytest -q
```

## MIMIC-IV run

Set the MIMIC-IV v3.1 root. Your current local path can be used directly:

```bash
export MIMIC_SOURCE=/home/asadr/datasets/MIMIC/physionet.org/files/mimiciv/3.1
```

Run a validation/preflight first:

```bash
python scripts/validate_mimic.py --config config/mimic.yaml
```

Then run the analysis:

```bash
python scripts/run_mimic.py --config config/mimic.yaml
```

or:

```bash
make mimic
```

Each run creates a timestamped directory under `outputs/mimic/`, including a manifest, logs, cohort-flow counts, estimates, bootstrap CIs, balance diagnostics, sensitivity analyses, figure-ready CSVs, and a ZIP archive for review.

## PSU / PCORnet external validation

Start by copying the mapping template and adapting only the source paths/fields that differ at Penn State:

```bash
cp config/pcornet_psu.example.yaml config/pcornet_psu.local.yaml
python scripts/validate_pcornet.py --config config/pcornet_psu.local.yaml
python scripts/run_pcornet.py --config config/pcornet_psu.local.yaml
```

Local configuration files are ignored by git. The PSU code is designed to support both `PRESCRIBING` and `MED_ADMIN`. The primary harmonized replication can use the source closest to the frozen MIMIC exposure definition, while a second administration-based analysis can quantify treatment reclassification.

## Reproducibility rules

1. Run from a clean environment/process.
2. Do not edit generated CSV estimates manually.
3. Keep random seeds in configuration.
4. Store a machine-readable run manifest with package versions, git commit, configuration hash, and timestamps.
5. Freeze MIMIC definitions before changing the PSU implementation.
6. Site-specific differences must be documented in `docs/psu_crosswalk.md`.
