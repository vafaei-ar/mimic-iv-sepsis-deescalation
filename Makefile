PYTHON ?= python
MIMIC_CONFIG ?= config/mimic.yaml
PCORNET_CONFIG ?= config/pcornet_psu.local.yaml

.PHONY: install test lint validate-mimic mimic validate-pcornet pcornet

install:
	$(PYTHON) -m pip install -e '.[dev]'

test:
	$(PYTHON) -m pytest -q

lint:
	$(PYTHON) -m ruff check src scripts tests

validate-mimic:
	$(PYTHON) scripts/validate_mimic.py --config $(MIMIC_CONFIG)

mimic:
	$(PYTHON) scripts/run_mimic.py --config $(MIMIC_CONFIG)

validate-pcornet:
	$(PYTHON) scripts/validate_pcornet.py --config $(PCORNET_CONFIG)

pcornet:
	$(PYTHON) scripts/run_pcornet.py --config $(PCORNET_CONFIG)
