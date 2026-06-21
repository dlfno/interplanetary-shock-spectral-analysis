# Spectral analysis of an interplanetary shock — convenience targets
PY ?= python3
VENV := .venv
BIN := $(VENV)/bin

.PHONY: help setup run synthetic test clean clean-data

help:
	@echo "make setup      - create .venv and install requirements"
	@echo "make run        - download data (cached) and generate all figures"
	@echo "make synthetic  - run the pipeline with the synthetic fallback"
	@echo "make test       - run the test-suite (pytest)"
	@echo "make clean      - remove generated figures"
	@echo "make clean-data - also remove the downloaded data cache"

setup:
	$(PY) -m venv $(VENV)
	$(BIN)/pip install --upgrade pip
	$(BIN)/pip install -r requirements.txt

run:
	$(BIN)/python scripts/run_analysis.py

synthetic:
	$(BIN)/python scripts/run_analysis.py --synthetic

test:
	$(BIN)/python -m pytest -q

clean:
	rm -f figures/*.png figures/results.md results.json

clean-data:
	rm -f data/*.cdf data/*.npz
