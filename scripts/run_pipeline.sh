#!/usr/bin/env bash
# Rebuild 100 bp bins, peak-region figures, mean tracks, and correlation tables from the raw GEO files.
set -euo pipefail
cd "$(dirname "$0")/.."
if [ ! -d .venv ]; then uv venv .venv && uv pip install --python .venv/bin/python -r requirements.txt; fi
PY=.venv/bin/python
( cd analysis/peaks   && $PY 01_bin.py && $PY 02_select_regions.py && $PY 03_extract_windows.py && $PY 04_plot.py )
( cd analysis/targets && $PY 05_targets_and_correlations.py && $PY 06_make_tables.py )
