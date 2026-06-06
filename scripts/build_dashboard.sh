#!/usr/bin/env bash
set -euo pipefail

echo "Building dashboard tables from ${HYPHY_RESULTS_DIR:-results/ANDV_trees_aln-hyphy}"
python scripts/build_hyphy_dashboard_tables.py
