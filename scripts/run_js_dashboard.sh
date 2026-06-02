#!/usr/bin/env bash
set -euo pipefail

python scripts/build_hyphy_dashboard_tables.py
python -m http.server 8502 --bind 127.0.0.1
