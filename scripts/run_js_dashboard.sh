#!/usr/bin/env bash
set -euo pipefail

python scripts/build_hyphy_dashboard_tables.py
echo "Serving dashboard at http://127.0.0.1:8502/dashboard-js/"
python -m http.server 8502 --bind 127.0.0.1
