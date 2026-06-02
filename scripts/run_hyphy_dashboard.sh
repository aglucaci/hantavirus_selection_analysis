#!/usr/bin/env bash
set -euo pipefail

python scripts/build_hyphy_dashboard_tables.py
STREAMLIT_BROWSER_GATHER_USAGE_STATS=false streamlit run dashboard/selection_dashboard.py --server.address 127.0.0.1 --server.port 8501
