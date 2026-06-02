#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

echo "[dashboard] Building dashboard-ready TSV/JSON files..."
"${PYTHON:-python}" scripts/build_dashboard_data.py

echo "[dashboard] Starting Streamlit dashboard..."
echo "[dashboard] URL: http://localhost:8501"
"${PYTHON:-python}" -m streamlit run dashboard/andv_gpc_dashboard.py --server.address 0.0.0.0 --server.port 8501
