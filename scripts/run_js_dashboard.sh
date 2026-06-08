#!/usr/bin/env bash
set -euo pipefail

echo "Building dashboard tables from ${HYPHY_RESULTS_DIR:-results/ANDV_trees_aln-hyphy}"
python scripts/build_hyphy_dashboard_tables.py

HOST="${DASHBOARD_HOST:-127.0.0.1}"
PORT="${DASHBOARD_PORT:-8502}"

while python - "$HOST" "$PORT" <<'PY'
import socket
import sys

host = sys.argv[1]
port = int(sys.argv[2])
with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
    sys.exit(0 if sock.connect_ex((host, port)) == 0 else 1)
PY
do
  echo "Port $PORT is already in use; trying $((PORT + 1))"
  PORT=$((PORT + 1))
done

echo "Serving dashboard at http://$HOST:$PORT/dashboard-js/"
python -m http.server "$PORT" --bind "$HOST"
