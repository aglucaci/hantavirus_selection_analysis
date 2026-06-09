#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
REPO_ROOT="$(cd -- "$SCRIPT_DIR/.." >/dev/null 2>&1 && pwd)"
cd "$REPO_ROOT"

find_python() {
  if [[ -n "${PYTHON:-}" ]] && command -v "$PYTHON" >/dev/null 2>&1; then
    printf '%s\n' "$PYTHON"
    return 0
  fi
  if command -v python3 >/dev/null 2>&1; then
    printf '%s\n' "python3"
    return 0
  fi
  if command -v python >/dev/null 2>&1; then
    printf '%s\n' "python"
    return 0
  fi
  if command -v conda >/dev/null 2>&1 && conda env list | awk '{print $1}' | grep -qx "hantavirus-snake"; then
    printf '%s\n' "conda run -n hantavirus-snake python"
    return 0
  fi
  return 1
}

PYTHON_CMD="$(find_python || true)"
if [[ -z "$PYTHON_CMD" ]]; then
  cat >&2 <<'EOF'
Could not find Python.

On Ubuntu/WSL, install Python or create the conda environment:
  sudo apt-get update && sudo apt-get install -y python3
  conda env create -f environment.yml

Then rerun:
  bash scripts/run_js_dashboard.sh
EOF
  exit 1
fi

run_python() {
  # shellcheck disable=SC2086
  $PYTHON_CMD "$@"
}

echo "Using Python: $PYTHON_CMD"
run_python - <<'PY'
import sys
if sys.version_info < (3, 10):
    raise SystemExit(f"Python >= 3.10 is required; found {sys.version.split()[0]}")
PY

for required_path in \
  "dashboard-js/index.html" \
  "dashboard-js/src/app.js" \
  "dashboard-js/src/styles.css" \
  "scripts/build_hyphy_dashboard_tables.py"
do
  if [[ ! -f "$required_path" ]]; then
    echo "Missing required dashboard file: $required_path" >&2
    exit 1
  fi
done

RESULTS_DIR="${HYPHY_RESULTS_DIR:-results/ANDV_trees_aln-hyphy}"
if [[ ! -d "$RESULTS_DIR" && "$RESULTS_DIR" == "results/ANDV_trees_aln-hyphy" && -d "results/ANDV_trees-aln-hyphy" ]]; then
  RESULTS_DIR="results/ANDV_trees-aln-hyphy"
  export HYPHY_RESULTS_DIR="$RESULTS_DIR"
elif [[ ! -d "$RESULTS_DIR" && "$RESULTS_DIR" == "results/ANDV_trees-aln-hyphy" && -d "results/ANDV_trees_aln-hyphy" ]]; then
  RESULTS_DIR="results/ANDV_trees_aln-hyphy"
  export HYPHY_RESULTS_DIR="$RESULTS_DIR"
fi

if [[ ! -d "$RESULTS_DIR" ]]; then
  echo "Results directory not found: $RESULTS_DIR" >&2
  echo "Set HYPHY_RESULTS_DIR to the HyPhy results folder, or run the analyses first." >&2
  exit 1
fi

echo "Building dashboard tables from $RESULTS_DIR"
run_python scripts/build_hyphy_dashboard_tables.py

HOST="${DASHBOARD_HOST:-127.0.0.1}"
PORT="${DASHBOARD_PORT:-8502}"

while true
do
  PORT_STATUS="$(run_python - "$HOST" "$PORT" <<'PY'
import socket
import sys

host = sys.argv[1]
port = int(sys.argv[2])
with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
    sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
    try:
        sock.bind((host, port))
    except OSError as exc:
        if exc.errno in (48, 98):
            print("in_use")
        else:
            print(f"unavailable:{exc}")
        sys.exit(0)
    print("available")
PY
)"
  if [[ "$PORT_STATUS" == "available" ]]; then
    break
  fi
  if [[ "$PORT_STATUS" == "in_use" ]]; then
    echo "Port $PORT is already in use; trying $((PORT + 1))"
    PORT=$((PORT + 1))
    continue
  fi
  echo "Cannot bind dashboard server to $HOST:$PORT ($PORT_STATUS)" >&2
  echo "Try another port or host, for example:" >&2
  echo "  DASHBOARD_PORT=8510 bash scripts/run_js_dashboard.sh" >&2
  echo "  DASHBOARD_HOST=0.0.0.0 bash scripts/run_js_dashboard.sh" >&2
  exit 1
done

if [[ "$HOST" == "0.0.0.0" || "$HOST" == "::" ]]; then
  URL_HOST="127.0.0.1"
else
  URL_HOST="$HOST"
fi

echo "Serving dashboard at http://$URL_HOST:$PORT/dashboard-js/"
echo "Press Ctrl-C to stop the server."
run_python -m http.server "$PORT" --bind "$HOST"
