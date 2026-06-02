#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Run HyPhy RELAX comparing human-derived vs reservoir-derived terminal branches.
#
# Defaults:
#   Test      = Homo sapiens tips labeled {Test}
#   Reference = known non-human host tips labeled {Reference}
#   Unknown/missing host tips are left unlabeled.
###############################################################################

ALIGNMENT="${ALIGNMENT:-results/hyphy_global_selection/ANDV_GPC_global.hyphy_ready.fasta}"
TREE="${TREE:-results/hyphy_global_selection/ANDV_GPC_global.hyphy_ready.treefile}"
METADATA="${METADATA:-data/andv_metadata_2026-05-11T1953.tsv/andv_metadata_2026-05-11T1953.tsv}"
DUPLICATE_MAP="${DUPLICATE_MAP:-results/GPC_CDS.mafft_codon_aligned.duplicates.tsv}"
OUTDIR="${OUTDIR:-results/hyphy_host_relax}"
PREFIX="${PREFIX:-ANDV_GPC_human_vs_reservoir_RELAX}"
TEST_HOST="${TEST_HOST:-Homo sapiens}"
TEST_LABEL="${TEST_LABEL:-Test}"
REFERENCE_LABEL="${REFERENCE_LABEL:-Reference}"

RELAX_TREE="${RELAX_TREE:-$OUTDIR/${PREFIX}.treefile}"
RELAX_QC="${RELAX_QC:-$OUTDIR/${PREFIX}.qc.tsv}"
RELAX_OUT="${RELAX_OUT:-$OUTDIR/${PREFIX}.json}"

find_hyphy_executable() {
  if [[ -n "${HYPHY_BIN:-}" ]]; then
    echo "$HYPHY_BIN"
    return
  fi
  if command -v hyphy-avx >/dev/null 2>&1; then
    command -v hyphy-avx
    return
  fi
  if [[ -x /usr/lib/hyphy/bin/hyphy-avx ]]; then
    echo "/usr/lib/hyphy/bin/hyphy-avx"
    return
  fi
  if command -v hyphy >/dev/null 2>&1; then
    command -v hyphy
    return
  fi
  return 1
}

find_hyphy_batch_file() {
  local analysis="$1"
  local candidate
  for candidate in \
    "/usr/share/hyphy/TemplateBatchFiles/SelectionAnalyses/${analysis}.bf" \
    "/usr/lib/hyphy/TemplateBatchFiles/SelectionAnalyses/${analysis}.bf" \
    "/usr/local/share/hyphy/TemplateBatchFiles/SelectionAnalyses/${analysis}.bf"
  do
    if [[ -s "$candidate" ]]; then
      echo "$candidate"
      return
    fi
  done
  find /usr/share/hyphy /usr/lib/hyphy /usr/local/share/hyphy \
    -path "*/SelectionAnalyses/${analysis}.bf" \
    -print -quit 2>/dev/null || true
}

echo "[setup] Checking RELAX inputs..."

if [[ ! -s "$ALIGNMENT" ]]; then
  echo "ERROR: HyPhy-ready alignment not found: $ALIGNMENT"
  echo "Run scripts/run_global_hyphy_selection_scan.sh first."
  exit 1
fi

if [[ ! -s "$TREE" ]]; then
  echo "ERROR: HyPhy-ready tree not found: $TREE"
  echo "Run scripts/run_global_hyphy_selection_scan.sh first."
  exit 1
fi

if [[ ! -s "$METADATA" ]]; then
  echo "ERROR: metadata not found: $METADATA"
  exit 1
fi

if ! HYPHY_EXECUTABLE="$(find_hyphy_executable)"; then
  echo "ERROR: HyPhy executable not found."
  exit 1
fi

RELAX_BF="$(find_hyphy_batch_file RELAX)"
if [[ -z "$RELAX_BF" ]]; then
  echo "ERROR: could not find HyPhy RELAX.bf."
  exit 1
fi

mkdir -p "$OUTDIR"

echo "[label] Building RELAX host-group tree..."
python scripts/label_tree_for_relax_host_groups.py \
  --tree "$TREE" \
  --metadata "$METADATA" \
  --output "$RELAX_TREE" \
  --qc "$RELAX_QC" \
  --test-host "$TEST_HOST" \
  --test-label "$TEST_LABEL" \
  --reference-label "$REFERENCE_LABEL" \
  --duplicate-map "$DUPLICATE_MAP" \
  --verbose

echo
echo "[run] HyPhy RELAX"
echo "Alignment: $ALIGNMENT"
echo "Tree:      $RELAX_TREE"
echo "Test:      $TEST_LABEL ($TEST_HOST)"
echo "Reference: $REFERENCE_LABEL (known non-human hosts)"
echo "Output:    $RELAX_OUT"

"$HYPHY_EXECUTABLE" "$RELAX_BF" \
  --alignment "$ALIGNMENT" \
  --tree "$RELAX_TREE" \
  --test "$TEST_LABEL" \
  --reference "$REFERENCE_LABEL" \
  --output "$RELAX_OUT"

echo
echo "[done]"
echo "RELAX JSON: $RELAX_OUT"
echo "Labeled tree: $RELAX_TREE"
echo "QC report: $RELAX_QC"
