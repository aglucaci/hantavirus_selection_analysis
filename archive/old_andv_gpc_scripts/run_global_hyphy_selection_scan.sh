#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Global ANDV GPC codon selection scan with HyPhy.
#
# First-pass, no-host-label baseline analysis:
#   - FEL
#   - MEME
#   - BUSTED
#   - aBSREL
#
# Expected inputs:
#   results/GPC_CDS.mafft_codon_aligned.deduplicated.fasta
#   results/iqtree2/GPC_CDS.mafft_codon_aligned.deduplicated.treefile
###############################################################################

ALIGNMENT="${ALIGNMENT:-results/GPC_CDS.mafft_codon_aligned.deduplicated.fasta}"
TREE="${TREE:-results/iqtree2/GPC_CDS.mafft_codon_aligned.deduplicated.treefile}"
OUTDIR="${OUTDIR:-results/hyphy_global_selection}"
PREFIX="${PREFIX:-ANDV_GPC_global}"
CLEAN_INPUTS="${CLEAN_INPUTS:-1}"
REFERENCE_ID="${REFERENCE_ID:-PP_006W0E7_1}"

FEL_OUT="${FEL_OUT:-$OUTDIR/${PREFIX}_FEL.json}"
MEME_OUT="${MEME_OUT:-$OUTDIR/${PREFIX}_MEME.json}"
BUSTED_OUT="${BUSTED_OUT:-$OUTDIR/${PREFIX}_BUSTED.json}"
ABSREL_OUT="${ABSREL_OUT:-$OUTDIR/${PREFIX}_aBSREL.json}"
CLEAN_ALIGNMENT="${CLEAN_ALIGNMENT:-$OUTDIR/${PREFIX}.hyphy_ready.fasta}"
CLEAN_TREE="${CLEAN_TREE:-$OUTDIR/${PREFIX}.hyphy_ready.treefile}"
CLEAN_QC="${CLEAN_QC:-$OUTDIR/${PREFIX}.hyphy_ready.qc.tsv}"

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

run_hyphy_analysis() {
  local label="$1"
  local batch_file="$2"
  local output="$3"

  echo "Running $label with batch file: $batch_file"
  "$HYPHY_EXECUTABLE" "$batch_file" \
    --alignment "$ALIGNMENT" \
    --tree "$TREE" \
    --output "$output"
}

echo "[setup] Checking inputs and tools..."

if [[ ! -s "$ALIGNMENT" ]]; then
  echo "ERROR: codon alignment not found or empty: $ALIGNMENT"
  echo "Run scripts/generate_mafft_codon_alignment.sh first."
  exit 1
fi

if [[ ! -s "$TREE" ]]; then
  echo "ERROR: tree not found or empty: $TREE"
  echo "Run scripts/build_iqtree2_tree_from_mafft_codon_alignment.sh first."
  exit 1
fi

if ! HYPHY_EXECUTABLE="$(find_hyphy_executable)"; then
  echo "ERROR: HyPhy executable not found."
  echo "Install with: conda install -c bioconda hyphy"
  exit 1
fi

FEL_BF="$(find_hyphy_batch_file FEL)"
MEME_BF="$(find_hyphy_batch_file MEME)"
BUSTED_BF="$(find_hyphy_batch_file BUSTED)"
ABSREL_BF="$(find_hyphy_batch_file aBSREL)"

if [[ -z "$FEL_BF" || -z "$MEME_BF" || -z "$BUSTED_BF" || -z "$ABSREL_BF" ]]; then
  echo "ERROR: could not find one or more HyPhy SelectionAnalyses batch files."
  echo "FEL:    ${FEL_BF:-missing}"
  echo "MEME:   ${MEME_BF:-missing}"
  echo "BUSTED: ${BUSTED_BF:-missing}"
  echo "aBSREL: ${ABSREL_BF:-missing}"
  exit 1
fi

mkdir -p "$OUTDIR"

if [[ "$CLEAN_INPUTS" == "1" ]]; then
  echo "[clean] Preparing HyPhy-ready codon alignment and tree..."
  python scripts/filter_hyphy_codon_inputs.py \
    --alignment "$ALIGNMENT" \
    --tree "$TREE" \
    --out-alignment "$CLEAN_ALIGNMENT" \
    --out-tree "$CLEAN_TREE" \
    --qc "$CLEAN_QC" \
    --require-id "$REFERENCE_ID" \
    --verbose

  ALIGNMENT="$CLEAN_ALIGNMENT"
  TREE="$CLEAN_TREE"
else
  echo "[clean] Skipping HyPhy input cleaning because CLEAN_INPUTS=0"
fi

echo "[run] Global ANDV GPC selection scan"
echo "Alignment: $ALIGNMENT"
echo "Tree:      $TREE"
echo "Output:    $OUTDIR"
echo "HyPhy:     $HYPHY_EXECUTABLE"
echo
echo "This is the unlabeled baseline scan: FEL + MEME + BUSTED + aBSREL."
echo

echo "[1/4] Running FEL..."
run_hyphy_analysis FEL "$FEL_BF" "$FEL_OUT"

echo
echo "[2/4] Running MEME..."
run_hyphy_analysis MEME "$MEME_BF" "$MEME_OUT"

echo
echo "[3/4] Running BUSTED..."
run_hyphy_analysis BUSTED "$BUSTED_BF" "$BUSTED_OUT"

echo
echo "[4/4] Running aBSREL..."
run_hyphy_analysis aBSREL "$ABSREL_BF" "$ABSREL_OUT"

echo
echo "[done]"
echo "FEL:     $FEL_OUT"
echo "MEME:    $MEME_OUT"
echo "BUSTED:  $BUSTED_OUT"
echo "aBSREL:  $ABSREL_OUT"
echo
echo "First report framing:"
echo "We first performed a global codon-level selection scan of high-quality ANDV"
echo "GPC sequences across host sources. Because GPC mediates viral entry and"
echo "contains the Gn/Gc envelope glycoproteins, this scan tested whether the"
echo "protein is dominated by purifying selection while retaining evidence of"
echo "episodic diversifying selection at specific codons or branches."
