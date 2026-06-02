#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Build a phylogenetic tree from the MAFFT codon alignment using IQ-TREE2.
#
# Expected input:
#   results/GPC_CDS.mafft_codon_aligned.deduplicated.fasta
#
# Main output:
#   results/iqtree2/GPC_CDS.mafft_codon_aligned.deduplicated.treefile
###############################################################################

ALIGNMENT="${ALIGNMENT:-results/GPC_CDS.mafft_codon_aligned.deduplicated.fasta}"
OUTDIR="${OUTDIR:-results/iqtree2}"
PREFIX="${PREFIX:-GPC_CDS.mafft_codon_aligned.deduplicated}"
THREADS="${THREADS:-AUTO}"
MODEL="${MODEL:-MFP}"
BOOTSTRAPS="${BOOTSTRAPS:-1000}"
ALRT="${ALRT:-1000}"
EXTRA_ARGS="${EXTRA_ARGS:-}"
REDO="${REDO:-1}"

echo "[setup] Checking inputs and tools..."

if [[ ! -s "$ALIGNMENT" ]]; then
  echo "ERROR: deduplicated MAFFT codon alignment not found or empty: $ALIGNMENT"
  echo "Run scripts/generate_mafft_codon_alignment.sh first."
  exit 1
fi

if ! command -v iqtree2 >/dev/null 2>&1; then
  echo "ERROR: iqtree2 not found in PATH."
  echo "Install with: conda install -c bioconda iqtree"
  exit 1
fi

mkdir -p "$OUTDIR"

echo "[run] Building tree with IQ-TREE2..."
echo "Alignment:  $ALIGNMENT"
echo "Output dir: $OUTDIR"
echo "Prefix:     $PREFIX"
echo "Model:      $MODEL"
echo "Threads:    $THREADS"
echo "UFBoot:     $BOOTSTRAPS"
echo "SH-aLRT:    $ALRT"
echo "Redo:       $REDO"

REDO_ARGS=()
if [[ "$REDO" == "1" ]]; then
  REDO_ARGS=(-redo)
fi

iqtree2 \
  -s "$ALIGNMENT" \
  -m "$MODEL" \
  -B "$BOOTSTRAPS" \
  -alrt "$ALRT" \
  -T "$THREADS" \
  --prefix "$OUTDIR/$PREFIX" \
  "${REDO_ARGS[@]}" \
  $EXTRA_ARGS

echo
echo "[done]"
echo "Tree file:        $OUTDIR/$PREFIX.treefile"
echo "IQ-TREE report:   $OUTDIR/$PREFIX.iqtree"
echo "Log file:         $OUTDIR/$PREFIX.log"
echo "Model file:       $OUTDIR/$PREFIX.model.gz"
