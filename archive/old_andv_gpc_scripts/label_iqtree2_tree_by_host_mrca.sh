#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Label an IQ-TREE2 Newick tree at the MRCA of hostNameScientific groups.
###############################################################################

TREE="${TREE:-results/iqtree2/GPC_CDS.mafft_codon_aligned.deduplicated.treefile}"
METADATA="${METADATA:-data/andv_metadata_2026-05-11T1953.tsv/andv_metadata_2026-05-11T1953.tsv}"
OUT_TREE="${OUT_TREE:-results/iqtree2/GPC_CDS.mafft_codon_aligned.deduplicated.host_mrca_labeled.treefile}"
QC="${QC:-results/iqtree2/GPC_CDS.mafft_codon_aligned.deduplicated.host_mrca_labeled.qc.tsv}"
MIN_TIPS="${MIN_TIPS:-2}"
BRANCH_LABEL="${BRANCH_LABEL:-}"

echo "[setup] Checking inputs..."

if [[ ! -s "$TREE" ]]; then
  echo "ERROR: tree not found or empty: $TREE"
  echo "Run scripts/build_iqtree2_tree_from_mafft_codon_alignment.sh first."
  exit 1
fi

if [[ ! -s "$METADATA" ]]; then
  echo "ERROR: metadata TSV not found or empty: $METADATA"
  exit 1
fi

mkdir -p "$(dirname "$OUT_TREE")" "$(dirname "$QC")"

echo "[run] Labeling host MRCA branches..."

BRANCH_LABEL_ARGS=()
if [[ -n "$BRANCH_LABEL" ]]; then
  BRANCH_LABEL_ARGS=(--branch-label "$BRANCH_LABEL")
fi

python scripts/label_tree_mrca_by_metadata.py \
  -t "$TREE" \
  -m "$METADATA" \
  -o "$OUT_TREE" \
  --qc "$QC" \
  --id-column accessionVersion \
  --label-column hostNameScientific \
  "${BRANCH_LABEL_ARGS[@]}" \
  --min-tips "$MIN_TIPS" \
  --skip-label Unknown_host \
  --verbose

echo
echo "[done]"
echo "Input tree:    $TREE"
echo "Metadata TSV:  $METADATA"
echo "Labeled tree:  $OUT_TREE"
echo "QC report:     $QC"
