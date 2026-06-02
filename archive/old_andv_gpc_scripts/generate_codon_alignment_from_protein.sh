#!/usr/bin/env bash
set -euo pipefail

python scripts/generate_codon_alignment_from_protein.py \
  -p data/andv_aligned-aa-GPC_2026-05-11T2017.fasta/andv_aligned-aa-GPC_2026-05-11T2017.fasta \
  -c results/GPC_CDS.recovered.raw.fasta \
  -o results/GPC_CDS.codon_aligned.fasta \
  --qc results/GPC_CDS.codon_aligned.qc.tsv \
  --min-translation-identity 0.70 \
  --missing-as gap \
  --verbose
