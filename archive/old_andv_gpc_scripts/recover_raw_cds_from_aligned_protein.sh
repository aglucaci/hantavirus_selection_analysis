#!/usr/bin/env bash
set -euo pipefail

python scripts/recover_raw_cds_from_aligned_protein.py \
  -p data/andv_aligned-aa-GPC_2026-05-11T2017.fasta/andv_aligned-aa-GPC_2026-05-11T2017.fasta \
  -n data/andv_nuc-M_2026-05-11T1953.fasta/andv_nuc-M_2026-05-11T1953.fasta \
  -o results/GPC_CDS.recovered.raw.fasta \
  --qc results/GPC_CDS.recovered.raw.qc.tsv \
  --min-identity 0.60 \
  --verbose
