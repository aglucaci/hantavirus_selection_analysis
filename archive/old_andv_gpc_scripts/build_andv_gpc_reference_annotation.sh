#!/usr/bin/env bash
set -euo pipefail

python scripts/build_andv_gpc_reference_annotation.py \
  --m-fasta data/andv_nuc-M_2026-05-11T1953.fasta/andv_nuc-M_2026-05-11T1953.fasta \
  --cds-fasta results/GPC_CDS.recovered.raw.fasta \
  --ref-id PP_006W0E7.1 \
  --outdir data/reference \
  --verbose
