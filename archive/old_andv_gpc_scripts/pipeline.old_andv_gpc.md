# ANDV GPC Analysis Pipeline

Run all commands from the repository root.

## 0. Inputs

The current inputs are:

```text
data/andv_aligned-aa-GPC_2026-05-11T2017.fasta/andv_aligned-aa-GPC_2026-05-11T2017.fasta
data/andv_nuc-M_2026-05-11T1953.fasta/andv_nuc-M_2026-05-11T1953.fasta
data/andv_metadata_2026-05-11T1953.tsv/andv_metadata_2026-05-11T1953.tsv
```

Required command-line tools for the full workflow:

```text
mafft
iqtree2
hyphy
```

## 1. Recover Raw GPC CDS

```bash
bash scripts/recover_raw_cds_from_aligned_protein.sh
```

Outputs:

```text
results/GPC_CDS.recovered.raw.fasta
results/GPC_CDS.recovered.raw.qc.tsv
```

Purpose:

```text
Use the aligned GPC protein records and raw M-segment nucleotide records to
recover raw ungapped GPC CDS sequences. Nucleotide IDs ending in |M are matched
to protein IDs without that suffix.
```

## 2. Build Reference Annotation

```bash
bash scripts/build_andv_gpc_reference_annotation.sh
```

Outputs:

```text
data/reference/ANDV_GPC_reference_annotation.tsv
data/reference/NC_003467.2_M.fasta
data/reference/NC_003467.2_GPC_CDS.no_stop.fasta
data/reference/NC_003467.2_GPC_protein.no_stop.fasta
```

Reference sequence:

```text
original ID:      PP_006W0E7.1
safe ID:          PP_006W0E7_1
RefSeq M segment: NC_003467.2
strain/sample:    Chile-9717869
```

Important coordinates:

```text
GPC CDS without stop: nt 52-3465
terminal stop codon:  nt 3466-3468, TAA
WAASA motif:          aa 647-651
Gn first pass:        aa 1-651
Gc first pass:        aa 652-1138
```

## 3. Generate MAFFT Codon Alignment

```bash
bash scripts/generate_mafft_codon_alignment.sh
```

This step:

```text
1. Sanitizes sequence IDs for MAFFT/IQ-TREE/HyPhy.
2. Translates raw CDS to protein.
3. Aligns proteins with MAFFT.
4. Back-translates the MAFFT protein alignment to codons.
5. Deduplicates the codon alignment.
```

Outputs:

```text
results/GPC_CDS.recovered.raw.hyphy_ids.fasta
results/GPC_CDS.recovered.raw.hyphy_ids.map.tsv
results/GPC_protein.from_recovered_CDS.fasta
results/GPC_protein.from_recovered_CDS.mafft.aligned.fasta
results/GPC_CDS.mafft_codon_aligned.fasta
results/GPC_CDS.mafft_codon_aligned.qc.tsv
results/GPC_CDS.mafft_codon_aligned.deduplicated.fasta
results/GPC_CDS.mafft_codon_aligned.duplicates.tsv
```

Downstream analyses use:

```text
results/GPC_CDS.mafft_codon_aligned.deduplicated.fasta
```

Reference handling:

```text
PP_006W0E7_1 is forced to remain as the representative if it is in a duplicate
group.
```

## 4. Build IQ-TREE2 Tree

```bash
bash scripts/build_iqtree2_tree_from_mafft_codon_alignment.sh
```

Default input:

```text
results/GPC_CDS.mafft_codon_aligned.deduplicated.fasta
```

Outputs:

```text
results/iqtree2/GPC_CDS.mafft_codon_aligned.deduplicated.treefile
results/iqtree2/GPC_CDS.mafft_codon_aligned.deduplicated.iqtree
results/iqtree2/GPC_CDS.mafft_codon_aligned.deduplicated.log
```

The script uses `-redo` by default so the tree is rebuilt from the current
deduplicated alignment.

## 5. Run Global HyPhy Baseline Scan

```bash
bash scripts/run_global_hyphy_selection_scan.sh
```

Default inputs:

```text
results/GPC_CDS.mafft_codon_aligned.deduplicated.fasta
results/iqtree2/GPC_CDS.mafft_codon_aligned.deduplicated.treefile
```

Before running HyPhy, this step creates HyPhy-ready files:

```text
results/hyphy_global_selection/ANDV_GPC_global.hyphy_ready.fasta
results/hyphy_global_selection/ANDV_GPC_global.hyphy_ready.treefile
results/hyphy_global_selection/ANDV_GPC_global.hyphy_ready.qc.tsv
```

The HyPhy-ready cleaner:

```text
removes sequences with stop codons
keeps PP_006W0E7_1 as a required reference
keeps bare safe FASTA IDs only
prunes the tree to match the kept alignment
```

HyPhy outputs:

```text
results/hyphy_global_selection/ANDV_GPC_global_FEL.json
results/hyphy_global_selection/ANDV_GPC_global_MEME.json
results/hyphy_global_selection/ANDV_GPC_global_BUSTED.json
results/hyphy_global_selection/ANDV_GPC_global_aBSREL.json
```

This is the first-pass unlabeled baseline:

```text
FEL + MEME + BUSTED + aBSREL
```

Do this before host-labeled tests.

## 6. Optional Host MRCA Tree Labels

Only after the global scan is clean:

```bash
bash scripts/label_iqtree2_tree_by_host_mrca.sh
```

Inputs:

```text
results/iqtree2/GPC_CDS.mafft_codon_aligned.deduplicated.treefile
data/andv_metadata_2026-05-11T1953.tsv/andv_metadata_2026-05-11T1953.tsv
```

Outputs:

```text
results/iqtree2/GPC_CDS.mafft_codon_aligned.deduplicated.host_mrca_labeled.treefile
results/iqtree2/GPC_CDS.mafft_codon_aligned.deduplicated.host_mrca_labeled.qc.tsv
```

This labels branches at the MRCA of hostNameScientific groups.

## 7. Run RELAX By Host Group

After the global baseline scan, compare human-derived terminal branches against
known non-human host terminal branches:

```bash
bash scripts/run_relax_by_host_group.sh
```

Default labels:

```text
Test      = Homo sapiens
Reference = known non-human hosts
Unlabeled = unknown/missing hosts or ambiguous mixed duplicate groups
```

Outputs:

```text
results/hyphy_host_relax/ANDV_GPC_human_vs_reservoir_RELAX.treefile
results/hyphy_host_relax/ANDV_GPC_human_vs_reservoir_RELAX.qc.tsv
results/hyphy_host_relax/ANDV_GPC_human_vs_reservoir_RELAX.json
```

## 8. Optional HyPhy Site Annotation

After extracting selected sites from HyPhy JSON into a TSV with a `hyphy_site`
column:

```bash
python scripts/annotate_hyphy_sites_with_reference.py \
  -i results/hyphy_global_selection/selected_sites.tsv \
  -o results/hyphy_global_selection/selected_sites.annotated.tsv
```

This adds:

```text
reference_codon
reference_nt_start
reference_nt_end
aa_ref
region
```

## 9. Build And Run The Dashboard

After HyPhy has finished, build the dashboard-ready summary tables:

```bash
python scripts/build_dashboard_data.py
```

Outputs:

```text
results/dashboard/overview_metrics.json
results/dashboard/alignment_sequence_qc.tsv
results/dashboard/alignment_site_qc.tsv
results/dashboard/pairwise_identity.tsv
results/dashboard/ANDV_GPC_selection_annotated_sites.tsv
results/dashboard/ANDV_GPC_FEL.all_sites.tsv
results/dashboard/ANDV_GPC_MEME.all_sites.tsv
results/dashboard/ANDV_GPC_BUSTED.summary.tsv
results/dashboard/ANDV_GPC_ABSREL.branches.tsv
results/dashboard/ANDV_GPC_ABSREL.all_branches.tsv
results/dashboard/ANDV_GPC_RELAX.summary.tsv
results/dashboard/metadata_host_counts.tsv
results/dashboard/metadata_country_counts.tsv
results/dashboard/metadata_year_counts.tsv
results/dashboard/selection_input_host_group_counts.tsv
results/dashboard/selection_input_host_group_details.tsv
```

`selection_input_host_group_counts.tsv` summarizes host composition for the
actual HyPhy-ready FASTA/tree used by FEL/MEME/BUSTED/aBSREL.

The dashboard builder computes Benjamini-Hochberg FDR q-values separately for
FEL and MEME. By default, selected sites are:

```text
FEL q <= 0.1
MEME q <= 0.1
```

Override these thresholds with:

```bash
python scripts/build_dashboard_data.py --fel-q 0.05 --meme-q 0.05
```

BUSTED summary output includes both `p_value` and `q_value`. In the current
single global BUSTED baseline, `q_value == p_value`; when multiple BUSTED
foreground tests are added, apply FDR across those BUSTED tests.

Run the local Streamlit dashboard:

```bash
bash scripts/run_dashboard.sh
```

On Windows PowerShell:

```powershell
powershell -ExecutionPolicy Bypass -File scripts/run_dashboard.ps1
```

On Windows Command Prompt:

```bat
scripts\run_dashboard.bat
```

If Streamlit is not installed:

```bash
python -m pip install -r requirements-dashboard.txt
```

Then open:

```text
http://localhost:8501
```

Dashboard pages:

```text
Overview
QC
Alignment
Selection Sites
Structure
Host/Metadata
Tree
Downloads
```

The Structure page uses the AlphaFold PDB model:

```text
results/colabfold/PP_006W0E7_1_97bc3/PP_006W0E7_1_97bc3_unrelaxed_rank_005_alphafold2_ptm_model_1_seed_000.pdb
```

Residues in this model are numbered 1-1138 and are annotated with the FEL/MEME
selection results by reference codon.

The dashboard tree page uses:

```text
results/hyphy_global_selection/ANDV_GPC_global.hyphy_ready.treefile
```

The host MRCA-labeled tree is optional and appears only after running
`scripts/label_iqtree2_tree_by_host_mrca.sh`.

## Recommended Command Block

For a full clean rerun:

```bash
bash scripts/recover_raw_cds_from_aligned_protein.sh
bash scripts/build_andv_gpc_reference_annotation.sh
bash scripts/generate_mafft_codon_alignment.sh
bash scripts/build_iqtree2_tree_from_mafft_codon_alignment.sh
bash scripts/run_global_hyphy_selection_scan.sh
bash scripts/run_relax_by_host_group.sh
python scripts/build_dashboard_data.py
```

Then, only if the global scan looks clean:

```bash
bash scripts/label_iqtree2_tree_by_host_mrca.sh
bash scripts/run_dashboard.sh
```
