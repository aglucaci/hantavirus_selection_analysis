Hantavirus ANDV M Segment / GPC CDS Recovery
============================================

For the current end-to-end command order, see:

```text
pipeline.md
```

Environment setup
-----------------

Create the Conda environment from the repository root:

```bash
conda env create -f environment.yml
```

Activate the environment:

```bash
conda activate hantavirus-snake
```

If the environment already exists, update it from the YAML instead:

```bash
conda env update -f environment.yml --prune
```

Confirm that Snakemake is available:

```bash
snakemake --version
```

Data source
-----------

Pathoplexus ANDV search:
https://pathoplexus.org/andv/search?orderBy=earliestReleaseDate&page=1&order=descending

Analysis date: 2026-05-11

Dataset contents
----------------

The working data are for the ANDV M segment / GPC region and include:

1. Aligned protein sequences
2. Raw nucleotide M-segment sequences
3. Metadata with all available columns

Current input files
-------------------

Aligned GPC protein FASTA:
data/andv_aligned-aa-GPC_2026-05-11T2017.fasta/andv_aligned-aa-GPC_2026-05-11T2017.fasta

Raw M-segment nucleotide FASTA:
data/andv_nuc-M_2026-05-11T1953.fasta/andv_nuc-M_2026-05-11T1953.fasta

Metadata:
data/andv_metadata_2026-05-11T1953.tsv/andv_metadata_2026-05-11T1953.tsv

Current workflow
----------------

Snakemake workflow
------------------

The workflow has been converted into a Snakemake pipeline driven by:

```text
Snakefile
config/config.yaml
```

The current configured inputs are the 2026-05-12 top-level files:

```text
data/andv_aligned-aa-GPC_2026-05-12T1516.fasta
data/andv_nuc-M_2026-05-12T1515.fasta
data/andv_metadata_2026-05-12T1515.tsv
```

By default, the Snakemake workflow filters the recovered CDS records to keep
only records with a known `hostNameScientific` value before MAFFT, IQ-TREE,
HyPhy, RELAX, and dashboard table generation. The full unfiltered recovery is
kept for audit here:

```text
results/GPC_CDS.recovered.raw.all_hosts.fasta
results/GPC_CDS.recovered.raw.qc.tsv
```

The filtered analysis FASTA and host-filter QC are:

```text
results/GPC_CDS.recovered.raw.fasta
results/GPC_CDS.recovered.raw.host_filter.qc.tsv
```

The configured coordinate reference is retained even though its May 12 metadata
row has a blank host name, because downstream coordinate annotation and HyPhy
input validation use it as the reference anchor. To make the host filter strict,
set `host_filter.keep_reference: false` in `config/config.yaml`.

Launch the Snakefile
--------------------

Run Snakemake from the repository root, the directory that contains
`Snakefile`:

```bash
cd /path/to/hantavirus_selection_analysis
```

Before running the workflow, install and activate the Conda environment as
described in "Environment setup" above.

Preview the workflow without running jobs:

```bash
snakemake -n --printshellcmds
```

Run the full workflow with 8 local cores:

```bash
snakemake --cores 8 --printshellcmds
```

Resume after an interrupted run with the same command. Snakemake will skip
outputs that are already complete and rebuild only missing or outdated targets.

On Windows, launch this from WSL or another shell where `bash`, `mafft`,
`iqtree2`, and `hyphy`/`hyphy-avx` are on `PATH`; the Snakefile explicitly uses
`bash` for its shell commands.

Important final targets include:

```text
results/hyphy_global_selection/ANDV_GPC_global_FEL.json
results/hyphy_global_selection/ANDV_GPC_global_MEME.json
results/hyphy_global_selection/ANDV_GPC_global_BUSTED.json
results/hyphy_global_selection/ANDV_GPC_global_aBSREL.json
results/hyphy_host_relax/ANDV_GPC_human_vs_reservoir_RELAX.json
results/dashboard/overview_metrics.json
results/dashboard/ANDV_GPC_selection_annotated_sites.tsv
```

The older shell launchers are still kept under `scripts/`, but the Snakefile is
now the preferred reproducible command graph.

The current active step recovers raw, ungapped coding sequences from the aligned
GPC protein FASTA and the raw nucleotide M-segment FASTA.

Script:
scripts/recover_raw_cds_from_aligned_protein.py

Launcher:
scripts/recover_raw_cds_from_aligned_protein.sh

Run from the repository root:

```bash
bash scripts/recover_raw_cds_from_aligned_protein.sh
```

Equivalent direct command:

```bash
python scripts/recover_raw_cds_from_aligned_protein.py \
  -p data/andv_aligned-aa-GPC_2026-05-11T2017.fasta/andv_aligned-aa-GPC_2026-05-11T2017.fasta \
  -n data/andv_nuc-M_2026-05-11T1953.fasta/andv_nuc-M_2026-05-11T1953.fasta \
  -o results/GPC_CDS.recovered.raw.fasta \
  --qc results/GPC_CDS.recovered.raw.qc.tsv \
  --min-identity 0.70 \
  --verbose
```

What the recovery script does
-----------------------------

For each aligned protein record, the script:

1. Removes protein alignment gaps.
2. Normalizes IDs so a nucleotide ID ending in "|M" matches the protein ID
   without that suffix.
3. Translates the raw nucleotide sequence in all six frames.
4. Finds the best contiguous translated region matching the ungapped protein.
5. Writes the corresponding raw, ungapped CDS sequence.
6. Writes a QC table for all records.

The script is standalone and does not require Biopython.

Use --verbose to print input paths, loaded record counts, per-record recovery
status, periodic progress summaries, and final QC counts.

Current outputs
---------------

Recovered raw CDS FASTA:
results/GPC_CDS.recovered.raw.fasta

QC report:
results/GPC_CDS.recovered.raw.qc.tsv

Codon-aligned CDS FASTA:
results/GPC_CDS.codon_aligned.fasta

Codon alignment QC report:
results/GPC_CDS.codon_aligned.qc.tsv

Most recent raw CDS recovery run:

```text
written: 127
no_alignment: 223
low_identity_skipped: 10
```

All written CDS records in the most recent run were multiples of 3.

Codon alignment workflow
------------------------

After raw CDS recovery, generate a codon alignment from the downloaded aligned
protein FASTA with:

```bash
bash scripts/generate_codon_alignment_from_protein.sh
```

Equivalent direct command:

```bash
python scripts/generate_codon_alignment_from_protein.py \
  -p data/andv_aligned-aa-GPC_2026-05-11T2017.fasta/andv_aligned-aa-GPC_2026-05-11T2017.fasta \
  -c results/GPC_CDS.recovered.raw.fasta \
  -o results/GPC_CDS.codon_aligned.fasta \
  --qc results/GPC_CDS.codon_aligned.qc.tsv \
  --min-translation-identity 0.70 \
  --missing-as gap \
  --verbose
```

The codon alignment script maps the aligned protein columns onto each raw CDS:

1. Protein gap columns become codon gaps: ---
2. Protein residue columns consume the next CDS codon.
3. Terminal protein stop symbols are skipped when the CDS has no terminal stop
   codon.
4. IDs are normalized the same way as in CDS recovery, so a trailing "|M" suffix
   is ignored.
5. A QC table reports missing CDS records, translation identity, used/unused
   codons, inserted missing codons, skipped terminal stops, and codon-alignment
   length.

Most recent codon alignment run:

```text
written: 127
missing_cds_record: 233
```

All written codon-aligned CDS records in the most recent run were length 3414 nt.

MAFFT codon alignment workflow
------------------------------

To create a fresh protein alignment from the recovered CDS using MAFFT, then
back-translate that MAFFT protein alignment to codons, run:

```bash
bash scripts/generate_mafft_codon_alignment.sh
```

This script:

1. Writes a sanitized raw-CDS FASTA with HyPhy-safe IDs.
2. Translates the sanitized CDS records to protein.
3. Aligns the translated proteins with MAFFT.
4. Back-translates the MAFFT protein alignment to codon-aligned CDS.

ID sanitization replaces characters such as "." and spaces with "_", for
example:

```text
PP_006VLSY.1 -> PP_006VLSY_1
```

Default outputs:

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

The full codon alignment is retained, but tree inference and downstream HyPhy
selection analyses use the exact-sequence deduplicated FASTA:

```text
results/GPC_CDS.mafft_codon_aligned.deduplicated.fasta
```

The duplicate map keeps track of every collapsed sequence. Current
deduplication result:

```text
input records: 127
unique sequences: 67
duplicate records removed: 60
```

The coordinate reference is never removed during deduplication. The MAFFT
launcher forces the reference ID to remain as the representative if it is in a
duplicate group:

```text
reference representative kept: PP_006W0E7_1
duplicate group size: 7
```

MAFFT must be available on PATH. One common install route is:

```bash
conda install -c bioconda mafft
```

IQ-TREE2 phylogeny workflow
---------------------------

After generating the MAFFT codon alignment and deduplicated FASTA, build a
phylogenetic tree with:

```bash
bash scripts/build_iqtree2_tree_from_mafft_codon_alignment.sh
```

Default input:

```text
results/GPC_CDS.mafft_codon_aligned.deduplicated.fasta
```

Default outputs:

```text
results/iqtree2/GPC_CDS.mafft_codon_aligned.deduplicated.treefile
results/iqtree2/GPC_CDS.mafft_codon_aligned.deduplicated.iqtree
results/iqtree2/GPC_CDS.mafft_codon_aligned.deduplicated.log
```

The script uses IQ-TREE2 model finding by default with ultrafast bootstrap and
SH-aLRT support:

```text
MODEL=MFP
BOOTSTRAPS=1000
ALRT=1000
THREADS=AUTO
```

You can override these settings inline, for example:

```bash
MODEL=GTR+G BOOTSTRAPS=2000 THREADS=8 bash scripts/build_iqtree2_tree_from_mafft_codon_alignment.sh
```

IQ-TREE2 must be available on PATH. One common install route is:

```bash
conda install -c bioconda iqtree
```

Tree host MRCA labeling workflow
--------------------------------

After building the IQ-TREE2 tree, label branches at the MRCA of each
hostNameScientific group with:

```bash
bash scripts/label_iqtree2_tree_by_host_mrca.sh
```

Default inputs:

```text
results/iqtree2/GPC_CDS.mafft_codon_aligned.deduplicated.treefile
data/andv_metadata_2026-05-11T1953.tsv/andv_metadata_2026-05-11T1953.tsv
```

Default outputs:

```text
results/iqtree2/GPC_CDS.mafft_codon_aligned.deduplicated.host_mrca_labeled.treefile
results/iqtree2/GPC_CDS.mafft_codon_aligned.deduplicated.host_mrca_labeled.qc.tsv
```

The script uses accessionVersion to match metadata rows to tree tips, then
groups tips by hostNameScientific. A trailing "|M" suffix is ignored during ID
matching. Labels are written as HyPhy-style Newick branch tags using the
sanitized hostNameScientific value by default:

```text
{Homo_sapiens}
```

The launcher skips Unknown_host by default and labels only groups with at least
two tree tips. The QC table reports how many tree tips were present, how many
metadata IDs were missing from the tree, and how many descendant tips are under
each labeled MRCA. Large MRCA descendant counts can indicate that a host group
is not monophyletic in the current tree. If multiple host groups map to the
same MRCA branch, their labels are joined with "+" in the Newick tag.
The labeler can match metadata accessions such as PP_006VLSY.1 to sanitized
tree tips such as PP_006VLSY_1.

To force a single branch label instead of host names:

```bash
BRANCH_LABEL=Foreground bash scripts/label_iqtree2_tree_by_host_mrca.sh
```

Global HyPhy selection scan
---------------------------

The first selection analysis should be global and unlabeled. Use the same
high-quality ANDV GPC codon alignment and the same tree, but do not use host
labels yet.

Run:

```bash
bash scripts/run_global_hyphy_selection_scan.sh
```

Default inputs:

```text
results/GPC_CDS.mafft_codon_aligned.deduplicated.fasta
results/iqtree2/GPC_CDS.mafft_codon_aligned.deduplicated.treefile
```

Default outputs:

```text
results/hyphy_global_selection/ANDV_GPC_global_FEL.json
results/hyphy_global_selection/ANDV_GPC_global_MEME.json
results/hyphy_global_selection/ANDV_GPC_global_BUSTED.json
results/hyphy_global_selection/ANDV_GPC_global_aBSREL.json
```

Before running HyPhy, the launcher prepares HyPhy-ready inputs by dropping
alignment sequences with in-frame stop codons or invalid codon lengths and
pruning the tree to the kept tips. It also renames sequence/tree tip IDs to
HyPhy-safe labels by replacing characters such as "." and spaces with "_".
FASTA headers are bare HyPhy-safe IDs only. Original-to-safe ID mapping is
preserved in the QC and duplicate-map tables:

```text
results/hyphy_global_selection/ANDV_GPC_global.hyphy_ready.fasta
results/hyphy_global_selection/ANDV_GPC_global.hyphy_ready.treefile
results/hyphy_global_selection/ANDV_GPC_global.hyphy_ready.qc.tsv
```

Current cleaning result:

```text
kept: 126
dropped: 1
dropped sequence: PP_006VXSM.1, stop codon at codon 1132 (TGA)
example rename: PP_006VLSY.1 -> PP_006VLSY_1
```

To bypass this cleaning step, run with CLEAN_INPUTS=0, but HyPhy will fail if
any sequence contains an in-frame stop codon.

This baseline scan asks:

```text
Is Andes virus GPC mostly under purifying selection, and are there specific
Gn/Gc codons or branches showing episodic positive selection?
```

The script runs:

```text
FEL    - pervasive purifying/diversifying selection at sites
MEME   - episodic diversifying selection at sites
BUSTED - gene-wide evidence of episodic positive selection
aBSREL - branches showing episodic selection
```

HyPhy must be available on PATH. One common install route is:

```bash
conda install -c bioconda hyphy
```

On some Linux/WSL installs, the `hyphy` wrapper mishandles arguments. The
launcher therefore prefers `hyphy-avx` or `/usr/lib/hyphy/bin/hyphy-avx` and
runs explicit HyPhy batch files. If your HyPhy executable is somewhere else,
override it with:

```bash
HYPHY_BIN=/path/to/hyphy-avx bash scripts/run_global_hyphy_selection_scan.sh
```

First report framing:

```text
We first performed a global codon-level selection scan of high-quality ANDV GPC
sequences across host sources. Because GPC mediates viral entry and contains
the Gn/Gc envelope glycoproteins, this scan tested whether the protein is
dominated by purifying selection while retaining evidence of episodic
diversifying selection at specific codons or branches.
```

Reference annotation workflow
-----------------------------

Build the first-pass ANDV GPC reference annotation with:

```bash
bash scripts/build_andv_gpc_reference_annotation.sh
```

This uses the local Pathoplexus RefSeq row:

```text
PP_006W0E7.1 == NC_003467.2.M
```

Reference sequence used for coordinate annotation:

```text
local Pathoplexus accessionVersion: PP_006W0E7.1
HyPhy/IQ-TREE-safe ID:              PP_006W0E7_1
M-segment RefSeq accession:         NC_003467.2
source database:                    RefSeq
specimen / strain label:            Chile-9717869
virus name in metadata:             Orthohantavirus andesense
M segment length:                   3671 nt
metadata completeness_M:            1.0
```

Default outputs:

```text
data/reference/ANDV_GPC_reference_annotation.tsv
data/reference/NC_003467.2_M.fasta
data/reference/NC_003467.2_GPC_CDS.no_stop.fasta
data/reference/NC_003467.2_GPC_protein.no_stop.fasta
```

Current NC_003467.2 / CHI-7913 coordinate layer:

```text
M segment length:       3671 nt
GPC CDS without stop:   nt 52-3465
terminal stop codon:    nt 3466-3468, TAA
GPC protein length:     1138 aa
WAASA motif:            aa 647-651
Gn first-pass region:   aa 1-651
Gc first-pass region:   aa 652-1138
```

The master annotation table columns are:

```text
feature
region_type
ref_accession
segment
protein
start_nt
end_nt
start_codon
end_codon
start_aa
end_aa
evidence
notes
```

Known first-pass rows include:

```text
GPC_CDS
GPC
Gn
Gc
WAASA_cleavage_motif
```

Rows for Gn/Gc ectodomains, transmembrane regions, and cytoplasmic tails are
included as pending placeholders until coordinates are filled from a topology
tool such as DeepTMHMM, TMHMM, or Phobius.

To annotate a selected-site table after extracting sites from HyPhy JSON:

```bash
python scripts/annotate_hyphy_sites_with_reference.py \
  -i results/hyphy_global_selection/selected_sites.tsv \
  -o results/hyphy_global_selection/selected_sites.annotated.tsv
```

The selected-site input should contain a site column named hyphy_site by
default, plus any result fields you want to preserve, for example:

```text
hyphy_site
test
p_value
q_value
omega
evidence_type
```

The annotator adds:

```text
reference_codon
reference_nt_start
reference_nt_end
aa_ref
region
```

Dashboard workflow
------------------

After the global HyPhy scan finishes, build dashboard-ready tables with:

```bash
python scripts/build_dashboard_data.py
```

Default dashboard data directory:

```text
results/dashboard/
```

The dashboard consumes:

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

The selection-input host summary reports host composition for the actual
HyPhy-ready FASTA/tree after deduplication and stop-codon filtering. It uses
the duplicate map so retained representatives can inherit host metadata from
collapsed identical sequences.

Selected-site calls in the dashboard are FDR-controlled with
Benjamini-Hochberg q-values computed separately for FEL and MEME. The default
threshold is:

```text
FEL q <= 0.1
MEME q <= 0.1
```

The all-sites FEL and MEME tables retain raw p-values and FDR q-values so plots
can show the full background while coloring only FDR-significant sites as
selected.

RELAX by host group
-------------------

After the global baseline, test selection intensity in human-derived versus
reservoir-derived branches:

```bash
bash scripts/run_relax_by_host_group.sh
```

Default comparison:

```text
Test:      Homo sapiens terminal branches
Reference: known non-human host terminal branches
Skipped:   unknown/missing host branches
```

The RELAX labeling step uses the deduplication map so retained representative
tips can inherit host evidence from identical collapsed sequences. Mixed
human/non-human duplicate groups are left unlabeled to avoid ambiguous host
assignment.

Outputs:

```text
results/hyphy_host_relax/ANDV_GPC_human_vs_reservoir_RELAX.treefile
results/hyphy_host_relax/ANDV_GPC_human_vs_reservoir_RELAX.qc.tsv
results/hyphy_host_relax/ANDV_GPC_human_vs_reservoir_RELAX.json
```

After RELAX finishes, rebuild dashboard data:

```bash
python scripts/build_dashboard_data.py
```

BUSTED is reported with both a raw p-value and an FDR q-value. For the current
single global BUSTED baseline test, the q-value equals the p-value. If multiple
BUSTED foreground tests are added later, apply FDR across those BUSTED tests.

The dashboard tree page uses the HyPhy-ready tree by default:

```text
results/hyphy_global_selection/ANDV_GPC_global.hyphy_ready.treefile
```

This is the tree actually used for the global FEL/MEME/BUSTED/aBSREL scan. The
host MRCA-labeled tree is optional and is shown only if it has been generated.

Install dashboard dependencies if needed:

```bash
python -m pip install -r requirements-dashboard.txt
```

Start the local dashboard:

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

Structure page input:

```text
results/colabfold/PP_006W0E7_1_97bc3/PP_006W0E7_1_97bc3_unrelaxed_rank_005_alphafold2_ptm_model_1_seed_000.pdb
```

This AlphaFold model is residue-numbered 1-1138, matching the NC_003467.2 /
CHI-7913 GPC reference codon coordinate system used by the selection plots.

The main dashboard table is:

```text
results/dashboard/ANDV_GPC_selection_annotated_sites.tsv
```

It maps HyPhy-selected sites onto the NC_003467.2 / CHI-7913 coordinate layer,
including reference codon, nucleotide coordinates, reference amino acid, and
first-pass GPC region annotation.

QC status meanings
------------------

written:
The sequence passed the identity threshold and was written to the output FASTA.

no_alignment:
No valid translated nucleotide window was found. In the current data this often
corresponds to nucleotide records that are much shorter than the full GPC
protein.

low_identity_skipped:
A candidate coding region was found, but its protein/translation identity was
below the selected --min-identity cutoff.

missing_nucleotide_record:
The protein record did not have a matching nucleotide record after ID
normalization.

empty_sequence:
The protein or nucleotide sequence was empty after cleaning.

Notes
-----

Archived exploratory scripts and earlier workflow attempts are under:
scripts/archive/

Author:
Alexander G. Lucaci
