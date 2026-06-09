# Hantavirus HyPhy Selection Analyses

This repository runs HyPhy FEL, MEME, MEME with imputed states, aBSREL, and
RELAX on prepared hantavirus codon alignments and labeled trees. MSS-GA is
available but skipped by default.

## Environment

Create and activate the Conda environment:

```bash
conda env create -f environment.yml
conda activate hantavirus-snake
```

If the environment already exists:

```bash
conda env update -f environment.yml --prune
conda activate hantavirus-snake
```

Confirm HyPhy is available:

```bash
hyphy --version
```

## Input Data

The workflow reads analyses from:

```text
sample_sheet.tsv
```

For the current ANDV tree/alignment analyses, use:

```text
andv_trees_sample_sheet.tsv
```

The sample sheet has one row per segment and foreground-label group. The
`source_tag` field controls the default results folder under `results/`:

```text
segment	group	alignment	tree	source_tag
L	ANDV	data/hantavirus/sequences_L_ex_lab_cells_2026_05_27v7.cds.linsi.fasta	data/hantavirus/sequences_L_ex_lab_cells_2026_05_27v7.cds.linsi.rt.annot-ANDV.nwk	hantavirus_hyphy
L	clade3	data/hantavirus/sequences_L_ex_lab_cells_2026_05_27v7.cds.linsi.fasta	data/hantavirus/sequences_L_ex_lab_cells_2026_05_27v7.cds.linsi.rt.annot-clade3.nwk	hantavirus_hyphy
S	ANDV	data/hantavirus/sequences_S_ex_lab_cells_2026_05_27v5.cds.linsi.fasta	data/hantavirus/sequences_S_ex_lab_cells_2026_05_27v5.cds.linsi.rt.annot-ANDV.nwk	hantavirus_hyphy
S	clade3	data/hantavirus/sequences_S_ex_lab_cells_2026_05_27v5.cds.linsi.fasta	data/hantavirus/sequences_S_ex_lab_cells_2026_05_27v5.cds.linsi.rt.annot-clade3.nwk	hantavirus_hyphy
```

To run with another sheet:

```bash
SAMPLE_SHEET=path/to/sample_sheet.tsv bash scripts/run_hantavirus_hyphy_selection.sh
```

To regenerate the repository sample sheets:

```bash
python scripts/generate_sample_sheets.py
```

The current `andv_trees_sample_sheet.tsv` analyses use these prepared input
files:

```text
data/ANDV_trees_aln-hyphy/sequences_L_ex_lab_cells_2026_05_27v7.cds.linsi.fasta
data/ANDV_trees_aln-hyphy/sequences_L_ex_lab_cells_2026_05_27v7.cds.linsi.rt.hyphy-ANDVall.nwk
data/ANDV_trees_aln-hyphy/sequences_L_ex_lab_cells_2026_05_27v7.cds.linsi.rt.hyphy-ANDVstems.nwk
data/ANDV_trees_aln-hyphy/sequences_L_ex_lab_cells_2026_05_27v7.cds.linsi.rt.hyphy-humanout.nwk
data/ANDV_trees_aln-hyphy/sequences_M_ex_lab_cells_concat_2026_05_27v6.cds.linsi.fasta
data/ANDV_trees_aln-hyphy/sequences_M_ex_lab_cells_concat_2026_05_27v6.cds.linsi.rt.hyphy-ANDVall.nwk
data/ANDV_trees_aln-hyphy/sequences_M_ex_lab_cells_concat_2026_05_27v6.cds.linsi.rt.hyphy-ANDVstems.nwk
data/ANDV_trees_aln-hyphy/sequences_M_ex_lab_cells_concat_2026_05_27v6.cds.linsi.rt.hyphy-humanout.nwk
data/ANDV_trees_aln-hyphy/sequences_S_ex_lab_cells_2026_05_27v5.cds.linsi.fasta
data/ANDV_trees_aln-hyphy/sequences_S_ex_lab_cells_2026_05_27v5.cds.linsi.rt.hyphy-ANDVall.nwk
data/ANDV_trees_aln-hyphy/sequences_S_ex_lab_cells_2026_05_27v5.cds.linsi.rt.hyphy-ANDVstems.nwk
data/ANDV_trees_aln-hyphy/sequences_S_ex_lab_cells_2026_05_27v5.cds.linsi.rt.hyphy-humanout.nwk
```

The tree files use HyPhy `{Foreground}` labels. RELAX tests those foreground
branches against unlabeled branches, and FEL, MEME, and aBSREL test the labeled
foreground branches by default.

## Run

From the repository root:

```bash
bash scripts/run_hantavirus_hyphy_selection.sh
```

To rerun every default step and stream verbose HyPhy output to the terminal:

```bash
FORCE=1 VERBOSE=1 bash scripts/run_hantavirus_hyphy_selection.sh
```

For the current ANDV tree/alignment sample sheet:

```bash
SAMPLE_SHEET=andv_trees_sample_sheet.tsv FORCE=1 VERBOSE=1 bash scripts/run_hantavirus_hyphy_selection.sh
```

The default run includes FEL, two MEME runs, aBSREL, and RELAX. The standard
MEME run writes `*_MEME.json` and is the only MEME output used for selection
q-value summaries. The second MEME run uses HyPhy `--impute-states Yes` and
writes `meme_imputed/*_MEME_imputed.json` plus matching
`*_MEME_imputed.log` files for ancestral/imputed-state export only. It does not
run MSS-GA unless you explicitly opt in with `RUN_MSS=1` or `ONLY=MSS`.

Important: `meme_imputed/*_MEME_imputed.json` files are not used for MEME
q-value calculations, dashboard significance counts, or evidence-tier selection
summaries. They are stored in a separate folder and generated only to provide
imputed/reconstructed states for the `ancestral_sequences/` export.

After HyPhy finishes, the workflow extracts full codon and amino-acid FASTA
sequences from the MEME imputed-state JSON files into `ancestral_sequences/`.
Disable this with `EXTRACT_ANCESTRAL=0`, or include standard MEME and other
methods with substitution maps using
`ANCESTRAL_METHODS="FEL MEME MEME_imputed aBSREL RELAX"`.

To rerun only RELAX with more random starting points:

```bash
ONLY=RELAX FORCE=1 VERBOSE=1 RELAX_STARTING_POINTS=10 bash scripts/run_hantavirus_hyphy_selection.sh
```

`RELAX_STARTING_POINTS` is passed to HyPhy as `--starting-points`. Increasing
it can help check whether RELAX convergence warnings or local maxima are
affecting the inferred K value, but it will make RELAX slower.

You can also set this in the workflow config file:

```text
config/hantavirus_hyphy_selection.env
```

For example:

```bash
: "${RELAX_STARTING_POINTS:=10}"
```

The config file is sourced automatically by
`scripts/run_hantavirus_hyphy_selection.sh`. Inline environment variables still
override config values, so this still works:

```bash
RELAX_STARTING_POINTS=25 ONLY=RELAX FORCE=1 bash scripts/run_hantavirus_hyphy_selection.sh
```

To run only FEL for the ANDV tree/alignment sample sheet:

```bash
ONLY=FEL VERBOSE=0 SAMPLE_SHEET=andv_trees_sample_sheet.tsv bash scripts/run_hantavirus_hyphy_selection.sh
```

To run only the MEME imputed-state analysis for the ANDV tree/alignment sample
sheet:

```bash
ONLY=MEME_IMPUTED FORCE=1 VERBOSE=1 SAMPLE_SHEET=andv_trees_sample_sheet.tsv bash scripts/run_hantavirus_hyphy_selection.sh
```

FEL, MEME, and aBSREL default to HyPhy's `Foreground` branch selector, i.e.
branches labeled `{Foreground}` in the input tree. To reproduce whole-tree
branch-tested runs, set `SELECTION_BRANCHES=All`. Individual methods can be
overridden with `FEL_BRANCHES`, `MEME_BRANCHES`, or `ABSREL_BRANCHES`.

To run only MSS-GA for the ANDV tree/alignment sample sheet:

```bash
ONLY=MSS VERBOSE=0 SAMPLE_SHEET=andv_trees_sample_sheet.tsv bash scripts/run_hantavirus_hyphy_selection.sh
```

To include MSS-GA in a full run:

```bash
RUN_MSS=1 FORCE=1 VERBOSE=1 SAMPLE_SHEET=andv_trees_sample_sheet.tsv bash scripts/run_hantavirus_hyphy_selection.sh
```

The MSS command generated by the workflow follows this shape:

```bash
hyphy mss-ga --filelist files.txt --output outfile.json
```

This uses HyPhy's MSS-GA command-line analysis, the genetic-algorithm search
for synonymous codon-selection models. If you need to run a specific batch file
instead, set `HYPHY_MSS_BF=/path/to/MSS-selector.bf`. `MSS_CLASSES` is only
passed when the selected batch file exposes a `--classes` option; the installed
`mss-ga` alias does not.

To run in the background:

```bash
FORCE=1 VERBOSE=1 nohup bash scripts/run_hantavirus_hyphy_selection.sh > results/logs/hantavirus_hyphy_batch.log 2>&1 &
echo $! > results/logs/hantavirus_hyphy_batch.pid
```

Monitor progress:

```bash
tail -f results/logs/hantavirus_hyphy_batch.log
```

Before rerunning after a failed or interrupted job, remove zero-byte JSON files:

```bash
find results/hantavirus_hyphy -name '*.json' -size 0 -delete
```

## Outputs

Cleaned HyPhy-ready inputs are written to:

```text
results/hantavirus_hyphy/inputs/
```

For `andv_trees_sample_sheet.tsv`, cleaned inputs are written to:

```text
results/ANDV_trees_aln-hyphy/inputs/
```

HyPhy JSON outputs are written to:

```text
results/hantavirus_hyphy/
```

For `andv_trees_sample_sheet.tsv`, outputs are written to:

```text
results/ANDV_trees_aln-hyphy/
```

MEME imputed-state JSON files are kept separate from the selection-test JSONs:

```text
results/ANDV_trees_aln-hyphy/meme_imputed/*_MEME_imputed.json
```

Per-step logs are written to:

```text
results/logs/
```

For `andv_trees_sample_sheet.tsv`, per-step logs are written to:

```text
results/ANDV_trees_aln-hyphy/logs/
```

## JavaScript Dashboard

### Run Locally

From the repository root:

```bash
cd /Users/agl4001/Documents/hantavirus_selection_analysis
conda activate hantavirus-snake
```

Build the dashboard summary tables for the current ANDV results:

```bash
bash scripts/build_dashboard.sh
```

The builder reads HyPhy JSON outputs from `results/ANDV_trees_aln-hyphy/` by
default and writes normalized TSV files to:

```text
results/ANDV_trees_aln-hyphy/dashboard_tables/
```

To build tables for a different results folder:

```bash
HYPHY_RESULTS_DIR=results/hantavirus_hyphy bash scripts/build_dashboard.sh
```

Start the local JavaScript dashboard server:

```bash
python -m http.server 8502 --bind 127.0.0.1
```

Then open:

```text
http://127.0.0.1:8502/dashboard-js/
```

Shortcut: build the tables and start the local server in one command:

```bash
bash scripts/run_js_dashboard.sh
```

The shortcut is portable across macOS and Ubuntu/WSL. It changes to the
repository root, finds `python3`, `python`, or the `hantavirus-snake` conda
environment, checks the dashboard files, rebuilds the normalized TSV tables,
and then starts a local static server. The dashboard builder uses the Python
standard library only, so no npm install is required for the JavaScript
dashboard.

The shortcut starts at port `8502` and automatically tries the next port if
that one is already in use. On WSL 2, open the printed `127.0.0.1` URL in your
Windows or Ubuntu browser, for example:

```text
http://127.0.0.1:8502/dashboard-js/
```

To choose a different starting port:

```bash
DASHBOARD_PORT=8510 bash scripts/run_js_dashboard.sh
```

For a different results folder with the shortcut:

```bash
HYPHY_RESULTS_DIR=results/hantavirus_hyphy bash scripts/run_js_dashboard.sh
```

If Ubuntu/WSL does not have Python available, install it or create the conda
environment first:

```bash
sudo apt-get update
sudo apt-get install -y python3
# or
conda env create -f environment.yml
```

To rebuild the MSS manuscript tables and relaunch the JavaScript dashboard:

```bash
python scripts/build_mss_interpretation_tables.py
bash scripts/run_js_dashboard.sh
```

This dashboard is a static JavaScript app that reads the normalized TSV tables
from `results/ANDV_trees_aln-hyphy/dashboard_tables/`. It summarizes the methods
present in this run: FEL, MEME, aBSREL, RELAX, and any opt-in MSS-GA outputs.
BUSTED, CFEL, and GARD are shown as not run until those outputs are added to
the workflow.

Each displayed table can be exported as TSV, and each chart can be exported as
a publication-ready SVG from the dashboard controls.

The gene browser links site-level and branch-level evidence: foreground FEL and
foreground MEME sites are shown on multi-track codon maps with FDR q-value
thresholds, tested-codon rug marks, selected-site labels, and an
alignment-quality track; branch views show RELAX Test branches and foreground
aBSREL-selected branches on the prepared labeled tree, plus RELAX
effect/reliability plots and foreground aBSREL branch-evidence plots.

The dashboard includes an **MSS Charts** tab for MSS-GA plots and tables. The
MSS Context table is sortable by column and can be exported as TSV.

Foreground FEL and foreground MEME site calls are Benjamini-Hochberg adjusted
within each method/run.
Dashboard site counts and the site-threshold slider use `q <= 0.10` by default;
`q <= 0.05` is the stricter support tier. Raw `p <= 0.10` rows are retained in
the tables as exploratory context.

The MEME tab also includes a branch EBF table reconstructed from MEME's stored
branch posterior annotations. It reports branch, codon, posterior positive-class
support, reconstructed EBF, FDR status, and reconstructed codon state so the
branches contributing to each MEME site can be inspected directly.

## MSS Interpretation Tables

Build manuscript-ready MSS-GA interpretation tables:

```bash
python scripts/build_mss_interpretation_tables.py
```

Outputs are written to:

```text
results/ANDV_trees_aln-hyphy/mss_interpretation_tables/
```

These include figure-ready tables for the MSS support heatmap, median active
parameter bars, top-parameter recurrence dot matrix, plain-language callouts,
radar summary, and figure legend text.

## Streamlit Dashboard

Build the dashboard summary tables:

```bash
bash scripts/build_dashboard.sh
```

Open the interactive Streamlit dashboard:

```bash
streamlit run dashboard/selection_dashboard.py
```

Or run both steps with:

```bash
bash scripts/run_hyphy_dashboard.sh
```

Then open:

```text
http://127.0.0.1:8501
```

The dashboard currently summarizes the methods present in this run: FEL, MEME,
aBSREL, RELAX, and any opt-in MSS-GA outputs. BUSTED, CFEL, and GARD are shown
as not run until those outputs are added to the workflow.

Dashboard-ready TSV tables are written to:

```text
results/ANDV_trees_aln-hyphy/dashboard_tables/
```

## Ancestral State Exports

HyPhy MEME with `--impute-states Yes` stores site-wise imputed codon
probabilities in `MLE -> Imputed States`. The extractor converts those
probabilities into full-length FASTA files by taking the highest-probability
codon at each site. The prepared HyPhy-ready alignment is used only as a
fallback when HyPhy emits no imputed state for a site.

Extract full MEME imputed-state codon and amino-acid FASTA files:

```bash
python scripts/extract_hyphy_ancestral_states.py --results-dir results/ANDV_trees_aln-hyphy --methods MEME_imputed
```

The `MEME_imputed` outputs in this folder should be interpreted as
ancestral/imputed-state exports only. Use the standard `*_MEME.json` outputs for
MEME site q-values and selection calls.

The source JSON files for this export are stored separately from selection-test
JSONs:

```text
results/ANDV_trees_aln-hyphy/meme_imputed/*_MEME_imputed.json
```

The default output files are full FASTA files, not sparse codon maps:

```text
results/ANDV_trees_aln-hyphy/ancestral_sequences/*_MEME_imputed.ancestral_codons.fasta
results/ANDV_trees_aln-hyphy/ancestral_sequences/*_MEME_imputed.ancestral_amino_acids.fasta
```

These FASTA files include terminal taxa, the `root`, and HyPhy internal
`Node*` records. The matching trees are exported into the same folder:

```text
results/ANDV_trees_aln-hyphy/ancestral_sequences/*_MEME_imputed.hyphy_ready.treefile
results/ANDV_trees_aln-hyphy/ancestral_sequences/*_MEME_imputed.hyphy_node_labeled.nwk
```

Use the `*.hyphy_node_labeled.nwk` files when tracking mutations by HyPhy
internal node IDs. The long TSV files include `node_type` and `state_source`
columns so you can distinguish root states, explicit branch substitutions,
inherited states, and imputed tip states.

If you explicitly include FEL, MEME, aBSREL, or RELAX in `--methods`, the
extractor can still write sparse substitution-map exports for those methods.

The main HyPhy runner performs that extraction automatically by default after
analyses finish:

```bash
SAMPLE_SHEET=andv_trees_sample_sheet.tsv FORCE=1 VERBOSE=1 bash scripts/run_hantavirus_hyphy_selection.sh
```

To skip automatic ancestral-state extraction:

```bash
EXTRACT_ANCESTRAL=0 SAMPLE_SHEET=andv_trees_sample_sheet.tsv FORCE=1 VERBOSE=1 bash scripts/run_hantavirus_hyphy_selection.sh
```

Extract all current methods with substitution maps:

```bash
python scripts/extract_hyphy_ancestral_states.py --results-dir results/ANDV_trees_aln-hyphy --methods FEL MEME MEME_imputed aBSREL RELAX
```

Outputs are written to:

```text
results/ANDV_trees_aln-hyphy/ancestral_sequences/
```

## Active Scripts

Only these scripts are part of the current workflow:

```text
scripts/filter_hyphy_codon_inputs.py
scripts/run_hantavirus_hyphy_selection.sh
scripts/build_mss_filelist.py
scripts/build_hyphy_dashboard_tables.py
scripts/build_dashboard.sh
scripts/extract_hyphy_ancestral_states.py
scripts/build_mss_interpretation_tables.py
scripts/run_js_dashboard.sh
scripts/run_hyphy_dashboard.sh
```

Old ANDV GPC recovery, Snakemake, dashboard, and reference-annotation scripts
were moved to:

```text
archive/old_andv_gpc_scripts/
```
