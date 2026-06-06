cd /Users/agl4001/Documents/hantavirus_selection_analysis
conda activate hantavirus-snake

SAMPLE_SHEET=andv_trees_sample_sheet.tsv \
ONLY=MSS \
FORCE=1 \
VERBOSE=1 \
bash scripts/run_hantavirus_hyphy_selection.sh
