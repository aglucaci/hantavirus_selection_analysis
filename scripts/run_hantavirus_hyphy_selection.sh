#!/usr/bin/env bash
set -euo pipefail

HYPHY_BIN="${HYPHY_BIN:-/opt/anaconda3/envs/hantavirus-snake/bin/hyphy}"
HYPHY_BF_DIR="${HYPHY_BF_DIR:-/opt/anaconda3/envs/hantavirus-snake/share/hyphy/TemplateBatchFiles/SelectionAnalyses}"
HYPHY_MSS_ANALYSIS="${HYPHY_MSS_ANALYSIS:-mss-ga}"
if [[ -n "${HYPHY_MSS_BF:-}" ]]; then
  HYPHY_MSS_ANALYSIS="$HYPHY_MSS_BF"
fi

HYPHY_CONFIG="${HYPHY_CONFIG:-config/hantavirus_hyphy_selection.env}"
if [[ -f "$HYPHY_CONFIG" ]]; then
  # shellcheck source=/dev/null
  source "$HYPHY_CONFIG"
fi

RESULTS_ROOT="${RESULTS_ROOT:-results}"
OUTDIR_OVERRIDE="${OUTDIR:-}"
LOG_DIR_OVERRIDE="${LOG_DIR:-}"
OUTDIR="${OUTDIR_OVERRIDE:-$RESULTS_ROOT/hantavirus_hyphy}"
INPUT_DIR="$OUTDIR/inputs"
MEME_IMPUTED_DIR="$OUTDIR/meme_imputed"
LOG_DIR="${LOG_DIR_OVERRIDE:-$RESULTS_ROOT/logs}"
SAMPLE_SHEET="${SAMPLE_SHEET:-sample_sheet.tsv}"
FORCE="${FORCE:-0}"
VERBOSE="${VERBOSE:-1}"
ONLY="${ONLY:-all}"
RUN_MSS="${RUN_MSS:-0}"
RELAX_STARTING_POINTS="${RELAX_STARTING_POINTS:-1}"
MSS_CLASSES="${MSS_CLASSES:-}"
SELECTION_BRANCHES="${SELECTION_BRANCHES:-Foreground}"
FEL_BRANCHES="${FEL_BRANCHES:-$SELECTION_BRANCHES}"
MEME_BRANCHES="${MEME_BRANCHES:-$SELECTION_BRANCHES}"
ABSREL_BRANCHES="${ABSREL_BRANCHES:-$SELECTION_BRANCHES}"
EXTRACT_ANCESTRAL="${EXTRACT_ANCESTRAL:-1}"
ANCESTRAL_METHODS="${ANCESTRAL_METHODS:-MEME_imputed}"
PROCESSED_OUTDIRS=""

print_command() {
  printf '[cmd]'
  printf ' %q' "$@"
  printf '\n'
}

run_cmd() {
  local name="$1"
  local output="$2"
  shift 2

  if [[ "$FORCE" != "1" && -s "$output" ]]; then
    echo "[skip] $name -> $output"
    return
  fi

  if [[ "$FORCE" == "1" && -e "$output" ]]; then
    echo "[force] removing previous output: $output"
    rm -f "$output"
  fi

  echo "[run] $name"
  echo "[out] $output"
  echo "[log] $LOG_DIR/${name}.log"
  print_command "$@"

  if [[ "$VERBOSE" == "1" ]]; then
    "$@" 2>&1 | tee "$LOG_DIR/${name}.log"
  else
    "$@" > "$LOG_DIR/${name}.log" 2>&1
  fi

  echo "[done] $name -> $output"
}

prepare_inputs() {
  local segment="$1"
  local group="$2"
  local alignment="$3"
  local tree="$4"
  local out_alignment="$INPUT_DIR/${segment}.hyphy_ready.fasta"
  local out_tree="$INPUT_DIR/${segment}_${group}.hyphy_ready.treefile"
  local qc="$INPUT_DIR/${segment}_${group}.hyphy_ready.qc.tsv"

  run_cmd "prepare_${segment}_${group}_hyphy_inputs" "$out_tree" \
    python scripts/filter_hyphy_codon_inputs.py \
      --alignment "$alignment" \
      --tree "$tree" \
      --out-alignment "$out_alignment" \
      --out-tree "$out_tree" \
      --qc "$qc" \
      --verbose
}

run_meme() {
  local segment="$1"
  local group="$2"
  local tree="$3"
  local output="$OUTDIR/${segment}_${group}_MEME.json"

  run_cmd "${segment}_${group}_MEME" "$output" \
    "$HYPHY_BIN" "$HYPHY_BF_DIR/MEME.bf" \
      --alignment "$INPUT_DIR/${segment}.hyphy_ready.fasta" \
      --tree "$tree" \
      --branches "$MEME_BRANCHES" \
      --output "$output"
}

run_meme_imputed() {
  local segment="$1"
  local group="$2"
  local tree="$3"
  local output="$MEME_IMPUTED_DIR/${segment}_${group}_MEME_imputed.json"

  mkdir -p "$MEME_IMPUTED_DIR"

  run_cmd "${segment}_${group}_MEME_imputed" "$output" \
    "$HYPHY_BIN" "$HYPHY_BF_DIR/MEME.bf" \
      --alignment "$INPUT_DIR/${segment}.hyphy_ready.fasta" \
      --tree "$tree" \
      --branches "$MEME_BRANCHES" \
      --impute-states Yes \
      --output "$output"
}

run_fel() {
  local segment="$1"
  local group="$2"
  local tree="$3"
  local output="$OUTDIR/${segment}_${group}_FEL.json"

  run_cmd "${segment}_${group}_FEL" "$output" \
    "$HYPHY_BIN" "$HYPHY_BF_DIR/FEL.bf" \
      --alignment "$INPUT_DIR/${segment}.hyphy_ready.fasta" \
      --tree "$tree" \
      --branches "$FEL_BRANCHES" \
      --output "$output"
}

run_absrel() {
  local segment="$1"
  local group="$2"
  local tree="$3"
  local output="$OUTDIR/${segment}_${group}_aBSREL.json"

  run_cmd "${segment}_${group}_aBSREL" "$output" \
    "$HYPHY_BIN" "$HYPHY_BF_DIR/aBSREL.bf" \
      --alignment "$INPUT_DIR/${segment}.hyphy_ready.fasta" \
      --tree "$tree" \
      --branches "$ABSREL_BRANCHES" \
      --output "$output"
}

run_relax() {
  local segment="$1"
  local group="$2"
  local tree="$3"
  local output="$OUTDIR/${segment}_${group}_RELAX.json"

  run_cmd "${segment}_${group}_RELAX" "$output" \
    "$HYPHY_BIN" "$HYPHY_BF_DIR/RELAX.bf" \
      --alignment "$INPUT_DIR/${segment}.hyphy_ready.fasta" \
      --tree "$tree" \
      --test Foreground \
      --reference "Unlabeled branches" \
      --starting-points "$RELAX_STARTING_POINTS" \
      --output "$output"
}

run_mss() {
  local segment="$1"
  local group="$2"
  local tree="$3"
  local output="$OUTDIR/${segment}_${group}_MSS.json"
  local nexus="$INPUT_DIR/${segment}_${group}.mss.nex"
  local filelist="$INPUT_DIR/${segment}_${group}.mss.files.txt"
  local mss_cmd=("$HYPHY_BIN" "$HYPHY_MSS_ANALYSIS")

  python scripts/build_mss_filelist.py \
    --alignment "$INPUT_DIR/${segment}.hyphy_ready.fasta" \
    --tree "$tree" \
    --out-nexus "$nexus" \
    --out-filelist "$filelist"

  if [[ -n "$MSS_CLASSES" ]]; then
    if [[ -f "$HYPHY_MSS_ANALYSIS" ]] && grep -q 'KeywordArgument ("classes"' "$HYPHY_MSS_ANALYSIS"; then
      mss_cmd+=(--classes "$MSS_CLASSES")
    else
      echo "[warn] $HYPHY_MSS_ANALYSIS does not expose --classes; ignoring MSS_CLASSES=$MSS_CLASSES"
    fi
  fi

  mss_cmd+=(--filelist "$filelist" --output "$output")
  run_cmd "${segment}_${group}_MSS" "$output" "${mss_cmd[@]}"
}

run_requested() {
  local method="$1"
  [[ "$ONLY" == "all" || "$ONLY" == "$method" ]]
}

record_output_dir() {
  local outdir="$1"
  case $'\n'"$PROCESSED_OUTDIRS"$'\n' in
    *$'\n'"$outdir"$'\n'*) ;;
    *) PROCESSED_OUTDIRS="${PROCESSED_OUTDIRS}"$'\n'"$outdir" ;;
  esac
}

extract_ancestral_states() {
  local outdir="$1"
  local output="$outdir/ancestral_sequences/ancestral_extraction_summary.tsv"
  local methods=()

  # shellcheck disable=SC2206
  methods=($ANCESTRAL_METHODS)

  if [[ "${#methods[@]}" -eq 0 ]]; then
    echo "[skip] ancestral state export; ANCESTRAL_METHODS is empty"
    return
  fi

  run_cmd "extract_ancestral_states_$(basename "$outdir")" "$output" \
    python scripts/extract_hyphy_ancestral_states.py \
      --results-dir "$outdir" \
      --methods "${methods[@]}"
}

configure_output_dirs() {
  local source_tag="$1"

  if [[ -n "$OUTDIR_OVERRIDE" ]]; then
    OUTDIR="$OUTDIR_OVERRIDE"
  elif [[ -n "$source_tag" ]]; then
    OUTDIR="$RESULTS_ROOT/$source_tag"
  else
    OUTDIR="$RESULTS_ROOT/hantavirus_hyphy"
  fi

  INPUT_DIR="$OUTDIR/inputs"
  MEME_IMPUTED_DIR="$OUTDIR/meme_imputed"

  if [[ -n "$LOG_DIR_OVERRIDE" ]]; then
    LOG_DIR="$LOG_DIR_OVERRIDE"
  elif [[ -n "$source_tag" ]]; then
    LOG_DIR="$OUTDIR/logs"
  else
    LOG_DIR="$RESULTS_ROOT/logs"
  fi

  mkdir -p "$OUTDIR" "$INPUT_DIR" "$MEME_IMPUTED_DIR" "$LOG_DIR"
}

if [[ ! -f "$SAMPLE_SHEET" ]]; then
  echo "ERROR: sample sheet not found: $SAMPLE_SHEET" >&2
  exit 1
fi

echo "[sample-sheet] $SAMPLE_SHEET"

line_number=0
has_source_tag=0
while IFS=$'\t' read -r segment group alignment tree source_tag extra || [[ -n "${segment:-}" ]]; do
  line_number=$((line_number + 1))

  if [[ "$line_number" -eq 1 ]]; then
    expected_header=$'segment\tgroup\talignment\ttree'
    expected_header_with_tag=$'segment\tgroup\talignment\ttree\tsource_tag'
    actual_header="${segment}"$'\t'"${group}"$'\t'"${alignment}"$'\t'"${tree}"
    actual_header_with_tag="$actual_header"$'\t'"${source_tag}"
    if [[ "$actual_header" == "$expected_header" && -z "${source_tag:-}" ]]; then
      has_source_tag=0
    elif [[ "$actual_header_with_tag" == "$expected_header_with_tag" ]]; then
      has_source_tag=1
    else
      echo "ERROR: sample sheet header must be either:" >&2
      echo "  $expected_header" >&2
      echo "  $expected_header_with_tag" >&2
      exit 1
    fi
    continue
  fi

  if [[ -z "${segment//[[:space:]]/}" || "${segment:0:1}" == "#" ]]; then
    continue
  fi

  if [[ "$has_source_tag" -eq 0 && -n "${source_tag:-}" ]]; then
    echo "ERROR: sample sheet line $line_number has too many columns." >&2
    exit 1
  fi

  if [[ -n "${extra:-}" ]]; then
    echo "ERROR: sample sheet line $line_number has too many columns." >&2
    exit 1
  fi

  if [[ -z "$segment" || -z "$group" || -z "$alignment" || -z "$tree" ]]; then
    echo "ERROR: sample sheet line $line_number has an empty required field." >&2
    exit 1
  fi

  if [[ "$has_source_tag" -eq 1 ]]; then
    if [[ -z "$source_tag" ]]; then
      echo "ERROR: sample sheet line $line_number has an empty source_tag field." >&2
      exit 1
    fi
    if [[ ! "$source_tag" =~ ^[A-Za-z0-9_.-]+$ ]]; then
      echo "ERROR: sample sheet line $line_number has an invalid source_tag: $source_tag" >&2
      exit 1
    fi
  else
    source_tag=""
  fi

  if [[ ! -f "$alignment" ]]; then
    echo "ERROR: alignment not found on sample sheet line $line_number: $alignment" >&2
    exit 1
  fi

  if [[ ! -f "$tree" ]]; then
    echo "ERROR: tree not found on sample sheet line $line_number: $tree" >&2
    exit 1
  fi

  configure_output_dirs "$source_tag"
  record_output_dir "$OUTDIR"
  echo "[results] $OUTDIR"

  prepare_inputs "$segment" "$group" "$alignment" "$tree"

  ready_tree="$INPUT_DIR/${segment}_${group}.hyphy_ready.treefile"
  if run_requested FEL; then run_fel "$segment" "$group" "$ready_tree"; fi
  if run_requested MEME; then
    run_meme "$segment" "$group" "$ready_tree"
    run_meme_imputed "$segment" "$group" "$ready_tree"
  fi
  if [[ "$ONLY" == "MEME_imputed" || "$ONLY" == "MEME_IMPUTED" ]]; then run_meme_imputed "$segment" "$group" "$ready_tree"; fi
  if run_requested aBSREL; then run_absrel "$segment" "$group" "$ready_tree"; fi
  if run_requested RELAX; then run_relax "$segment" "$group" "$ready_tree"; fi
  if run_requested MSS && [[ "$RUN_MSS" == "1" || "$ONLY" == "MSS" ]]; then run_mss "$segment" "$group" "$ready_tree"; fi
done < "$SAMPLE_SHEET"

if [[ "$EXTRACT_ANCESTRAL" == "1" ]]; then
  while IFS= read -r processed_outdir; do
    if [[ -n "$processed_outdir" ]]; then
      extract_ancestral_states "$processed_outdir"
    fi
  done <<< "$PROCESSED_OUTDIRS"
else
  echo "[skip] ancestral state export; EXTRACT_ANCESTRAL=$EXTRACT_ANCESTRAL"
fi

echo "[complete] HyPhy selection batch finished."
