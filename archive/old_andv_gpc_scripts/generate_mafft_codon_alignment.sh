#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Generate a codon alignment using MAFFT
#
# Steps:
#   1. Sanitize recovered raw CDS IDs.
#   2. Translate recovered raw CDS to protein.
#   3. Align translated proteins with MAFFT.
#   4. Back-translate the MAFFT protein alignment to codon-aligned CDS.
#   5. Deduplicate the codon alignment for tree/HyPhy downstream use.
###############################################################################

RAW_CDS="${RAW_CDS:-results/GPC_CDS.recovered.raw.fasta}"
RAW_CDS_SANITIZED="${RAW_CDS_SANITIZED:-results/GPC_CDS.recovered.raw.hyphy_ids.fasta}"
ID_MAP="${ID_MAP:-results/GPC_CDS.recovered.raw.hyphy_ids.map.tsv}"
PROTEIN_RAW="${PROTEIN_RAW:-results/GPC_protein.from_recovered_CDS.fasta}"
PROTEIN_ALIGNED="${PROTEIN_ALIGNED:-results/GPC_protein.from_recovered_CDS.mafft.aligned.fasta}"
CODON_ALIGNED="${CODON_ALIGNED:-results/GPC_CDS.mafft_codon_aligned.fasta}"
CODON_QC="${CODON_QC:-results/GPC_CDS.mafft_codon_aligned.qc.tsv}"
CODON_ALIGNED_DEDUP="${CODON_ALIGNED_DEDUP:-results/GPC_CDS.mafft_codon_aligned.deduplicated.fasta}"
DEDUP_MAP="${DEDUP_MAP:-results/GPC_CDS.mafft_codon_aligned.duplicates.tsv}"
REFERENCE_ID="${REFERENCE_ID:-PP_006W0E7_1}"
THREADS="${THREADS:--1}"
MIN_TRANSLATION_IDENTITY="${MIN_TRANSLATION_IDENTITY:-0.70}"

echo "[setup] Checking inputs and tools..."

if [[ ! -s "$RAW_CDS" ]]; then
  echo "ERROR: raw CDS FASTA not found or empty: $RAW_CDS"
  echo "Run scripts/recover_raw_cds_from_aligned_protein.sh first."
  exit 1
fi

if [[ ! -s scripts/generate_codon_alignment_from_protein.py ]]; then
  echo "ERROR: missing script: scripts/generate_codon_alignment_from_protein.py"
  exit 1
fi

if [[ ! -s scripts/deduplicate_fasta_by_sequence.py ]]; then
  echo "ERROR: missing script: scripts/deduplicate_fasta_by_sequence.py"
  exit 1
fi

if ! command -v mafft >/dev/null 2>&1; then
  echo "ERROR: MAFFT not found in PATH."
  echo "Install with: conda install -c bioconda mafft"
  exit 1
fi

mkdir -p "$(dirname "$RAW_CDS_SANITIZED")" "$(dirname "$PROTEIN_RAW")" "$(dirname "$PROTEIN_ALIGNED")" "$(dirname "$CODON_ALIGNED")" "$(dirname "$CODON_ALIGNED_DEDUP")"

echo "[1/4] Sanitizing raw CDS sequence IDs for MAFFT/IQ-TREE/HyPhy..."

python - "$RAW_CDS" "$RAW_CDS_SANITIZED" "$ID_MAP" <<'PY'
import re
import sys

inp, out_fasta, out_map = sys.argv[1], sys.argv[2], sys.argv[3]

def hyphy_safe_id(seq_id):
    safe = re.sub(r"[^A-Za-z0-9_]+", "_", seq_id.strip())
    safe = re.sub(r"_+", "_", safe).strip("_")
    if not safe:
        safe = "sequence"
    if safe[0].isdigit():
        safe = "seq_" + safe
    return safe

def read_fasta(path):
    records = []
    current_id = None
    current_desc = ""
    current_seq = []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    records.append((current_id, current_desc, "".join(current_seq)))
                current_desc = line[1:].strip()
                current_id = current_desc.split()[0]
                current_seq = []
            else:
                current_seq.append(line)
    if current_id is not None:
        records.append((current_id, current_desc, "".join(current_seq)))
    return records

records = read_fasta(inp)
safe_seen = {}

with open(out_fasta, "w") as fasta, open(out_map, "w") as id_map:
    id_map.write("original_id\thyphy_id\toriginal_description\n")
    for original_id, original_desc, seq in records:
        safe_id = hyphy_safe_id(original_id)
        previous = safe_seen.get(safe_id)
        if previous is not None and previous != original_id:
            raise SystemExit(
                f"ERROR: sanitized ID collision: {previous!r} and {original_id!r} both become {safe_id!r}"
            )
        safe_seen[safe_id] = original_id
        id_map.write(f"{original_id}\t{safe_id}\t{original_desc}\n")
        fasta.write(f">{safe_id}\n")
        seq = re.sub(r"\s+", "", seq)
        for i in range(0, len(seq), 70):
            fasta.write(seq[i:i + 70] + "\n")

print(f"Saved sanitized raw CDS FASTA: {out_fasta}")
print(f"Saved ID map: {out_map}")
print(f"Sequences written: {len(records)}")
PY

echo "[2/4] Translating sanitized raw CDS to protein..."

python - "$RAW_CDS_SANITIZED" "$PROTEIN_RAW" <<'PY'
import re
import sys

inp, out = sys.argv[1], sys.argv[2]

codon_table = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

def read_fasta(path):
    records = []
    current_id = None
    current_seq = []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    records.append((current_id, "".join(current_seq)))
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
    if current_id is not None:
        records.append((current_id, "".join(current_seq)))
    return records

def clean_nt(seq):
    return re.sub(r"[^ACGTN]", "", seq.upper().replace("U", "T"))

def translate(seq):
    seq = clean_nt(seq)
    seq = seq[: len(seq) - (len(seq) % 3)]
    aa = []
    for i in range(0, len(seq), 3):
        aa.append(codon_table.get(seq[i:i + 3], "X"))
    aa = "".join(aa)
    if aa.endswith("*"):
        aa = aa[:-1]
    return aa.replace("*", "X")

records = []
for seq_id, seq in read_fasta(inp):
    aa = translate(seq)
    if aa:
        records.append((seq_id, aa))

with open(out, "w") as handle:
    for seq_id, aa in records:
        handle.write(f">{seq_id} translated_from_recovered_CDS\n")
        for i in range(0, len(aa), 70):
            handle.write(aa[i:i + 70] + "\n")

print(f"Saved translated protein FASTA: {out}")
print(f"Sequences written: {len(records)}")

if not records:
    raise SystemExit("ERROR: no protein sequences written.")
PY

echo "[3/4] Aligning translated proteins with MAFFT..."

mafft \
  --auto \
  --thread "$THREADS" \
  "$PROTEIN_RAW" > "$PROTEIN_ALIGNED"

echo "[4/5] Back-translating MAFFT protein alignment to codon alignment..."

python scripts/generate_codon_alignment_from_protein.py \
  -p "$PROTEIN_ALIGNED" \
  -c "$RAW_CDS_SANITIZED" \
  -o "$CODON_ALIGNED" \
  --qc "$CODON_QC" \
  --min-translation-identity "$MIN_TRANSLATION_IDENTITY" \
  --missing-as gap \
  --verbose

echo "[5/5] Deduplicating codon alignment for tree and HyPhy analyses..."

python scripts/deduplicate_fasta_by_sequence.py \
  -i "$CODON_ALIGNED" \
  -o "$CODON_ALIGNED_DEDUP" \
  --map "$DEDUP_MAP" \
  --keep-id "$REFERENCE_ID" \
  --verbose

echo
echo "[done]"
echo "Raw CDS:                   $RAW_CDS"
echo "Sanitized raw CDS:         $RAW_CDS_SANITIZED"
echo "ID map:                    $ID_MAP"
echo "Translated protein FASTA:  $PROTEIN_RAW"
echo "MAFFT protein alignment:   $PROTEIN_ALIGNED"
echo "Codon alignment:           $CODON_ALIGNED"
echo "Codon alignment QC:        $CODON_QC"
echo "Deduplicated codon aln:    $CODON_ALIGNED_DEDUP"
echo "Duplicate map:             $DEDUP_MAP"
echo "Reference kept as rep:     $REFERENCE_ID"
