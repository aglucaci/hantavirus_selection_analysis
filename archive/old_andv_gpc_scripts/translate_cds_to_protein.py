#!/usr/bin/env python3
"""Translate CDS FASTA records to protein FASTA records."""

import argparse
import re


CODON_TABLE = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "TGT": "C",
    "TGC": "C",
    "TGA": "*",
    "TGG": "W",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
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


def translate(seq):
    seq = re.sub(r"[^ACGTN]", "", seq.upper().replace("U", "T"))
    seq = seq[: len(seq) - (len(seq) % 3)]
    aa = "".join(CODON_TABLE.get(seq[i : i + 3], "X") for i in range(0, len(seq), 3))
    if aa.endswith("*"):
        aa = aa[:-1]
    return aa.replace("*", "X")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--output", required=True)
    args = parser.parse_args()

    records = [(seq_id, translate(seq)) for seq_id, seq in read_fasta(args.input)]
    records = [(seq_id, aa) for seq_id, aa in records if aa]
    if not records:
        raise SystemExit("ERROR: no protein sequences written.")

    with open(args.output, "w") as handle:
        for seq_id, aa in records:
            handle.write(f">{seq_id} translated_from_recovered_CDS\n")
            for i in range(0, len(aa), 70):
                handle.write(aa[i : i + 70] + "\n")

    print(f"Saved translated protein FASTA: {args.output}")
    print(f"Sequences written: {len(records)}")


if __name__ == "__main__":
    main()
