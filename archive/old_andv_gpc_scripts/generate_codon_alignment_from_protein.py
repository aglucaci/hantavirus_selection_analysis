#!/usr/bin/env python3
"""
Generate a codon-aligned CDS FASTA from:
  1. an aligned protein FASTA
  2. matching raw, ungapped CDS FASTA

For each protein alignment column, amino-acid gaps become codon gaps (---), and
amino-acid residues consume the next codon from the matching CDS sequence.
"""

import argparse
import csv
import gzip
import re
import sys
from collections import Counter
from pathlib import Path


CODON_TABLE = {
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


def log(message, verbose=True):
    if verbose:
        print(message, flush=True)


def open_text(path, mode="rt"):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, mode)
    return open(path, mode, newline="" if "w" in mode else None)


def clean_id(seq_id):
    """Normalize IDs so matching records like PP_123 and PP_123|M pair up."""
    seq_id = seq_id.strip().split()[0]

    if seq_id.endswith("|M"):
        seq_id = seq_id[:-2]

    return seq_id.split("|")[0]


def clean_nt(seq):
    seq = str(seq).upper().replace("U", "T")
    return re.sub(r"[^ACGTN-]", "", seq).replace("-", "")


def clean_aa(seq):
    seq = str(seq).upper().replace(".", "-")
    return re.sub(r"[^A-ZX*\-]", "", seq)


def read_fasta(path):
    records = []
    with open_text(path, "rt") as handle:
        current_id = None
        current_desc = ""
        current_seq = []

        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    records.append(
                        {
                            "id": current_id,
                            "description": current_desc,
                            "seq": "".join(current_seq),
                        }
                    )
                current_desc = line[1:].strip()
                current_id = current_desc.split()[0]
                current_seq = []
            else:
                current_seq.append(line)

        if current_id is not None:
            records.append(
                {
                    "id": current_id,
                    "description": current_desc,
                    "seq": "".join(current_seq),
                }
            )

    return records


def write_fasta(records, path):
    with open_text(path, "wt") as handle:
        for record in records:
            description = f" {record['description']}" if record.get("description") else ""
            handle.write(f">{record['id']}{description}\n")
            seq = record["seq"]
            for i in range(0, len(seq), 70):
                handle.write(seq[i : i + 70] + "\n")


def write_qc(rows, path):
    fieldnames = [
        "id",
        "status",
        "protein_alignment_length_aa",
        "ungapped_protein_length_aa",
        "cds_length_nt",
        "available_codons",
        "used_codons",
        "unused_codons",
        "missing_codons_inserted",
        "terminal_stop_skipped",
        "padded_terminal_nt",
        "translation_matches",
        "translation_mismatches",
        "translation_identity",
        "codon_alignment_length_nt",
    ]

    with open_text(path, "wt") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def pad_to_codon_length(nt_seq):
    remainder = len(nt_seq) % 3
    if remainder == 0:
        return nt_seq, 0

    n_pad = 3 - remainder
    return nt_seq + ("N" * n_pad), n_pad


def translate_codon(codon):
    if "-" in codon:
        return "-"
    return CODON_TABLE.get(codon, "X")


def score_translation(aligned_protein, codon_alignment):
    matches = 0
    mismatches = 0

    for i, aa in enumerate(aligned_protein):
        codon = codon_alignment[i * 3 : i * 3 + 3]
        translated = translate_codon(codon).replace("*", "X")
        aa = "X" if aa == "*" else aa

        if aa == "-":
            continue
        if codon == "---":
            continue
        if aa == translated or aa == "X" or translated == "X":
            matches += 1
        else:
            mismatches += 1

    compared = matches + mismatches
    return matches, mismatches, matches / compared if compared else 0.0


def codon_align(aligned_protein, cds_seq, missing_as):
    aligned_protein = clean_aa(aligned_protein)
    cds_seq = clean_nt(cds_seq)
    cds_seq, padded_terminal_nt = pad_to_codon_length(cds_seq)
    codons = [cds_seq[i : i + 3] for i in range(0, len(cds_seq), 3)]

    codon_index = 0
    codon_alignment = []
    missing_codons_inserted = 0
    terminal_stop_skipped = 0

    for aa_index, aa in enumerate(aligned_protein):
        if aa == "-":
            codon_alignment.append("---")
            continue

        is_terminal_stop = aa == "*" and not aligned_protein[aa_index + 1 :].replace("-", "")
        if is_terminal_stop and codon_index >= len(codons):
            terminal_stop_skipped += 1
            continue

        if codon_index < len(codons):
            codon_alignment.append(codons[codon_index])
            codon_index += 1
        else:
            missing_codons_inserted += 1
            if missing_as == "NNN":
                codon_alignment.append("NNN")
            elif missing_as == "gap":
                codon_alignment.append("---")
            elif missing_as == "skip":
                continue

    codon_alignment = "".join(codon_alignment)
    matches, mismatches, identity = score_translation(aligned_protein, codon_alignment)

    return {
        "codon_alignment": codon_alignment,
        "protein_alignment_length_aa": len(aligned_protein),
        "ungapped_protein_length_aa": len(aligned_protein.replace("-", "")),
        "cds_length_nt": len(clean_nt(cds_seq)),
        "available_codons": len(codons),
        "used_codons": codon_index,
        "unused_codons": len(codons) - codon_index,
        "missing_codons_inserted": missing_codons_inserted,
        "terminal_stop_skipped": terminal_stop_skipped,
        "padded_terminal_nt": padded_terminal_nt,
        "translation_matches": matches,
        "translation_mismatches": mismatches,
        "translation_identity": identity,
        "codon_alignment_length_nt": len(codon_alignment),
    }


def blank_qc_row(seq_id, status, protein_length="", ungapped_length=""):
    return {
        "id": seq_id,
        "status": status,
        "protein_alignment_length_aa": protein_length,
        "ungapped_protein_length_aa": ungapped_length,
        "cds_length_nt": "",
        "available_codons": "",
        "used_codons": "",
        "unused_codons": "",
        "missing_codons_inserted": "",
        "terminal_stop_skipped": "",
        "padded_terminal_nt": "",
        "translation_matches": "",
        "translation_mismatches": "",
        "translation_identity": "",
        "codon_alignment_length_nt": "",
    }


def main():
    parser = argparse.ArgumentParser(
        description="Generate codon-aligned CDS FASTA from aligned protein FASTA and raw CDS FASTA."
    )
    parser.add_argument("-p", "--protein-alignment", required=True, help="Aligned protein FASTA")
    parser.add_argument("-c", "--cds-fasta", required=True, help="Raw ungapped CDS FASTA")
    parser.add_argument("-o", "--output", required=True, help="Output codon-aligned CDS FASTA")
    parser.add_argument("--qc", default="codon_alignment.qc.tsv", help="Output QC TSV")
    parser.add_argument(
        "--min-translation-identity",
        type=float,
        default=0.70,
        help="Minimum ungapped protein-vs-CDS translation identity to write. Default: 0.70",
    )
    parser.add_argument(
        "--missing-as",
        choices=["NNN", "gap", "skip"],
        default="gap",
        help="How to represent protein residues after CDS codons run out. Default: gap",
    )
    parser.add_argument("--verbose", action="store_true", help="Print progress and per-record details.")
    parser.add_argument(
        "--progress-every",
        type=int,
        default=25,
        help="With --verbose, report progress after this many protein records. Default: 25",
    )
    args = parser.parse_args()

    log("Starting codon alignment generation", args.verbose)
    log(f"Protein alignment: {args.protein_alignment}", args.verbose)
    log(f"Raw CDS FASTA:     {args.cds_fasta}", args.verbose)
    log(f"Output alignment:  {args.output}", args.verbose)
    log(f"QC report:         {args.qc}", args.verbose)
    log(f"Minimum translation identity: {args.min_translation_identity:.3f}", args.verbose)

    protein_records = read_fasta(args.protein_alignment)
    cds_records = read_fasta(args.cds_fasta)
    cds_by_id = {clean_id(record["id"]): record for record in cds_records}

    log(f"Loaded protein records: {len(protein_records)}", args.verbose)
    log(f"Loaded CDS records:     {len(cds_records)}", args.verbose)
    log(f"Unique normalized CDS IDs: {len(cds_by_id)}", args.verbose)

    output_records = []
    qc_rows = []
    status_counts = Counter()

    for record_index, protein_record in enumerate(protein_records, start=1):
        seq_id = clean_id(protein_record["id"])
        aligned_protein = clean_aa(protein_record["seq"])
        ungapped_length = len(aligned_protein.replace("-", ""))

        if seq_id not in cds_by_id:
            status = "missing_cds_record"
            status_counts[status] += 1
            qc_rows.append(blank_qc_row(seq_id, status, len(aligned_protein), ungapped_length))
            log(f"[{record_index}/{len(protein_records)}] {seq_id}: {status}", args.verbose)
            continue

        result = codon_align(aligned_protein, cds_by_id[seq_id]["seq"], args.missing_as)
        status = (
            "written"
            if result["translation_identity"] >= args.min_translation_identity
            else "low_translation_identity_skipped"
        )
        status_counts[status] += 1

        if status == "written":
            output_records.append(
                {
                    "id": seq_id,
                    "seq": result["codon_alignment"],
                    "description": (
                        "codon_aligned_CDS "
                        f"translation_identity={result['translation_identity']:.4f}"
                    ),
                }
            )

        qc_row = {
            "id": seq_id,
            "status": status,
            **{key: value for key, value in result.items() if key != "codon_alignment"},
        }
        qc_row["translation_identity"] = round(qc_row["translation_identity"], 6)
        qc_rows.append(qc_row)

        log(
            f"[{record_index}/{len(protein_records)}] {seq_id}: {status} "
            f"identity={result['translation_identity']:.4f} "
            f"used_codons={result['used_codons']} unused_codons={result['unused_codons']} "
            f"missing_codons={result['missing_codons_inserted']} "
            f"aligned_nt={result['codon_alignment_length_nt']}",
            args.verbose,
        )

        if (
            args.verbose
            and args.progress_every > 0
            and record_index % args.progress_every == 0
        ):
            log(
                "Progress summary: "
                + ", ".join(f"{status}={count}" for status, count in sorted(status_counts.items())),
                True,
            )

    write_fasta(output_records, args.output)
    write_qc(qc_rows, args.qc)

    print(f"Saved codon alignment FASTA: {args.output}")
    print(f"Saved QC report: {args.qc}")
    print(f"Sequences written: {len(output_records)}")
    print("QC status counts:")
    for status, count in sorted(status_counts.items()):
        print(f"  {status}: {count}")

    if not output_records:
        print("ERROR: no sequences written.", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
