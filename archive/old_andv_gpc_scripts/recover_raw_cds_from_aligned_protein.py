#!/usr/bin/env python3
"""
Recover raw, ungapped CDS FASTA records from:
  1. an aligned protein FASTA
  2. matching raw nucleotide FASTA records, such as full M segment sequences

For each protein record, this script removes alignment gaps, finds the best
matching region in the raw nucleotide record across all six reading frames,
and writes the corresponding coding nucleotide sequence.
"""

import argparse
import csv
import gzip
import re
import sys
from collections import Counter
from pathlib import Path


GENETIC_CODE_TABLE = 1
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
RC_TABLE = str.maketrans("ACGTN", "TGCAN")


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
    return re.sub(r"[^ACGTN]", "", seq)


def clean_protein(seq):
    seq = str(seq).upper()
    seq = seq.replace("-", "").replace(".", "").replace("*", "")
    return re.sub(r"[^A-ZX]", "", seq)


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
    return len(records)


def write_qc(rows, path):
    fieldnames = [
        "id",
        "status",
        "protein_length_aa",
        "nt_length",
        "strand",
        "frame_0based",
        "score",
        "identity",
        "matches",
        "mismatches",
        "codons_used",
        "missing_codons",
        "cds_length_nt",
    ]

    with open_text(path, "wt") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def codonize(nt_seq, frame):
    usable = nt_seq[frame:]
    usable = usable[: len(usable) - (len(usable) % 3)]
    return [usable[i : i + 3] for i in range(0, len(usable), 3)]


def translate_codons(codons):
    amino_acids = []
    for codon in codons:
        amino_acids.append(CODON_TABLE.get(codon, "X").replace("*", "X"))
    return "".join(amino_acids)


def reverse_complement(nt_seq):
    return nt_seq.translate(RC_TABLE)[::-1]


def six_frame_translations(nt_seq):
    nt_seq = clean_nt(nt_seq)
    rc_seq = reverse_complement(nt_seq)

    for strand, seq in [("+", nt_seq), ("-", rc_seq)]:
        for frame in [0, 1, 2]:
            codons = codonize(seq, frame)
            yield {
                "strand": strand,
                "frame": frame,
                "codons": codons,
                "translation": translate_codons(codons),
            }


def score_window(protein, translation_window):
    matches = 0
    mismatches = 0
    for p_aa, t_aa in zip(protein, translation_window):
        if p_aa == t_aa or p_aa == "X" or t_aa == "X":
            matches += 1
        else:
            mismatches += 1

    compared = matches + mismatches
    return matches, mismatches, matches / compared if compared else 0.0


def protein_seed_positions(protein, seed_len):
    seeds = []
    seen = set()
    max_positions = 80

    for i in range(0, max(0, len(protein) - seed_len + 1), seed_len):
        seed = protein[i : i + seed_len]
        if "X" in seed or seed in seen:
            continue
        seen.add(seed)
        seeds.append((seed, i))
        if len(seeds) >= max_positions:
            break

    return seeds


def best_contiguous_cds(protein, frame_info, seed_len=8):
    """Find the best same-length translated window in one frame."""
    translation = frame_info["translation"]
    protein_len = len(protein)
    if protein_len == 0 or len(translation) < protein_len:
        return None

    candidate_starts = set()
    for seed, protein_pos in protein_seed_positions(protein, seed_len):
        search_from = 0
        while True:
            hit = translation.find(seed, search_from)
            if hit == -1:
                break
            start = hit - protein_pos
            if 0 <= start <= len(translation) - protein_len:
                candidate_starts.add(start)
            search_from = hit + 1

    if not candidate_starts:
        candidate_starts = range(0, len(translation) - protein_len + 1)

    best = None
    for start in candidate_starts:
        window = translation[start : start + protein_len]
        matches, mismatches, identity = score_window(protein, window)
        cds = "".join(frame_info["codons"][start : start + protein_len])
        result = {
            **frame_info,
            "score": matches * 3 - mismatches * 2,
            "identity": identity,
            "matches": matches,
            "mismatches": mismatches,
            "codons_used": protein_len,
            "missing_codons": 0,
            "cds_length_nt": len(cds),
            "cds": cds,
        }
        rank = (result["identity"], result["score"], result["cds_length_nt"])
        if best is None or rank > best["_rank"]:
            best = result
            best["_rank"] = rank

    return best


def score_alignment(protein, translation):
    alignment = smith_waterman(protein, translation)
    if alignment is None:
        return None

    return {
        "protein_aln": alignment["protein_aln"],
        "translation_aln": alignment["translation_aln"],
        "score": alignment["score"],
    }


def smith_waterman(protein, translation):
    """Small standalone local aligner with linear gap scoring."""
    if not protein or not translation:
        return None

    match_score = 3
    mismatch_score = -2
    gap_score = -4

    n_cols = len(translation) + 1
    previous = [0] * n_cols
    trace = [bytearray(n_cols) for _ in range(len(protein) + 1)]
    best_score = 0
    best_i = 0
    best_j = 0

    for i, p_aa in enumerate(protein, start=1):
        current = [0] * n_cols
        row_trace = trace[i]

        for j, t_aa in enumerate(translation, start=1):
            diagonal_score = previous[j - 1] + (
                match_score if p_aa == t_aa or p_aa == "X" or t_aa == "X" else mismatch_score
            )
            up_score = previous[j] + gap_score
            left_score = current[j - 1] + gap_score
            score = max(0, diagonal_score, up_score, left_score)
            current[j] = score

            if score == 0:
                row_trace[j] = 0
            elif score == diagonal_score:
                row_trace[j] = 1
            elif score == up_score:
                row_trace[j] = 2
            else:
                row_trace[j] = 3

            if score > best_score:
                best_score = score
                best_i = i
                best_j = j

        previous = current

    if best_score <= 0:
        return None

    protein_aln = []
    translation_aln = []
    i = best_i
    j = best_j

    while i > 0 and j > 0:
        direction = trace[i][j]
        if direction == 0:
            break
        if direction == 1:
            protein_aln.append(protein[i - 1])
            translation_aln.append(translation[j - 1])
            i -= 1
            j -= 1
        elif direction == 2:
            protein_aln.append(protein[i - 1])
            translation_aln.append("-")
            i -= 1
        elif direction == 3:
            protein_aln.append("-")
            translation_aln.append(translation[j - 1])
            j -= 1

    return {
        "protein_aln": "".join(reversed(protein_aln)),
        "translation_aln": "".join(reversed(translation_aln)),
        "score": best_score,
    }


def reconstruct_cds(protein_aln, translation_aln, codons, missing_as):
    cds_codons = []
    translated_index = 0
    matches = 0
    mismatches = 0
    codons_used = 0
    missing_codons = 0

    for p_aa, t_aa in zip(protein_aln, translation_aln):
        codon = None
        if t_aa != "-":
            codon = codons[translated_index]
            translated_index += 1

        if p_aa == "-" and t_aa != "-":
            continue

        if p_aa != "-" and t_aa == "-":
            missing_codons += 1
            if missing_as == "NNN":
                cds_codons.append("NNN")
            elif missing_as == "gap":
                cds_codons.append("---")
            elif missing_as == "skip":
                pass
            continue

        if p_aa != "-" and t_aa != "-":
            cds_codons.append(codon)
            codons_used += 1
            if p_aa == t_aa or p_aa == "X" or t_aa == "X":
                matches += 1
            else:
                mismatches += 1

    compared = matches + mismatches
    cds = "".join(cds_codons)
    return {
        "cds": cds,
        "matches": matches,
        "mismatches": mismatches,
        "identity": matches / compared if compared else 0.0,
        "codons_used": codons_used,
        "missing_codons": missing_codons,
        "cds_length_nt": len(cds),
    }


def recover_best_cds(protein_seq, nt_seq, missing_as, local_fallback=False):
    best = None

    for frame_info in six_frame_translations(nt_seq):
        contiguous = best_contiguous_cds(protein_seq, frame_info)
        if contiguous is not None:
            rank = (contiguous["identity"], contiguous["score"], contiguous["cds_length_nt"])
            if best is None or rank > best["_rank"]:
                best = contiguous
                best["_rank"] = rank

    if best is not None or not local_fallback:
        return best

    for frame_info in six_frame_translations(nt_seq):
        alignment = score_alignment(protein_seq, frame_info["translation"])
        if alignment is None:
            continue

        reconstructed = reconstruct_cds(
            alignment["protein_aln"],
            alignment["translation_aln"],
            frame_info["codons"],
            missing_as,
        )

        result = {**frame_info, **alignment, **reconstructed}
        rank = (result["identity"], result["score"], result["cds_length_nt"])

        if best is None or rank > best["_rank"]:
            best = result
            best["_rank"] = rank

    return best


def blank_qc_row(seq_id, status, protein_length="", nt_length=""):
    return {
        "id": seq_id,
        "status": status,
        "protein_length_aa": protein_length,
        "nt_length": nt_length,
        "strand": "",
        "frame_0based": "",
        "score": "",
        "identity": "",
        "matches": "",
        "mismatches": "",
        "codons_used": "",
        "missing_codons": "",
        "cds_length_nt": "",
    }


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Recover raw ungapped CDS FASTA from aligned protein FASTA and "
            "matching raw nucleotide FASTA."
        )
    )
    parser.add_argument("-p", "--protein-alignment", required=True)
    parser.add_argument(
        "-n",
        "--nucleotide-fasta",
        required=True,
        help="Raw nucleotide FASTA. A trailing segment suffix such as '|M' is ignored for ID matching.",
    )
    parser.add_argument("-o", "--output-cds", required=True)
    parser.add_argument("--qc", default="recovered_raw_cds.qc.tsv")
    parser.add_argument(
        "--min-identity",
        type=float,
        default=0.70,
        help="Minimum protein/translation identity required to write a CDS. Default: 0.70",
    )
    parser.add_argument(
        "--missing-as",
        choices=["NNN", "gap", "skip"],
        default="NNN",
        help="How to represent protein residues without a nucleotide codon. Default: NNN",
    )
    parser.add_argument(
        "--local-fallback",
        action="store_true",
        help="Use slower Smith-Waterman fallback if no contiguous translated window is found.",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print progress and per-record recovery details.",
    )
    parser.add_argument(
        "--progress-every",
        type=int,
        default=25,
        help="With --verbose, report progress after this many protein records. Default: 25",
    )
    args = parser.parse_args()

    log("Starting raw CDS recovery", args.verbose)
    log(f"Protein alignment: {args.protein_alignment}", args.verbose)
    log(f"Nucleotide FASTA:   {args.nucleotide_fasta}", args.verbose)
    log(f"Output CDS FASTA:   {args.output_cds}", args.verbose)
    log(f"QC report:          {args.qc}", args.verbose)
    log(f"Minimum identity:   {args.min_identity:.3f}", args.verbose)
    log(f"Local fallback:     {args.local_fallback}", args.verbose)

    protein_records = read_fasta(args.protein_alignment)
    nucleotide_records = read_fasta(args.nucleotide_fasta)
    nucleotide_by_id = {clean_id(rec["id"]): rec for rec in nucleotide_records}

    log(f"Loaded protein records:    {len(protein_records)}", args.verbose)
    log(f"Loaded nucleotide records: {len(nucleotide_records)}", args.verbose)
    log(f"Unique normalized nucleotide IDs: {len(nucleotide_by_id)}", args.verbose)

    output_records = []
    qc_rows = []
    status_counts = Counter()

    for record_index, protein_record in enumerate(protein_records, start=1):
        seq_id = clean_id(protein_record["id"])
        protein_seq = clean_protein(protein_record["seq"])

        if seq_id not in nucleotide_by_id:
            status = "missing_nucleotide_record"
            status_counts[status] += 1
            qc_rows.append(blank_qc_row(seq_id, status, len(protein_seq)))
            log(f"[{record_index}/{len(protein_records)}] {seq_id}: {status}", args.verbose)
            continue

        nt_seq = clean_nt(nucleotide_by_id[seq_id]["seq"])
        if not protein_seq or not nt_seq:
            status = "empty_sequence"
            status_counts[status] += 1
            qc_rows.append(blank_qc_row(seq_id, status, len(protein_seq), len(nt_seq)))
            log(
                f"[{record_index}/{len(protein_records)}] {seq_id}: {status} "
                f"protein_aa={len(protein_seq)} nt={len(nt_seq)}",
                args.verbose,
            )
            continue

        best = recover_best_cds(protein_seq, nt_seq, args.missing_as, args.local_fallback)
        if best is None:
            status = "no_alignment"
            status_counts[status] += 1
            qc_rows.append(blank_qc_row(seq_id, status, len(protein_seq), len(nt_seq)))
            log(
                f"[{record_index}/{len(protein_records)}] {seq_id}: {status} "
                f"protein_aa={len(protein_seq)} nt={len(nt_seq)}",
                args.verbose,
            )
            continue

        status = "written" if best["identity"] >= args.min_identity else "low_identity_skipped"
        status_counts[status] += 1

        if status == "written":
            output_records.append(
                {
                    "id": seq_id,
                    "seq": best["cds"],
                    "description": (
                        f"raw_CDS strand={best['strand']} frame={best['frame']} "
                        f"identity={best['identity']:.4f}"
                    ),
                }
            )

        qc_rows.append(
            {
                "id": seq_id,
                "status": status,
                "protein_length_aa": len(protein_seq),
                "nt_length": len(nt_seq),
                "strand": best["strand"],
                "frame_0based": best["frame"],
                "score": round(best["score"], 3),
                "identity": round(best["identity"], 6),
                "matches": best["matches"],
                "mismatches": best["mismatches"],
                "codons_used": best["codons_used"],
                "missing_codons": best["missing_codons"],
                "cds_length_nt": best["cds_length_nt"],
            }
        )

        log(
            f"[{record_index}/{len(protein_records)}] {seq_id}: {status} "
            f"identity={best['identity']:.4f} strand={best['strand']} "
            f"frame={best['frame']} cds_nt={best['cds_length_nt']} "
            f"matches={best['matches']} mismatches={best['mismatches']}",
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

    write_fasta(output_records, args.output_cds)
    write_qc(qc_rows, args.qc)

    print(f"Saved raw CDS FASTA: {args.output_cds}")
    print(f"Saved QC report: {args.qc}")
    print(f"Sequences written: {len(output_records)}")
    print("QC status counts:")
    for status, count in sorted(status_counts.items()):
        print(f"  {status}: {count}")

    if not output_records:
        print("ERROR: no sequences written. Try lowering --min-identity.", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
