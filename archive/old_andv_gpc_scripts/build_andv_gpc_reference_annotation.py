#!/usr/bin/env python3
"""
Build the first-pass ANDV GPC reference annotation files.

Reference:
  NC_003467.2 / CHI-7913 Andes virus segment M

This script uses the local Pathoplexus RefSeq row when available
PP_006W0E7.1 == NC_003467.2.M, plus the recovered GPC CDS, to write:
  - data/reference/NC_003467.2_M.fasta
  - data/reference/NC_003467.2_GPC_CDS.no_stop.fasta
  - data/reference/NC_003467.2_GPC_protein.no_stop.fasta
  - data/reference/ANDV_GPC_reference_annotation.tsv
"""

import argparse
import csv
import re
import sys
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


def clean_id(seq_id):
    seq_id = seq_id.strip().split()[0]
    if seq_id.endswith("|M"):
        seq_id = seq_id[:-2]
    return seq_id.split("|")[0]


def clean_nt(seq):
    return re.sub(r"[^ACGTN]", "", seq.upper().replace("U", "T"))


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
        records.append({"id": current_id, "description": current_desc, "seq": "".join(current_seq)})

    return records


def write_fasta(path, seq_id, description, seq):
    with open(path, "w") as handle:
        handle.write(f">{seq_id} {description}\n")
        for i in range(0, len(seq), 70):
            handle.write(seq[i : i + 70] + "\n")


def translate(seq):
    seq = clean_nt(seq)
    seq = seq[: len(seq) - (len(seq) % 3)]
    return "".join(CODON_TABLE.get(seq[i : i + 3], "X") for i in range(0, len(seq), 3))


def nt_start_for_aa(cds_start_nt, aa_pos):
    return cds_start_nt + ((aa_pos - 1) * 3)


def nt_end_for_aa(cds_start_nt, aa_pos):
    return nt_start_for_aa(cds_start_nt, aa_pos) + 2


def annotation_row(
    feature,
    region_type,
    start_nt,
    end_nt,
    start_codon,
    end_codon,
    start_aa,
    end_aa,
    evidence,
    notes,
):
    return {
        "feature": feature,
        "region_type": region_type,
        "ref_accession": "NC_003467.2",
        "segment": "M",
        "protein": "GPC",
        "start_nt": start_nt,
        "end_nt": end_nt,
        "start_codon": start_codon,
        "end_codon": end_codon,
        "start_aa": start_aa,
        "end_aa": end_aa,
        "evidence": evidence,
        "notes": notes,
    }


def main():
    parser = argparse.ArgumentParser(description="Build ANDV GPC reference annotation TSV and reference FASTA files.")
    parser.add_argument(
        "--m-fasta",
        default="data/andv_nuc-M_2026-05-11T1953.fasta/andv_nuc-M_2026-05-11T1953.fasta",
        help="Raw M-segment nucleotide FASTA containing PP_006W0E7.1 / NC_003467.2",
    )
    parser.add_argument(
        "--cds-fasta",
        default="results/GPC_CDS.recovered.raw.fasta",
        help="Recovered raw CDS FASTA containing PP_006W0E7.1",
    )
    parser.add_argument("--ref-id", default="PP_006W0E7.1", help="Local ID corresponding to NC_003467.2.M")
    parser.add_argument("--outdir", default="data/reference")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    m_records = {clean_id(record["id"]): record for record in read_fasta(args.m_fasta)}
    cds_records = {clean_id(record["id"]): record for record in read_fasta(args.cds_fasta)}
    ref_id = clean_id(args.ref_id)

    if ref_id not in m_records:
        raise SystemExit(f"ERROR: reference M segment record not found in {args.m_fasta}: {ref_id}")
    if ref_id not in cds_records:
        raise SystemExit(f"ERROR: reference recovered CDS record not found in {args.cds_fasta}: {ref_id}")

    m_seq = clean_nt(m_records[ref_id]["seq"])
    cds_no_stop = clean_nt(cds_records[ref_id]["seq"])
    cds_start_index = m_seq.find(cds_no_stop)

    if cds_start_index < 0:
        raise SystemExit("ERROR: recovered CDS was not found as a contiguous substring of the reference M segment.")

    cds_start_nt = cds_start_index + 1
    cds_end_nt_no_stop = cds_start_nt + len(cds_no_stop) - 1
    stop_start_nt = cds_end_nt_no_stop + 1
    stop_codon = m_seq[stop_start_nt - 1 : stop_start_nt + 2]
    cds_with_stop_end_nt = cds_end_nt_no_stop + 3 if stop_codon in {"TAA", "TAG", "TGA"} else cds_end_nt_no_stop
    protein = translate(cds_no_stop)

    waasa_start_aa = protein.find("WAASA") + 1
    if waasa_start_aa <= 0:
        raise SystemExit("ERROR: WAASA motif not found in reference GPC protein.")
    waasa_end_aa = waasa_start_aa + 4
    gpc_end_aa = len(protein)

    rows = [
        annotation_row(
            "GPC_CDS",
            "CDS",
            cds_start_nt,
            cds_with_stop_end_nt,
            1,
            gpc_end_aa + (1 if cds_with_stop_end_nt > cds_end_nt_no_stop else 0),
            1,
            gpc_end_aa,
            "NC_003467.2 local Pathoplexus RefSeq sequence; recovered CDS substring; terminal stop codon checked",
            f"coding region without stop is nt {cds_start_nt}-{cds_end_nt_no_stop}; stop_codon={stop_codon}",
        ),
        annotation_row(
            "GPC",
            "protein_region",
            cds_start_nt,
            cds_end_nt_no_stop,
            1,
            gpc_end_aa,
            1,
            gpc_end_aa,
            "NC_003467.2 GPC CDS translated from recovered no-stop CDS",
            "Full glycoprotein precursor, no terminal stop included in protein coordinates",
        ),
        annotation_row(
            "Gn",
            "protein_region",
            nt_start_for_aa(cds_start_nt, 1),
            nt_end_for_aa(cds_start_nt, waasa_end_aa),
            1,
            waasa_end_aa,
            1,
            waasa_end_aa,
            "WAASA cleavage motif split",
            "First-pass Gn assignment includes WAASA motif through its final alanine",
        ),
        annotation_row(
            "Gc",
            "protein_region",
            nt_start_for_aa(cds_start_nt, waasa_end_aa + 1),
            nt_end_for_aa(cds_start_nt, gpc_end_aa),
            waasa_end_aa + 1,
            gpc_end_aa,
            waasa_end_aa + 1,
            gpc_end_aa,
            "WAASA cleavage motif split",
            "First-pass Gc assignment begins after WAASA motif",
        ),
        annotation_row(
            "WAASA_cleavage_motif",
            "motif",
            nt_start_for_aa(cds_start_nt, waasa_start_aa),
            nt_end_for_aa(cds_start_nt, waasa_end_aa),
            waasa_start_aa,
            waasa_end_aa,
            waasa_start_aa,
            waasa_end_aa,
            "Motif search on NC_003467.2 GPC protein",
            "Conserved hantavirus GPC cleavage motif WAASA",
        ),
    ]

    for feature in [
        "Gn_ectodomain",
        "Gn_TM",
        "Gn_cytoplasmic_tail",
        "Gc_ectodomain",
        "Gc_TM",
        "Gc_cytoplasmic_tail",
    ]:
        rows.append(
            annotation_row(
                feature,
                "pending_region",
                "",
                "",
                "",
                "",
                "",
                "",
                "pending DeepTMHMM/TMHMM/Phobius prediction on NC_003467.2 GPC protein",
                "Fill these coordinates after transmembrane/topology prediction",
            )
        )

    annotation_path = outdir / "ANDV_GPC_reference_annotation.tsv"
    with open(annotation_path, "w", newline="") as handle:
        fieldnames = [
            "feature",
            "region_type",
            "ref_accession",
            "segment",
            "protein",
            "start_nt",
            "end_nt",
            "start_codon",
            "end_codon",
            "start_aa",
            "end_aa",
            "evidence",
            "notes",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    write_fasta(outdir / "NC_003467.2_M.fasta", "NC_003467.2", "Andes virus segment M complete genome from local PP_006W0E7.1", m_seq)
    write_fasta(outdir / "NC_003467.2_GPC_CDS.no_stop.fasta", "NC_003467.2_GPC_CDS_no_stop", f"nt {cds_start_nt}-{cds_end_nt_no_stop}", cds_no_stop)
    write_fasta(outdir / "NC_003467.2_GPC_protein.no_stop.fasta", "NC_003467.2_GPC_no_stop", f"aa 1-{gpc_end_aa}", protein)

    print(f"Saved annotation TSV: {annotation_path}")
    print(f"Saved reference FASTA files in: {outdir}")
    print(f"NC_003467.2 M length: {len(m_seq)} nt")
    print(f"GPC CDS no-stop: nt {cds_start_nt}-{cds_end_nt_no_stop} ({len(cds_no_stop)} nt)")
    print(f"Stop codon: nt {stop_start_nt}-{stop_start_nt + 2} {stop_codon}")
    print(f"GPC protein length: {gpc_end_aa} aa")
    print(f"WAASA motif: aa {waasa_start_aa}-{waasa_end_aa}")


if __name__ == "__main__":
    main()
