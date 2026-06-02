#!/usr/bin/env python3
"""
Annotate HyPhy selected-site tables with NC_003467.2 GPC reference coordinates.

Input should be a TSV/CSV with a site column, for example:
  hyphy_site,test,p_value,q_value,omega,evidence_type

Output adds:
  reference_codon, reference_nt_start, reference_nt_end, aa_ref, region
"""

import argparse
import csv
import re


def sniff_delimiter(path):
    if path.lower().endswith(".csv"):
        return ","
    return "\t"


def read_fasta_seq(path):
    seq = []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            seq.append(line)
    return "".join(seq)


def read_annotation(path):
    rows = []
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            try:
                row["_start_aa"] = int(row["start_aa"])
                row["_end_aa"] = int(row["end_aa"])
            except (TypeError, ValueError):
                continue
            rows.append(row)
    return rows


def find_cds_start(annotation_rows):
    for row in annotation_rows:
        if row["feature"] == "GPC":
            return int(row["start_nt"])
    for row in annotation_rows:
        if row["feature"] == "GPC_CDS":
            return int(row["start_nt"])
    raise ValueError("Could not find GPC or GPC_CDS start_nt in annotation table.")


def region_for_site(site, annotation_rows):
    hits = []
    priority = {"motif": 0, "protein_region": 1, "CDS": 2}

    for row in annotation_rows:
        if row["_start_aa"] <= site <= row["_end_aa"]:
            hits.append(row)

    if not hits:
        return ""

    hits.sort(key=lambda row: (priority.get(row["region_type"], 9), row["_end_aa"] - row["_start_aa"]))
    return ";".join(row["feature"] for row in hits if row["feature"] != "GPC_CDS")


def parse_site(value):
    match = re.search(r"\d+", str(value))
    if not match:
        raise ValueError(f"Could not parse site value: {value}")
    return int(match.group(0))


def main():
    parser = argparse.ArgumentParser(description="Annotate HyPhy selected sites with ANDV GPC reference coordinates.")
    parser.add_argument("-i", "--input", required=True, help="Input selected-site TSV/CSV")
    parser.add_argument("-o", "--output", required=True, help="Output annotated TSV")
    parser.add_argument(
        "--annotation",
        default="data/reference/ANDV_GPC_reference_annotation.tsv",
        help="Reference annotation TSV",
    )
    parser.add_argument(
        "--reference-protein",
        default="data/reference/NC_003467.2_GPC_protein.no_stop.fasta",
        help="Reference GPC protein FASTA",
    )
    parser.add_argument("--site-column", default="hyphy_site")
    args = parser.parse_args()

    annotation_rows = read_annotation(args.annotation)
    cds_start_nt = find_cds_start(annotation_rows)
    protein = read_fasta_seq(args.reference_protein)
    delimiter = sniff_delimiter(args.input)

    with open(args.input, newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if args.site_column not in reader.fieldnames:
            raise SystemExit(f"ERROR: site column not found: {args.site_column}")
        rows = list(reader)

    extra_fields = [
        "reference_codon",
        "reference_nt_start",
        "reference_nt_end",
        "aa_ref",
        "region",
    ]
    fieldnames = list(rows[0].keys()) + [field for field in extra_fields if field not in rows[0]]

    with open(args.output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()

        for row in rows:
            site = parse_site(row[args.site_column])
            nt_start = cds_start_nt + ((site - 1) * 3)
            nt_end = nt_start + 2
            row["reference_codon"] = site
            row["reference_nt_start"] = nt_start
            row["reference_nt_end"] = nt_end
            row["aa_ref"] = protein[site - 1] if 1 <= site <= len(protein) else ""
            row["region"] = region_for_site(site, annotation_rows)
            writer.writerow(row)

    print(f"Saved annotated selected-site table: {args.output}")


if __name__ == "__main__":
    main()
