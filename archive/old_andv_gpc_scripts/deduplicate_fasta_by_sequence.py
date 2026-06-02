#!/usr/bin/env python3
"""
Deduplicate a FASTA by exact sequence identity.

Writes:
  - representative FASTA
  - duplicate map TSV

The first occurrence of each unique sequence is kept as the representative.
"""

import argparse
import csv
from collections import OrderedDict


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
                            "seq": "".join(current_seq).upper(),
                        }
                    )
                current_desc = line[1:].strip()
                current_id = current_desc.split()[0]
                current_seq = []
            else:
                current_seq.append(line)

    if current_id is not None:
        records.append({"id": current_id, "description": current_desc, "seq": "".join(current_seq).upper()})

    return records


def write_fasta(records, path):
    with open(path, "w") as handle:
        for record in records:
            handle.write(f">{record['representative_id']}\n")
            seq = record["seq"]
            for i in range(0, len(seq), 70):
                handle.write(seq[i : i + 70] + "\n")


def main():
    parser = argparse.ArgumentParser(description="Deduplicate FASTA records by exact sequence.")
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("--map", required=True, help="Output duplicate map TSV")
    parser.add_argument(
        "--keep-id",
        action="append",
        default=[],
        help="Sequence ID to force as representative if it is present in a duplicate group. Can be repeated.",
    )
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    records = read_fasta(args.input)
    keep_ids = set(args.keep_id)
    by_sequence = OrderedDict()

    for record in records:
        seq = record["seq"]
        if seq not in by_sequence:
            by_sequence[seq] = {
                "representative_id": record["id"],
                "seq": seq,
                "duplicates": [],
            }
        by_sequence[seq]["duplicates"].append(record["id"])

    representatives = list(by_sequence.values())

    for group in representatives:
        preferred = [seq_id for seq_id in group["duplicates"] if seq_id in keep_ids]
        if preferred:
            group["representative_id"] = preferred[0]

    write_fasta(representatives, args.output)

    with open(args.map, "w", newline="") as handle:
        fieldnames = [
            "representative_id",
            "member_id",
            "member_role",
            "duplicate_count",
            "sequence_length_nt",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()

        for group in representatives:
            duplicate_count = len(group["duplicates"])
            for member_id in group["duplicates"]:
                writer.writerow(
                    {
                        "representative_id": group["representative_id"],
                        "member_id": member_id,
                        "member_role": "representative" if member_id == group["representative_id"] else "duplicate",
                        "duplicate_count": duplicate_count,
                        "sequence_length_nt": len(group["seq"]),
                    }
                )

    print(f"Saved deduplicated FASTA: {args.output}")
    print(f"Saved duplicate map: {args.map}")
    print(f"Input records: {len(records)}")
    print(f"Unique sequences: {len(representatives)}")
    print(f"Duplicate records removed: {len(records) - len(representatives)}")

    if args.verbose:
        groups_with_duplicates = [group for group in representatives if len(group["duplicates"]) > 1]
        print(f"Duplicate groups: {len(groups_with_duplicates)}")
        for group in groups_with_duplicates[:20]:
            print(f"{group['representative_id']}: {','.join(group['duplicates'])}")


if __name__ == "__main__":
    main()
