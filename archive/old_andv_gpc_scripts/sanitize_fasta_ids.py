#!/usr/bin/env python3
"""Write a FASTA with HyPhy-safe record IDs and a mapping table."""

import argparse
import re


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


def write_fasta(records, path):
    with open(path, "w") as handle:
        for seq_id, seq in records:
            handle.write(f">{seq_id}\n")
            seq = re.sub(r"\s+", "", seq)
            for i in range(0, len(seq), 70):
                handle.write(seq[i : i + 70] + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("--map", required=True, help="Output TSV mapping original IDs to safe IDs")
    args = parser.parse_args()

    records = read_fasta(args.input)
    safe_seen = {}
    safe_records = []

    with open(args.map, "w") as id_map:
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
            safe_records.append((safe_id, seq))

    write_fasta(safe_records, args.output)
    print(f"Saved sanitized FASTA: {args.output}")
    print(f"Saved ID map: {args.map}")
    print(f"Sequences written: {len(safe_records)}")


if __name__ == "__main__":
    main()
