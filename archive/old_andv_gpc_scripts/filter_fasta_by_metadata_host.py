#!/usr/bin/env python3
"""Filter FASTA records to sequences with known hostNameScientific metadata."""

import argparse
import csv
from pathlib import Path


UNKNOWN_VALUES = {"", "na", "nan", "none", "null", "unknown", "unknown_host", "not collected", "not provided"}


def clean_id(seq_id):
    seq_id = (seq_id or "").strip().split()[0]
    if seq_id.endswith("|M"):
        seq_id = seq_id[:-2]
    return seq_id.split("|")[0]


def safe_to_metadata_id(seq_id):
    parts = seq_id.rsplit("_", 1)
    if len(parts) == 2 and parts[1].isdigit():
        return f"{parts[0]}.{parts[1]}"
    return seq_id


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
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as handle:
        for seq_id, desc, seq in records:
            handle.write(f">{desc or seq_id}\n")
            for i in range(0, len(seq), 70):
                handle.write(seq[i : i + 70] + "\n")


def read_metadata(path, id_column, host_column):
    rows_by_id = {}
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if id_column not in reader.fieldnames:
            raise SystemExit(f"ERROR: metadata ID column not found: {id_column}")
        if host_column not in reader.fieldnames:
            raise SystemExit(f"ERROR: metadata host column not found: {host_column}")
        for row in reader:
            metadata_id = clean_id(row.get(id_column, ""))
            if metadata_id:
                rows_by_id[metadata_id] = row
    return rows_by_id


def known_host(value):
    host = (value or "").strip()
    return host and host.lower() not in UNKNOWN_VALUES


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--qc", required=True)
    parser.add_argument("--id-column", default="accessionVersion")
    parser.add_argument("--host-column", default="hostNameScientific")
    parser.add_argument(
        "--keep-id",
        action="append",
        default=[],
        help="Sequence ID to retain even if host metadata is unknown. Can be repeated.",
    )
    args = parser.parse_args()

    records = read_fasta(args.input)
    metadata_by_id = read_metadata(args.metadata, args.id_column, args.host_column)
    keep_ids = {clean_id(seq_id) for seq_id in args.keep_id}
    output_records = []
    qc_rows = []

    for seq_id, desc, seq in records:
        normalized_id = clean_id(seq_id)
        candidates = [normalized_id, safe_to_metadata_id(normalized_id)]
        metadata_row = None
        for candidate in candidates:
            metadata_row = metadata_by_id.get(clean_id(candidate))
            if metadata_row:
                break

        host = metadata_row.get(args.host_column, "").strip() if metadata_row else ""
        is_known = known_host(host)
        kept_as_reference = normalized_id in keep_ids
        status = "kept_known_host" if is_known else "dropped_unknown_host"
        if kept_as_reference and not is_known:
            status = "kept_reference_unknown_host"

        if is_known or kept_as_reference:
            output_records.append((seq_id, desc, seq))

        qc_rows.append(
            {
                "sequence_id": seq_id,
                "metadata_id": metadata_row.get(args.id_column, "") if metadata_row else "",
                "hostNameScientific": host,
                "status": status,
            }
        )

    write_fasta(output_records, args.output)
    Path(args.qc).parent.mkdir(parents=True, exist_ok=True)
    with open(args.qc, "w", newline="") as handle:
        fieldnames = ["sequence_id", "metadata_id", "hostNameScientific", "status"]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(qc_rows)

    counts = {}
    for row in qc_rows:
        counts[row["status"]] = counts.get(row["status"], 0) + 1
    print(f"Input records: {len(records)}")
    print(f"Output records: {len(output_records)}")
    for status, count in sorted(counts.items()):
        print(f"{status}: {count}")


if __name__ == "__main__":
    main()
