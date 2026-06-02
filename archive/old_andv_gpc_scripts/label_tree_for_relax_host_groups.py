#!/usr/bin/env python3
"""
Label HyPhy-ready tree terminal branches for RELAX host-group tests.

Default comparison:
  Test      = human-derived tips, hostNameScientific == Homo sapiens
  Reference = reservoir-derived non-human tips with known hostNameScientific

Unknown/missing hosts are left unlabeled.
"""

import argparse
import csv
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

from label_tree_mrca_by_metadata import (  # noqa: E402
    NewickParser,
    clean_id,
    count_descendant_leaves,
    hyphy_safe_id,
    iter_nodes,
    leaf_names,
    read_newick,
    safe_label,
    write_newick,
)


UNKNOWN_VALUES = {"", "na", "nan", "none", "null", "unknown", "not collected", "not provided"}


def normalize_host(value):
    value = (value or "").strip()
    if value.lower() in UNKNOWN_VALUES:
        return ""
    return value


def read_host_metadata(path, id_column, host_column):
    rows_by_id = {}
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if id_column not in reader.fieldnames:
            raise SystemExit(f"ERROR: metadata ID column not found: {id_column}")
        if host_column not in reader.fieldnames:
            raise SystemExit(f"ERROR: metadata host column not found: {host_column}")
        for row in reader:
            raw_id = row.get(id_column, "").strip()
            if not raw_id:
                continue
            rows_by_id[clean_id(raw_id)] = row
    return rows_by_id


def safe_to_metadata_id(seq_id):
    parts = seq_id.rsplit("_", 1)
    if len(parts) == 2 and parts[1].isdigit():
        return f"{parts[0]}.{parts[1]}"
    return seq_id


def read_duplicate_members(path):
    members_by_representative = {}
    if not path or not Path(path).exists():
        return members_by_representative
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            representative = row.get("representative_id", "").strip()
            member = row.get("member_id", "").strip()
            if not representative or not member:
                continue
            members_by_representative.setdefault(representative, set()).add(member)
    return members_by_representative


def metadata_row_for_id(seq_id, metadata_by_id):
    candidates = [seq_id, clean_id(seq_id), hyphy_safe_id(clean_id(seq_id)), safe_to_metadata_id(seq_id)]
    for candidate in candidates:
        if candidate in metadata_by_id:
            return candidate, metadata_by_id[candidate]
    return "", None


def main():
    parser = argparse.ArgumentParser(description="Label tree tips as Test/Reference for HyPhy RELAX.")
    parser.add_argument("-t", "--tree", required=True, help="Input HyPhy-ready Newick tree")
    parser.add_argument("-m", "--metadata", required=True, help="Metadata TSV")
    parser.add_argument("-o", "--output", required=True, help="Output RELAX-labeled Newick tree")
    parser.add_argument("--qc", required=True, help="Output QC TSV")
    parser.add_argument("--id-column", default="accessionVersion")
    parser.add_argument("--host-column", default="hostNameScientific")
    parser.add_argument("--test-host", default="Homo sapiens")
    parser.add_argument("--test-label", default="Test")
    parser.add_argument("--reference-label", default="Reference")
    parser.add_argument("--unknown-label", default="", help="Optional label for unknown hosts; default leaves unlabeled")
    parser.add_argument("--duplicate-map", default="", help="Optional deduplication map TSV to infer host groups for representatives")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    root = NewickParser(read_newick(args.tree)).parse()
    metadata_by_id = read_host_metadata(args.metadata, args.id_column, args.host_column)
    duplicate_members = read_duplicate_members(args.duplicate_map)

    test_label = safe_label(args.test_label)
    reference_label = safe_label(args.reference_label)
    unknown_label = safe_label(args.unknown_label) if args.unknown_label else ""

    counts = {"test": 0, "reference": 0, "unknown": 0, "missing_metadata": 0}
    qc_rows = []

    for node in [item for item in iter_nodes(root) if item.is_leaf()]:
        raw_tip = node.name
        cleaned_tip = clean_id(raw_tip)
        member_ids = duplicate_members.get(cleaned_tip) or duplicate_members.get(hyphy_safe_id(cleaned_tip)) or {cleaned_tip}

        metadata_ids = []
        hosts = []
        missing_members = []
        for member_id in sorted(member_ids):
            metadata_id, metadata_row = metadata_row_for_id(member_id, metadata_by_id)
            if metadata_row is None:
                missing_members.append(member_id)
                continue
            metadata_ids.append(metadata_id)
            host = normalize_host(metadata_row.get(args.host_column, ""))
            if host:
                hosts.append(host)

        known_hosts = sorted(set(hosts))
        has_test = args.test_host in known_hosts
        has_reference = any(host != args.test_host for host in known_hosts)

        if not metadata_ids:
            status = "missing_metadata"
            host = ""
            group = ""
            label = ""
            counts["missing_metadata"] += 1
        elif has_test and has_reference:
            status = "ambiguous_mixed_duplicate_hosts"
            host = ";".join(known_hosts)
            group = "ambiguous"
            label = ""
            counts["unknown"] += 1
        else:
            host = ";".join(known_hosts)
            if not known_hosts:
                status = "unknown_host"
                group = "unknown"
                label = unknown_label
                counts["unknown"] += 1
            elif has_test:
                status = "labeled_test"
                group = "test"
                label = test_label
                counts["test"] += 1
            else:
                status = "labeled_reference"
                group = "reference"
                label = reference_label
                counts["reference"] += 1

        if label and label not in node.labels:
            node.labels.append(label)

        qc_rows.append(
            {
                "tip": raw_tip,
                "metadata_id": ",".join(metadata_ids),
                "hostNameScientific": host,
                "relax_group": group,
                "branch_label": label,
                "status": status,
                "duplicate_members": ",".join(sorted(member_ids)),
                "missing_duplicate_members": ",".join(missing_members),
            }
        )

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    Path(args.qc).parent.mkdir(parents=True, exist_ok=True)

    with open(args.output, "w") as handle:
        handle.write(write_newick(root) + ";\n")

    with open(args.qc, "w", newline="") as handle:
        fieldnames = [
            "tip",
            "metadata_id",
            "hostNameScientific",
            "relax_group",
            "branch_label",
            "status",
            "duplicate_members",
            "missing_duplicate_members",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(qc_rows)

    print(f"Saved RELAX-labeled tree: {args.output}")
    print(f"Saved QC report: {args.qc}")
    print(f"Tree tips: {len(leaf_names(root))}")
    print(f"Test tips ({args.test_host}): {counts['test']}")
    print(f"Reference tips (known non-human hosts): {counts['reference']}")
    print(f"Unknown host tips: {counts['unknown']}")
    print(f"Missing metadata tips: {counts['missing_metadata']}")

    if counts["test"] == 0 or counts["reference"] == 0:
        raise SystemExit("ERROR: RELAX needs at least one Test and one Reference branch.")
    if args.verbose:
        print(f"Root descendant tips: {count_descendant_leaves(root)}")


if __name__ == "__main__":
    main()
