#!/usr/bin/env python3
"""
Label Newick tree branches by metadata groups at each group's MRCA.

Default use case:
  - Tree tip names are accessionVersion IDs.
  - Metadata has accessionVersion and hostNameScientific columns.
  - Each host group is labeled on the branch leading to its MRCA.
  - The branch label can be fixed, e.g. {Foreground}, for HyPhy-style analyses.
"""

import argparse
import csv
import re
import sys
from collections import defaultdict


class Node:
    def __init__(self):
        self.name = ""
        self.length = ""
        self.children = []
        self.parent = None
        self.labels = []

    def is_leaf(self):
        return not self.children


def clean_id(seq_id):
    seq_id = seq_id.strip().split()[0]
    if seq_id.endswith("|M"):
        seq_id = seq_id[:-2]
    return seq_id.split("|")[0]


def safe_label(value):
    value = str(value).strip()
    if not value or value.lower() in {"na", "nan", "none", "null", "unknown", "not collected", "not provided"}:
        value = "Unknown_host"
    value = re.sub(r"[^\w.-]+", "_", value)
    value = re.sub(r"_+", "_", value).strip("_")
    return value or "Unknown_host"


def hyphy_safe_id(seq_id):
    safe = re.sub(r"[^A-Za-z0-9_]+", "_", seq_id.strip())
    safe = re.sub(r"_+", "_", safe).strip("_")
    if not safe:
        safe = "sequence"
    if safe[0].isdigit():
        safe = "seq_" + safe
    return safe


def quote_name(name):
    if not name:
        return ""
    if re.search(r"[\s,:;()\[\]{}']", name):
        return "'" + name.replace("'", "''") + "'"
    return name


class NewickParser:
    def __init__(self, text):
        self.text = text.strip()
        self.index = 0

    def peek(self):
        if self.index >= len(self.text):
            return ""
        return self.text[self.index]

    def consume(self, expected=None):
        char = self.peek()
        if expected is not None and char != expected:
            raise ValueError(f"Expected '{expected}' at position {self.index}, found '{char}'")
        self.index += 1
        return char

    def skip_space(self):
        while self.peek() and self.peek().isspace():
            self.index += 1

    def parse(self):
        self.skip_space()
        root = self.parse_subtree()
        self.skip_space()
        if self.peek() == ";":
            self.consume(";")
        self.skip_space()
        if self.index != len(self.text):
            raise ValueError(f"Unexpected text at position {self.index}: {self.text[self.index:self.index + 40]}")
        return root

    def parse_subtree(self):
        self.skip_space()
        node = Node()

        if self.peek() == "(":
            self.consume("(")
            while True:
                child = self.parse_subtree()
                child.parent = node
                node.children.append(child)
                self.skip_space()
                if self.peek() == ",":
                    self.consume(",")
                    continue
                if self.peek() == ")":
                    self.consume(")")
                    break
                raise ValueError(f"Expected ',' or ')' at position {self.index}")

        node.name = self.parse_name()
        self.skip_space()
        if self.peek() == ":":
            self.consume(":")
            node.length = self.parse_length()

        return node

    def parse_name(self):
        self.skip_space()
        if self.peek() == "'":
            self.consume("'")
            pieces = []
            while True:
                char = self.consume()
                if char == "'":
                    if self.peek() == "'":
                        self.consume("'")
                        pieces.append("'")
                        continue
                    break
                pieces.append(char)
            return "".join(pieces)

        start = self.index
        while self.peek() and self.peek() not in ":,();":
            self.index += 1
        return self.text[start:self.index].strip()

    def parse_length(self):
        start = self.index
        while self.peek() and self.peek() not in ",();":
            self.index += 1
        return self.text[start:self.index].strip()


def read_newick(path):
    with open(path) as handle:
        return handle.read().strip()


def write_newick(node):
    if node.children:
        body = "(" + ",".join(write_newick(child) for child in node.children) + ")"
    else:
        body = ""

    name = quote_name(node.name)
    label = ""
    if node.labels:
        label = "{" + "+".join(node.labels) + "}"

    length = f":{node.length}" if node.length else ""
    return f"{body}{name}{label}{length}"


def iter_nodes(root):
    yield root
    for child in root.children:
        yield from iter_nodes(child)


def leaf_names(root):
    return [node.name for node in iter_nodes(root) if node.is_leaf()]


def ancestors(node):
    seen = []
    while node is not None:
        seen.append(node)
        node = node.parent
    return seen


def mrca(nodes):
    if not nodes:
        return None
    ancestor_sets = [set(ancestors(node)) for node in nodes]
    for candidate in ancestors(nodes[0]):
        if all(candidate in aset for aset in ancestor_sets[1:]):
            return candidate
    return None


def count_descendant_leaves(node):
    if node.is_leaf():
        return 1
    return sum(count_descendant_leaves(child) for child in node.children)


def read_metadata(path, id_column, label_column):
    groups = defaultdict(list)
    rows_by_id = {}

    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if id_column not in reader.fieldnames:
            raise SystemExit(f"ERROR: metadata ID column not found: {id_column}")
        if label_column not in reader.fieldnames:
            raise SystemExit(f"ERROR: metadata label column not found: {label_column}")

        for row in reader:
            raw_id = row.get(id_column, "").strip()
            if not raw_id:
                continue
            seq_id = clean_id(raw_id)
            label_raw = row.get(label_column, "")
            label = safe_label(label_raw)
            rows_by_id[seq_id] = row
            groups[label].append(seq_id)

    return groups, rows_by_id


def main():
    parser = argparse.ArgumentParser(
        description="Label Newick branches at the MRCA of metadata-defined tip groups."
    )
    parser.add_argument("-t", "--tree", required=True, help="Input Newick tree")
    parser.add_argument("-m", "--metadata", required=True, help="Metadata TSV")
    parser.add_argument("-o", "--output", required=True, help="Output labeled Newick tree")
    parser.add_argument("--qc", default="tree_mrca_labels.qc.tsv", help="Output QC TSV")
    parser.add_argument("--id-column", default="accessionVersion")
    parser.add_argument("--label-column", default="hostNameScientific")
    parser.add_argument(
        "--branch-label",
        default=None,
        help=(
            "Fixed Newick branch label to write at each selected MRCA, e.g. Foreground. "
            "If omitted, the sanitized metadata group value is used."
        ),
    )
    parser.add_argument("--min-tips", type=int, default=2, help="Minimum tree tips needed to label a group")
    parser.add_argument(
        "--skip-label",
        action="append",
        default=[],
        help="Sanitized label to skip. Can be provided multiple times.",
    )
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    if args.verbose:
        print(f"Reading tree: {args.tree}")
        print(f"Reading metadata: {args.metadata}")
        print(f"Metadata columns: {args.id_column}, {args.label_column}")

    root = NewickParser(read_newick(args.tree)).parse()
    tree_leaf_names = leaf_names(root)
    tree_nodes_by_id = {}
    for node in [n for n in iter_nodes(root) if n.is_leaf()]:
        cleaned = clean_id(node.name)
        tree_nodes_by_id[cleaned] = node
        tree_nodes_by_id[hyphy_safe_id(cleaned)] = node
    groups, _metadata_by_id = read_metadata(args.metadata, args.id_column, args.label_column)

    if args.verbose:
        print(f"Tree tips: {len(tree_leaf_names)}")
        print(f"Metadata groups: {len(groups)}")

    skip_labels = set(args.skip_label)
    qc_rows = []
    labeled_groups = 0

    for label, metadata_ids in sorted(groups.items()):
        present_nodes = []
        missing_ids = []

        for seq_id in sorted(set(metadata_ids)):
            cleaned_id = clean_id(seq_id)
            node = tree_nodes_by_id.get(cleaned_id) or tree_nodes_by_id.get(hyphy_safe_id(cleaned_id))
            if node is None:
                missing_ids.append(seq_id)
            else:
                present_nodes.append(node)

        if label in skip_labels:
            status = "skipped_by_label"
            target = None
        elif len(present_nodes) < args.min_tips:
            status = "skipped_below_min_tips"
            target = None
        else:
            target = mrca(present_nodes)
            if target is None:
                status = "no_mrca"
            else:
                branch_label = safe_label(args.branch_label) if args.branch_label else label
                if branch_label not in target.labels:
                    target.labels.append(branch_label)
                status = "labeled"
                labeled_groups += 1

        qc_rows.append(
            {
                "label": label,
                "status": status,
                "metadata_ids": len(set(metadata_ids)),
                "tree_tips_present": len(present_nodes),
                "tree_tips_missing": len(missing_ids),
                "mrca_node_name": target.name if target else "",
                "mrca_descendant_tips": count_descendant_leaves(target) if target else "",
                "missing_ids": ",".join(missing_ids[:50]),
            }
        )

        if args.verbose:
            print(
                f"{label}: {status}; present={len(present_nodes)} "
                f"missing={len(missing_ids)} "
                f"mrca_descendants={count_descendant_leaves(target) if target else ''}"
            )

    with open(args.output, "w") as handle:
        handle.write(write_newick(root) + ";\n")

    with open(args.qc, "w", newline="") as handle:
        fieldnames = [
            "label",
            "status",
            "metadata_ids",
            "tree_tips_present",
            "tree_tips_missing",
            "mrca_node_name",
            "mrca_descendant_tips",
            "missing_ids",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(qc_rows)

    print(f"Saved labeled tree: {args.output}")
    print(f"Saved QC report: {args.qc}")
    print(f"Groups labeled: {labeled_groups}")

    if labeled_groups == 0:
        print("ERROR: no groups labeled.", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
