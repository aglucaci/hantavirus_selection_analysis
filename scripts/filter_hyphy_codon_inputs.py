#!/usr/bin/env python3
"""
Prepare codon alignment and tree for HyPhy by dropping invalid sequences.

HyPhy codon analyses require:
  - alignment length divisible by 3
  - no in-frame stop codons
  - tree tips matching the alignment

This script writes a filtered alignment, a pruned Newick tree, and a QC TSV.
"""

import argparse
import csv
import re
import sys


STOP_CODONS = {"TAA", "TAG", "TGA"}


def hyphy_safe_id(seq_id):
    safe = re.sub(r"[^A-Za-z0-9_]+", "_", seq_id.strip())
    safe = re.sub(r"_+", "_", safe).strip("_")
    if not safe:
        safe = "sequence"
    if safe[0].isdigit():
        safe = "seq_" + safe
    return safe


class Node:
    def __init__(self):
        self.name = ""
        self.labels = ""
        self.length = ""
        self.children = []
        self.parent = None

    def is_leaf(self):
        return not self.children


class NewickParser:
    def __init__(self, text):
        self.text = text.strip()
        self.index = 0

    def peek(self):
        return self.text[self.index] if self.index < len(self.text) else ""

    def consume(self, expected=None):
        char = self.peek()
        if expected is not None and char != expected:
            raise ValueError(f"Expected {expected!r} at position {self.index}, found {char!r}")
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
            raise ValueError(f"Unexpected Newick text at position {self.index}")
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
        node.labels = self.parse_labels()
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
        while self.peek() and self.peek() not in ":,();{}":
            self.index += 1
        return self.text[start:self.index].strip()

    def parse_labels(self):
        labels = []
        while self.peek() == "{":
            start = self.index
            self.consume("{")
            while self.peek() and self.peek() != "}":
                self.index += 1
            self.consume("}")
            labels.append(self.text[start:self.index])
        return "".join(labels)

    def parse_length(self):
        start = self.index
        while self.peek() and self.peek() not in ",();":
            self.index += 1
        return self.text[start:self.index].strip()


def quote_name(name):
    if not name:
        return ""
    if re.search(r"[\s,:;()\[\]{}']", name):
        return "'" + name.replace("'", "''") + "'"
    return name


def write_newick(node):
    if node.children:
        body = "(" + ",".join(write_newick(child) for child in node.children) + ")"
    else:
        body = ""
    length = f":{node.length}" if node.length else ""
    return f"{body}{quote_name(node.name)}{node.labels}{length}"


def add_branch_lengths(child_length, parent_length):
    if not parent_length:
        return child_length
    if not child_length:
        return parent_length
    try:
        return f"{float(child_length) + float(parent_length):.16g}"
    except ValueError:
        return child_length


def prune_tree(node, keep_names):
    if node.is_leaf():
        return node if node.name in keep_names else None

    kept_children = []
    for child in node.children:
        kept = prune_tree(child, keep_names)
        if kept is not None:
            kept.parent = node
            kept_children.append(kept)

    node.children = kept_children

    if not node.children:
        return None
    if len(node.children) == 1 and node.parent is not None:
        child = node.children[0]
        child.length = add_branch_lengths(child.length, node.length)
        child.parent = node.parent
        return child
    return node


def rename_tree_leaves(node, id_map):
    if node.is_leaf():
        if node.name in id_map:
            node.name = id_map[node.name]
        return

    for child in node.children:
        rename_tree_leaves(child, id_map)


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
                        {"id": current_id, "description": current_desc, "seq": "".join(current_seq)}
                    )
                current_desc = line[1:].strip()
                current_id = current_desc.split()[0]
                current_seq = []
            else:
                current_seq.append(line)
    if current_id is not None:
        records.append({"id": current_id, "description": current_desc, "seq": "".join(current_seq)})
    return records


def read_clustal(path):
    sequences = {}
    order = []
    with open(path) as handle:
        for line in handle:
            stripped = line.rstrip()
            if not stripped or stripped.upper().startswith("CLUSTAL"):
                continue
            if stripped[0].isspace():
                continue
            parts = stripped.split()
            if len(parts) < 2:
                continue
            seq_id, chunk = parts[0], parts[1]
            if set(chunk) <= {"*", ":", "."}:
                continue
            if seq_id not in sequences:
                sequences[seq_id] = []
                order.append(seq_id)
            sequences[seq_id].append(chunk)
    return [{"id": seq_id, "description": seq_id, "seq": "".join(sequences[seq_id])} for seq_id in order]


def read_alignment(path):
    with open(path) as handle:
        first = ""
        for line in handle:
            if line.strip():
                first = line.strip()
                break
    if first.upper().startswith("CLUSTAL"):
        return read_clustal(path)
    return read_fasta(path)


def write_fasta(records, path):
    with open(path, "w") as handle:
        for record in records:
            desc = f" {record['description']}" if record.get("description") else ""
            handle.write(f">{record['id']}{desc}\n")
            seq = record["seq"]
            for i in range(0, len(seq), 70):
                handle.write(seq[i : i + 70] + "\n")


def stop_codons(seq):
    hits = []
    seq = seq.upper().replace("U", "T")
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i : i + 3]
        if codon in STOP_CODONS:
            hits.append(f"{i // 3 + 1}:{codon}")
    return hits


def iter_leaf_names(node):
    if node.is_leaf():
        yield node.name
    for child in node.children:
        yield from iter_leaf_names(child)


def main():
    parser = argparse.ArgumentParser(description="Filter codon alignment and prune tree for HyPhy.")
    parser.add_argument("-a", "--alignment", required=True)
    parser.add_argument("-t", "--tree", required=True)
    parser.add_argument("--out-alignment", required=True)
    parser.add_argument("--out-tree", required=True)
    parser.add_argument("--qc", required=True)
    parser.add_argument(
        "--require-id",
        action="append",
        default=[],
        help="Require this sanitized or original ID to be present in the kept alignment and output tree.",
    )
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    records = read_alignment(args.alignment)
    kept = []
    rows = []
    id_map = {}
    safe_to_original = {}

    for record in records:
        original_id = record["id"]
        seq = record["seq"].upper().replace("U", "T")
        reasons = []
        stops = stop_codons(seq)
        if len(seq) % 3 != 0:
            reasons.append("length_not_divisible_by_3")
        if stops:
            reasons.append("stop_codon")

        status = "kept" if not reasons else "dropped"
        safe_id = hyphy_safe_id(original_id)

        if status == "kept":
            previous = safe_to_original.get(safe_id)
            if previous is not None and previous != original_id:
                raise SystemExit(
                    "ERROR: HyPhy-safe ID collision after sanitizing names: "
                    f"{previous!r} and {original_id!r} both become {safe_id!r}"
                )
            safe_to_original[safe_id] = original_id
            id_map[original_id] = safe_id

        if status == "kept":
            record["seq"] = seq
            record["description"] = ""
            record["id"] = safe_id
            kept.append(record)

        rows.append(
            {
                "id": original_id,
                "hyphy_id": safe_id if status == "kept" else "",
                "status": status,
                "length_nt": len(seq),
                "reasons": ",".join(reasons),
                "stop_codons": ",".join(stops),
            }
        )

    root = NewickParser(open(args.tree).read()).parse()
    tree_tips_before = set(iter_leaf_names(root))

    tree_tip_safe_names = {hyphy_safe_id(name) for name in tree_tips_before}
    kept = [record for record in kept if record["id"] in tree_tip_safe_names]
    rows_by_id = {row["hyphy_id"] or row["id"]: row for row in rows}

    for row in rows:
        if row["status"] == "kept" and row["hyphy_id"] not in tree_tip_safe_names:
            row["status"] = "dropped"
            row["reasons"] = ",".join(filter(None, [row["reasons"], "missing_from_tree"]))

    keep_safe_names = {record["id"] for record in kept}
    pruned = prune_tree(root, {name for name in tree_tips_before if hyphy_safe_id(name) in keep_safe_names})

    if pruned is None:
        raise SystemExit("ERROR: pruning removed all tree tips.")

    tree_tips_after_original = set(iter_leaf_names(pruned))
    tree_rename_map = {name: hyphy_safe_id(name) for name in tree_tips_after_original}
    rename_tree_leaves(pruned, tree_rename_map)
    tree_tips_after = set(iter_leaf_names(pruned))
    kept_alignment_ids = {record["id"] for record in kept}
    missing_after_rename = kept_alignment_ids - tree_tips_after
    if missing_after_rename:
        raise SystemExit(f"ERROR: sanitized alignment IDs missing from tree: {sorted(missing_after_rename)[:10]}")

    for required_id in args.require_id:
        required_safe_id = hyphy_safe_id(required_id)
        if required_safe_id not in kept_alignment_ids:
            raise SystemExit(f"ERROR: required ID is missing from kept alignment: {required_id} -> {required_safe_id}")
        if required_safe_id not in tree_tips_after:
            raise SystemExit(f"ERROR: required ID is missing from output tree: {required_id} -> {required_safe_id}")

    write_fasta(kept, args.out_alignment)
    with open(args.out_tree, "w") as handle:
        handle.write(write_newick(pruned) + ";\n")

    with open(args.qc, "w", newline="") as handle:
        fieldnames = ["id", "hyphy_id", "status", "length_nt", "reasons", "stop_codons"]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"Saved HyPhy-ready alignment: {args.out_alignment}")
    print(f"Saved HyPhy-ready tree: {args.out_tree}")
    print(f"Saved QC report: {args.qc}")
    print(f"Alignment records before: {len(records)}")
    print(f"Alignment records kept:   {len(kept)}")
    print(f"Alignment records dropped:{len(records) - len(kept)}")
    print(f"Tree tips before:         {len(tree_tips_before)}")
    print(f"Tree tips after:          {len(tree_tips_after)}")

    if args.verbose:
        for row in rows:
            if row["status"] == "dropped":
                print(f"Dropped {row['id']}: {row['reasons']} {row['stop_codons']}")


if __name__ == "__main__":
    main()
