#!/usr/bin/env python3
"""Create an MSS-GA file list from a HyPhy-ready alignment and tree."""

from __future__ import annotations

import argparse
from pathlib import Path


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    name = ""
    chunks: list[str] = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name:
                records.append((name, "".join(chunks)))
            name = line[1:].split()[0]
            chunks = []
        elif name:
            chunks.append(line.strip())
    if name:
        records.append((name, "".join(chunks)))
    return records


def write_nexus(records: list[tuple[str, str]], tree: str, path: Path) -> None:
    if not records:
        raise ValueError("alignment has no records")
    nchar = len(records[0][1])
    if any(len(sequence) != nchar for _, sequence in records):
        raise ValueError("alignment records have inconsistent lengths")
    tree = tree.strip()
    if not tree.endswith(";"):
        tree += ";"
    lines = [
        "#NEXUS",
        "BEGIN DATA;",
        f"  DIMENSIONS NTAX={len(records)} NCHAR={nchar};",
        "  FORMAT DATATYPE=DNA MISSING=? GAP=-;",
        "  MATRIX",
    ]
    lines.extend(f"  {name} {sequence}" for name, sequence in records)
    lines.extend(
        [
            "  ;",
            "END;",
            "BEGIN TREES;",
            f"  TREE tree_1 = {tree}",
            "END;",
            "",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--alignment", required=True, type=Path)
    parser.add_argument("--tree", required=True, type=Path)
    parser.add_argument("--out-nexus", required=True, type=Path)
    parser.add_argument("--out-filelist", required=True, type=Path)
    args = parser.parse_args()

    records = read_fasta(args.alignment)
    write_nexus(records, args.tree.read_text(), args.out_nexus)
    args.out_filelist.parent.mkdir(parents=True, exist_ok=True)
    args.out_filelist.write_text(f"{args.out_nexus}\n")


if __name__ == "__main__":
    main()
