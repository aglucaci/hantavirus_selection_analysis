#!/usr/bin/env python3
"""Generate workflow sample sheets from the local data directories."""

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def write_tsv(path, rows):
    path = ROOT / path
    with path.open("w", newline="") as handle:
        handle.write("segment\tgroup\talignment\ttree\tsource_tag\n")
        for row in rows:
            handle.write("\t".join(row) + "\n")
    print(f"Wrote {path.relative_to(ROOT)} ({len(rows)} rows)")


def require_file(path):
    if not (ROOT / path).is_file():
        raise SystemExit(f"ERROR: required file not found: {path}")
    return path


def hantavirus_sheet_rows():
    data_dir = Path("data/hantavirus")
    specs = [
        (
            "L",
            "ANDV",
            "sequences_L_ex_lab_cells_2026_05_27v7.cds.linsi.fasta",
            "sequences_L_ex_lab_cells_2026_05_27v7.cds.linsi.rt.annot-ANDV.nwk",
        ),
        (
            "L",
            "clade3",
            "sequences_L_ex_lab_cells_2026_05_27v7.cds.linsi.fasta",
            "sequences_L_ex_lab_cells_2026_05_27v7.cds.linsi.rt.annot-clade3.nwk",
        ),
        (
            "S",
            "ANDV",
            "sequences_S_ex_lab_cells_2026_05_27v5.cds.linsi.fasta",
            "sequences_S_ex_lab_cells_2026_05_27v5.cds.linsi.rt.annot-ANDV.nwk",
        ),
        (
            "S",
            "clade3",
            "sequences_S_ex_lab_cells_2026_05_27v5.cds.linsi.fasta",
            "sequences_S_ex_lab_cells_2026_05_27v5.cds.linsi.rt.annot-clade3.nwk",
        ),
    ]
    return [
        (
            segment,
            group,
            require_file(str(data_dir / alignment)),
            require_file(str(data_dir / tree)),
            "hantavirus_hyphy",
        )
        for segment, group, alignment, tree in specs
    ]


def andv_trees_sheet_rows():
    data_dir = Path("data/ANDV_trees_aln-hyphy")
    alignments = {
        "L": "sequences_L_ex_lab_cells_2026_05_27v7.cds.linsi.fasta",
        "M": "sequences_M_ex_lab_cells_concat_2026_05_27v6.cds.linsi.fasta",
        "S": "sequences_S_ex_lab_cells_2026_05_27v5.cds.linsi.fasta",
    }
    groups = [
        ("ANDVall", "ANDVall"),
        ("ANDVstems", "ANDVstems"),
        ("humanout", "humanout"),
    ]

    rows = []
    for segment in ("L", "M", "S"):
        alignment = require_file(str(data_dir / alignments[segment]))
        tree_prefix = alignments[segment].replace(".fasta", ".rt.hyphy")
        for group, tree_label in groups:
            tree = require_file(str(data_dir / f"{tree_prefix}-{tree_label}.nwk"))
            rows.append((segment, group, alignment, tree, "ANDV_trees_aln-hyphy"))
    return rows


def main():
    write_tsv("sample_sheet.tsv", hantavirus_sheet_rows())
    write_tsv("andv_trees_sample_sheet.tsv", andv_trees_sheet_rows())


if __name__ == "__main__":
    main()
