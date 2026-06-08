#!/usr/bin/env python3
"""Extract HyPhy reconstructed/imputed states from HyPhy JSON files."""

from __future__ import annotations

import argparse
import csv
import json
import re
import shutil
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

GENETIC_CODE = {
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


def parse_run_name(path: Path) -> tuple[str, str, str]:
    match = re.match(r"^(?P<segment>[A-Za-z0-9]+)_(?P<group>.+)_(?P<method>FEL|MEME_imputed|MEME|aBSREL|RELAX|MSS)\.json$", path.name)
    if not match:
        return "", "", path.stem
    return match.group("segment"), match.group("group"), match.group("method")


def node_type(node: str) -> str:
    if node == "root":
        return "root"
    if node.startswith("Node") or node.startswith("internal_"):
        return "internal"
    return "tip"


@dataclass
class NewickNode:
    name: str = ""
    children: list["NewickNode"] = field(default_factory=list)


def parse_newick(text: str) -> NewickNode:
    text = text.strip().rstrip(";")
    index = 0

    def skip_ws() -> None:
        nonlocal index
        while index < len(text) and text[index].isspace():
            index += 1

    def parse_label() -> str:
        nonlocal index
        skip_ws()
        start = index
        while index < len(text) and text[index] not in ":,();":
            index += 1
        return text[start:index].strip()

    def skip_length() -> None:
        nonlocal index
        skip_ws()
        if index < len(text) and text[index] == ":":
            index += 1
            while index < len(text) and text[index] not in ",();":
                index += 1

    def parse_subtree() -> NewickNode:
        nonlocal index
        skip_ws()
        if index < len(text) and text[index] == "(":
            index += 1
            children = []
            while True:
                children.append(parse_subtree())
                skip_ws()
                if index < len(text) and text[index] == ",":
                    index += 1
                    continue
                if index < len(text) and text[index] == ")":
                    index += 1
                    break
                break
            node = NewickNode(name=parse_label(), children=children)
            skip_length()
            return node
        node = NewickNode(name=parse_label())
        skip_length()
        return node

    root = parse_subtree()
    if not root.name:
        root.name = "root"
    return root


def iter_tree_nodes(root: NewickNode):
    yield root
    for child in root.children:
        yield from iter_tree_nodes(child)


def classify_tree_node(node: NewickNode) -> str:
    if node.name == "root":
        return "root"
    if node.children:
        return "internal"
    return node_type(node.name)


def translate_codon(codon: str) -> str:
    codon = codon.upper()
    if codon == "---":
        return "-"
    if len(codon) != 3 or any(base not in "ACGT" for base in codon):
        return "X"
    return GENETIC_CODE.get(codon, "X")


def wrap(sequence: str, width: int = 80) -> str:
    return "\n".join(sequence[i : i + width] for i in range(0, len(sequence), width))


def read_json(path: Path) -> dict | None:
    try:
        with path.open() as handle:
            return json.load(handle)
    except (json.JSONDecodeError, OSError):
        return None


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, list[str]] = {}
    current = ""
    if not path.exists():
        return {}
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            current = line[1:].split()[0]
            records[current] = []
        elif current:
            records[current].append(line.strip())
    return {name: "".join(parts).upper() for name, parts in records.items()}


def clean_codon(codon: str) -> str:
    codon = str(codon).upper()
    if len(codon) != 3:
        return "NNN"
    if any(base not in "ACGT" for base in codon):
        return "NNN"
    return codon


def best_codon(states: dict | None) -> tuple[str, float | str]:
    if not isinstance(states, dict) or not states:
        return "", ""
    codon, probability = max(
        states.items(),
        key=lambda item: (float(item[1]) if isinstance(item[1], (int, float)) else -1.0, item[0]),
    )
    return clean_codon(codon), probability


def fallback_codons(sequence: str, codons: int) -> list[str]:
    return [clean_codon(sequence[index * 3 : index * 3 + 3]) for index in range(codons)]


def export_run_trees(data: dict, results_dir: Path, outdir: Path, run_id: str, segment: str, group: str) -> list[Path]:
    tree_files: list[Path] = []
    source_tree = results_dir / "inputs" / f"{segment}_{group}.hyphy_ready.treefile"
    if source_tree.exists():
        copied_tree = outdir / f"{run_id}.hyphy_ready.treefile"
        shutil.copyfile(source_tree, copied_tree)
        tree_files.append(copied_tree)

    trees_by_partition = data.get("input", {}).get("trees", {})
    if isinstance(trees_by_partition, dict):
        for partition, tree_text in sorted(trees_by_partition.items()):
            if not isinstance(tree_text, str) or not tree_text.strip():
                continue
            suffix = "hyphy_node_labeled.nwk" if len(trees_by_partition) == 1 else f"partition_{partition}.hyphy_node_labeled.nwk"
            node_tree = outdir / f"{run_id}.{suffix}"
            node_tree.write_text(tree_text.rstrip(";") + ";\n")
            tree_files.append(node_tree)
    return tree_files


def extract_file(path: Path, outdir: Path, include_tips: bool, results_dir: Path) -> dict:
    data = read_json(path)
    segment, group, method = parse_run_name(path)
    run_id = path.stem
    summary = {
        "run_id": run_id,
        "segment": segment,
        "group": group,
        "method": method,
        "json": str(path),
        "status": "ok",
        "partitions": 0,
        "rows": 0,
        "nodes": 0,
        "codons": 0,
        "codon_fasta": "",
        "amino_acid_fasta": "",
        "tree_files": "",
        "notes": "",
    }
    if data is None:
        summary.update(status="fail", notes="JSON could not be parsed")
        return summary

    if method == "MEME_imputed":
        return extract_meme_imputed_file(path, outdir, results_dir, data, summary, segment, group, method, run_id)

    substitutions = data.get("substitutions")
    if not isinstance(substitutions, dict):
        summary.update(status="missing", notes="No substitutions object in JSON")
        return summary

    outdir.mkdir(parents=True, exist_ok=True)
    long_rows: list[dict] = []
    node_states: dict[str, dict[int, str]] = defaultdict(dict)
    max_codon = int(data.get("input", {}).get("number of sites") or 0)
    partitions_seen = 0

    for partition, sites in substitutions.items():
        if not isinstance(sites, dict):
            continue
        partitions_seen += 1
        for site_key, states in sites.items():
            if not isinstance(states, dict):
                continue
            try:
                codon_position = int(site_key) + 1
            except ValueError:
                continue
            max_codon = max(max_codon, codon_position)
            for node, codon in states.items():
                kind = node_type(node)
                if kind == "tip" and not include_tips:
                    continue
                codon = str(codon).upper()
                aa = translate_codon(codon)
                node_key = f"{partition}|{node}"
                node_states[node_key][codon_position] = codon
                long_rows.append(
                    {
                        "run_id": run_id,
                        "segment": segment,
                        "group": group,
                        "method": method,
                        "partition": partition,
                        "codon_position": codon_position,
                        "node": node,
                        "node_type": kind,
                        "codon_state": codon,
                        "amino_acid_state": aa,
                    }
                )

    long_path = outdir / f"{run_id}.ancestral_states.long.tsv"
    with long_path.open("w", newline="") as handle:
        fields = [
            "run_id",
            "segment",
            "group",
            "method",
            "partition",
            "codon_position",
            "node",
            "node_type",
            "codon_state",
            "amino_acid_state",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(long_rows)

    codon_fasta = outdir / f"{run_id}.ancestral_codons.sparse.fasta"
    aa_fasta = outdir / f"{run_id}.ancestral_amino_acids.sparse.fasta"
    with codon_fasta.open("w") as codon_handle, aa_fasta.open("w") as aa_handle:
        for node_key in sorted(node_states):
            partition, node = node_key.split("|", 1)
            header = f"{run_id}|partition={partition}|node={node}|type={node_type(node)}"
            codons = [node_states[node_key].get(position, "NNN") for position in range(1, max_codon + 1)]
            aas = [translate_codon(codon) if codon != "NNN" else "X" for codon in codons]
            codon_handle.write(f">{header}\n{wrap(''.join(codons))}\n")
            aa_handle.write(f">{header}\n{wrap(''.join(aas))}\n")

    summary.update(
        partitions=partitions_seen,
        rows=len(long_rows),
        nodes=len(node_states),
        codons=max_codon,
        codon_fasta=str(codon_fasta),
        amino_acid_fasta=str(aa_fasta),
        notes="Sparse states from HyPhy substitutions map; NNN/X means not emitted for that node/site",
    )
    return summary


def extract_meme_imputed_file(
    path: Path,
    outdir: Path,
    results_dir: Path,
    data: dict,
    summary: dict,
    segment: str,
    group: str,
    method: str,
    run_id: str,
) -> dict:
    imputed_by_partition = data.get("MLE", {}).get("Imputed States")
    if not isinstance(imputed_by_partition, dict):
        summary.update(status="missing", notes="No MLE/Imputed States object in JSON")
        return summary

    outdir.mkdir(parents=True, exist_ok=True)
    max_codon = int(data.get("input", {}).get("number of sites") or 0)
    fallback_records = read_fasta(results_dir / "inputs" / f"{segment}.hyphy_ready.fasta")
    if fallback_records and max_codon <= 0:
        max_codon = max(len(sequence) // 3 for sequence in fallback_records.values())
    substitutions_by_partition = data.get("substitutions", {})
    if not isinstance(substitutions_by_partition, dict):
        substitutions_by_partition = {}
    trees_by_partition = data.get("input", {}).get("trees", {})
    if not isinstance(trees_by_partition, dict):
        trees_by_partition = {}

    long_rows: list[dict] = []
    sequence_count = 0
    partitions_seen = 0
    codon_fasta = outdir / f"{run_id}.ancestral_codons.fasta"
    aa_fasta = outdir / f"{run_id}.ancestral_amino_acids.fasta"
    tree_files = export_run_trees(data, results_dir, outdir, run_id, segment, group)

    with codon_fasta.open("w") as codon_handle, aa_fasta.open("w") as aa_handle:
        for partition, sites in sorted(imputed_by_partition.items()):
            if not isinstance(sites, dict):
                continue
            partitions_seen += 1
            partition_max = max([int(key) + 1 for key, value in sites.items() if isinstance(value, dict)] or [max_codon])
            codon_count = max(max_codon, partition_max)
            tree_text = trees_by_partition.get(partition)
            root = parse_newick(tree_text) if isinstance(tree_text, str) and tree_text.strip() else NewickNode("root")
            if root.name != "root":
                root = NewickNode("root", [root])
            tree_node_names = {node.name for node in iter_tree_nodes(root) if node.name}
            imputed_names = set(fallback_records)
            for site_records in sites.values():
                if isinstance(site_records, dict):
                    imputed_names.update(site_records)
            for name in sorted(imputed_names - tree_node_names):
                root.children.append(NewickNode(name))
                tree_node_names.add(name)

            node_codons: dict[str, list[str]] = {
                node.name: fallback_codons(fallback_records.get(node.name, ""), codon_count)
                for node in iter_tree_nodes(root)
                if node.name
            }
            for codons in node_codons.values():
                if len(codons) < codon_count:
                    codons.extend(["NNN"] * (codon_count - len(codons)))

            substitution_sites = substitutions_by_partition.get(partition, {})
            if not isinstance(substitution_sites, dict):
                substitution_sites = {}

            def state_from_imputed(site_records: dict, name: str) -> tuple[str, float | str, str, str]:
                payload = site_records.get(name)
                if not isinstance(payload, dict):
                    return "", "", "", ""
                codon, probability = best_codon(payload.get("imputed"))
                source = "imputed_tip"
                if not codon:
                    codon, probability = best_codon(payload.get("observed"))
                    source = "observed_tip"
                return codon, probability, source, payload.get("support", "")

            for codon_position in range(1, codon_count + 1):
                site_key = str(codon_position - 1)
                site_records = sites.get(site_key)
                if not isinstance(site_records, dict):
                    site_records = {}
                substitution_states = substitution_sites.get(site_key)
                if not isinstance(substitution_states, dict):
                    substitution_states = {}

                root_codon = clean_codon(substitution_states.get("root", "NNN"))
                if root_codon == "NNN":
                    root_codon, _, _, _ = state_from_imputed(site_records, "root")
                if not root_codon:
                    root_codon = "NNN"

                def fill_node(node: NewickNode, inherited_codon: str) -> None:
                    name = node.name
                    if not name:
                        return
                    explicit = substitution_states.get(name)
                    codon = clean_codon(explicit) if explicit is not None else inherited_codon
                    source = "root_state" if name == "root" else "inherited_from_parent"
                    probability: float | str = ""
                    support: float | str = ""

                    if explicit is not None:
                        source = "explicit_substitution" if name != "root" else "root_state"
                    elif not node.children:
                        imputed_codon, probability, imputed_source, support = state_from_imputed(site_records, name)
                        if imputed_codon:
                            codon = imputed_codon
                            source = imputed_source

                    node_codons[name][codon_position - 1] = codon
                    long_rows.append(
                        {
                            "run_id": run_id,
                            "segment": segment,
                            "group": group,
                            "method": method,
                            "partition": partition,
                            "codon_position": codon_position,
                            "sequence": name,
                            "node": name,
                            "node_type": classify_tree_node(node),
                            "codon_state": codon,
                            "amino_acid_state": translate_codon(codon),
                            "state_source": source,
                            "posterior_probability": probability,
                            "support": support,
                        }
                    )
                    for child in node.children:
                        fill_node(child, codon)

                fill_node(root, root_codon)

            for node in sorted(iter_tree_nodes(root), key=lambda item: (classify_tree_node(item), item.name)):
                name = node.name
                if not name or name not in node_codons:
                    continue
                codons = node_codons[name]
                header = f"{run_id}|partition={partition}|node={name}|type={classify_tree_node(node)}|source=MEME_imputed_full"
                aas = [translate_codon(codon) if codon != "NNN" else "X" for codon in codons]
                codon_handle.write(f">{header}\n{wrap(''.join(codons))}\n")
                aa_handle.write(f">{header}\n{wrap(''.join(aas))}\n")
                sequence_count += 1
                max_codon = max(max_codon, len(codons))

    long_path = outdir / f"{run_id}.imputed_states.long.tsv"
    with long_path.open("w", newline="") as handle:
        fields = [
            "run_id",
            "segment",
            "group",
            "method",
            "partition",
            "codon_position",
            "sequence",
            "node",
            "node_type",
            "codon_state",
            "amino_acid_state",
            "state_source",
            "posterior_probability",
            "support",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(long_rows)

    summary.update(
        partitions=partitions_seen,
        rows=len(long_rows),
        nodes=sequence_count,
        codons=max_codon,
        codon_fasta=str(codon_fasta),
        amino_acid_fasta=str(aa_fasta),
        tree_files=";".join(str(path) for path in tree_files),
        notes="Full tip and internal-node sequences from MEME --impute-states Yes plus HyPhy substitution-map propagation; state_source marks root, explicit substitutions, inherited states, and imputed tips",
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extract full MEME-imputed FASTA sequences or sparse substitution-map states from HyPhy JSON."
    )
    parser.add_argument("--results-dir", default="results/ANDV_trees_aln-hyphy", help="Directory containing HyPhy JSON files")
    parser.add_argument("--outdir", default=None, help="Output folder; default: <results-dir>/ancestral_sequences")
    parser.add_argument(
        "--methods",
        nargs="+",
        default=["MEME_imputed"],
        help="Methods to extract. Default: MEME_imputed. Use FEL MEME MEME_imputed aBSREL RELAX to also extract sparse substitution maps.",
    )
    parser.add_argument("--include-tips", action="store_true", help="Also include terminal/tip states")
    args = parser.parse_args()

    results_dir = Path(args.results_dir)
    outdir = Path(args.outdir) if args.outdir else results_dir / "ancestral_sequences"
    methods = set(args.methods)
    summaries = []
    outdir.mkdir(parents=True, exist_ok=True)
    for pattern in ("*.ancestral_codons.sparse.fasta", "*.ancestral_amino_acids.sparse.fasta", "*.ancestral_states.long.tsv"):
        for stale_path in outdir.glob(pattern):
            stale_path.unlink()

    for path in sorted(results_dir.rglob("*.json")):
        if "ancestral_sequences" in path.parts or "dashboard_tables" in path.parts:
            continue
        segment, group, method = parse_run_name(path)
        if method not in methods:
            continue
        summaries.append(extract_file(path, outdir, args.include_tips, results_dir))

    summary_path = outdir / "ancestral_extraction_summary.tsv"
    with summary_path.open("w", newline="") as handle:
        fields = [
            "run_id",
            "segment",
            "group",
            "method",
            "json",
            "status",
            "partitions",
            "rows",
            "nodes",
            "codons",
            "codon_fasta",
            "amino_acid_fasta",
            "tree_files",
            "notes",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(summaries)

    print(f"Wrote ancestral state exports to {outdir}")
    print(f"Summary: {summary_path}")


if __name__ == "__main__":
    main()
