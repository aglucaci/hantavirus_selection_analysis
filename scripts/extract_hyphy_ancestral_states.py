#!/usr/bin/env python3
"""Extract HyPhy reconstructed/imputed states from HyPhy JSON files."""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import defaultdict
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

    long_rows: list[dict] = []
    sequence_count = 0
    partitions_seen = 0
    codon_fasta = outdir / f"{run_id}.ancestral_codons.fasta"
    aa_fasta = outdir / f"{run_id}.ancestral_amino_acids.fasta"

    with codon_fasta.open("w") as codon_handle, aa_fasta.open("w") as aa_handle:
        for partition, sites in sorted(imputed_by_partition.items()):
            if not isinstance(sites, dict):
                continue
            partitions_seen += 1
            partition_max = max([int(key) + 1 for key, value in sites.items() if isinstance(value, dict)] or [max_codon])
            codon_count = max(max_codon, partition_max)
            taxa = set(fallback_records)
            for site_records in sites.values():
                if isinstance(site_records, dict):
                    taxa.update(site_records)

            for taxon in sorted(taxa):
                codons = fallback_codons(fallback_records.get(taxon, ""), codon_count)
                if len(codons) < codon_count:
                    codons.extend(["NNN"] * (codon_count - len(codons)))

                for site_key, site_records in sites.items():
                    if not isinstance(site_records, dict):
                        continue
                    try:
                        codon_position = int(site_key) + 1
                    except ValueError:
                        continue
                    payload = site_records.get(taxon)
                    if not isinstance(payload, dict):
                        continue
                    codon, probability = best_codon(payload.get("imputed"))
                    source = "imputed"
                    if not codon:
                        codon, probability = best_codon(payload.get("observed"))
                        source = "observed"
                    if not codon:
                        codon = codons[codon_position - 1] if codon_position <= len(codons) else "NNN"
                        probability = ""
                        source = "alignment_fallback"
                    codons[codon_position - 1] = codon
                    long_rows.append(
                        {
                            "run_id": run_id,
                            "segment": segment,
                            "group": group,
                            "method": method,
                            "partition": partition,
                            "codon_position": codon_position,
                            "sequence": taxon,
                            "codon_state": codon,
                            "amino_acid_state": translate_codon(codon),
                            "state_source": source,
                            "posterior_probability": probability,
                            "support": payload.get("support", ""),
                        }
                    )

                header = f"{run_id}|partition={partition}|sequence={taxon}|type=imputed_full"
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
        notes="Full sequences from MEME --impute-states Yes; highest posterior codon per sequence/site",
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
        fields = ["run_id", "segment", "group", "method", "json", "status", "partitions", "rows", "nodes", "codons", "codon_fasta", "amino_acid_fasta", "notes"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(summaries)

    print(f"Wrote ancestral state exports to {outdir}")
    print(f"Summary: {summary_path}")


if __name__ == "__main__":
    main()
