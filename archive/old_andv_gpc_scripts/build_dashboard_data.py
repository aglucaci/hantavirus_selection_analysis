#!/usr/bin/env python3
"""
Build dashboard-ready TSV/JSON files from the ANDV GPC pipeline outputs.
"""

import argparse
import csv
import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path


STOP_CODONS = {"TAA", "TAG", "TGA"}


def read_fasta(path):
    records = []
    if not Path(path).exists():
        return records
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
                    records.append({"id": current_id, "description": current_desc, "seq": "".join(current_seq)})
                current_desc = line[1:].strip()
                current_id = current_desc.split()[0]
                current_seq = []
            else:
                current_seq.append(line)
    if current_id is not None:
        records.append({"id": current_id, "description": current_desc, "seq": "".join(current_seq)})
    return records


def read_tsv(path):
    if not Path(path).exists():
        return []
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, fieldnames):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path, data):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as handle:
        json.dump(data, handle, indent=2)


def numeric(value, default=None):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def first_existing(paths):
    for path in paths:
        if Path(path).exists():
            return Path(path)
    return None


def clean_metadata_id(seq_id):
    seq_id = (seq_id or "").strip().split()[0]
    if seq_id.endswith("|M"):
        seq_id = seq_id[:-2]
    return seq_id.split("|")[0]


def read_reference_annotation(path):
    rows = read_tsv(path)
    parsed = []
    for row in rows:
        try:
            row["_start_aa"] = int(row["start_aa"])
            row["_end_aa"] = int(row["end_aa"])
            row["_start_nt"] = int(row["start_nt"])
        except (KeyError, TypeError, ValueError):
            continue
        parsed.append(row)
    return parsed


def reference_context(site, annotation_rows, protein_seq):
    gpc_row = next((row for row in annotation_rows if row["feature"] == "GPC"), None)
    cds_start_nt = int(gpc_row["start_nt"]) if gpc_row else 52
    nt_start = cds_start_nt + ((site - 1) * 3)
    nt_end = nt_start + 2
    aa_ref = protein_seq[site - 1] if 1 <= site <= len(protein_seq) else ""

    hits = []
    priority = {"motif": 0, "protein_region": 1, "CDS": 2}
    for row in annotation_rows:
        if row["_start_aa"] <= site <= row["_end_aa"] and row["feature"] != "GPC_CDS":
            hits.append(row)
    hits.sort(key=lambda row: (priority.get(row["region_type"], 9), row["_end_aa"] - row["_start_aa"]))
    region = ";".join(row["feature"] for row in hits) if hits else ""
    return {
        "reference_codon": site,
        "reference_nt_start": nt_start,
        "reference_nt_end": nt_end,
        "aa_ref": aa_ref,
        "region": region,
    }


def load_protein(path):
    records = read_fasta(path)
    return records[0]["seq"] if records else ""


def hyphy_headers(mle):
    return [item[0] for item in mle.get("headers", [])]


def add_bh_fdr(rows, q_threshold):
    valid = [(index, row["p_value"]) for index, row in enumerate(rows) if row.get("p_value") is not None]
    valid.sort(key=lambda item: item[1])
    n = len(valid)
    previous_q = 1.0
    q_values = {}
    for rank_from_end, (index, p_value) in enumerate(reversed(valid), start=1):
        rank = n - rank_from_end + 1
        q_value = min(previous_q, p_value * n / rank)
        previous_q = q_value
        q_values[index] = min(q_value, 1.0)

    for index, row in enumerate(rows):
        q_value = q_values.get(index, "")
        row["q_value"] = q_value
        significant = q_value != "" and q_value <= q_threshold
        row["significant"] = "yes" if significant else "no"
        row["direction"] = row["_selection_direction"] if significant else "non-significant"
        row.pop("_selection_direction", None)
    return rows


def parse_fel(path, annotation_rows, protein_seq, q_threshold, significant_only=True):
    if not Path(path).exists():
        return []
    data = json.loads(Path(path).read_text())
    headers = hyphy_headers(data["MLE"])
    rows = data["MLE"]["content"].get("0", [])
    parsed = []
    for index, values in enumerate(rows, start=1):
        record = dict(zip(headers, values))
        alpha = numeric(record.get("alpha"), 0)
        beta = numeric(record.get("beta"), 0)
        p_value = numeric(record.get("p-value"), 1)
        direction = "purifying" if beta < alpha else "diversifying"
        safe_p = max(p_value if p_value is not None else 1, 1e-300)
        parsed.append(
            {
                "hyphy_site": index,
                "test": "FEL",
                "direction": "",
                "_selection_direction": direction,
                "p_value": p_value,
                "minus_log10_p": -math.log10(safe_p),
                "significant": "",
                "alpha": alpha,
                "beta": beta,
                "q_value": "",
                "omega": beta / alpha if alpha else "",
                "effect_size": beta - alpha,
                "evidence_type": "pervasive selection",
                **reference_context(index, annotation_rows, protein_seq),
            }
        )
    parsed = add_bh_fdr(parsed, q_threshold)
    if significant_only:
        parsed = [row for row in parsed if row["significant"] == "yes"]
    return parsed


def parse_meme(path, annotation_rows, protein_seq, q_threshold, significant_only=True):
    if not Path(path).exists():
        return []
    data = json.loads(Path(path).read_text())
    headers = hyphy_headers(data["MLE"])
    rows = data["MLE"]["content"].get("0", [])
    parsed = []
    for index, values in enumerate(rows, start=1):
        record = dict(zip(headers, values))
        p_value = numeric(record.get("p-value"), 1)
        alpha = numeric(record.get("&alpha;"), 0)
        beta_plus = numeric(record.get("&beta;<sup>+</sup>"), 0)
        p_plus = numeric(record.get("p<sup>+</sup>"), 0)
        branch_count = numeric(record.get("# branches under selection"), 0)
        safe_p = max(p_value if p_value is not None else 1, 1e-300)
        parsed.append(
            {
                "hyphy_site": index,
                "test": "MEME",
                "direction": "",
                "_selection_direction": "episodic diversifying",
                "p_value": p_value,
                "minus_log10_p": -math.log10(safe_p),
                "significant": "",
                "alpha": alpha,
                "beta_plus": beta_plus,
                "p_plus": p_plus,
                "branches_under_selection": branch_count,
                "q_value": "",
                "omega": beta_plus / alpha if alpha else "",
                "effect_size": beta_plus - alpha,
                "evidence_type": "episodic diversifying selection",
                **reference_context(index, annotation_rows, protein_seq),
            }
        )
    parsed = add_bh_fdr(parsed, q_threshold)
    if significant_only:
        parsed = [row for row in parsed if row["significant"] == "yes"]
    return parsed


def parse_busted(path):
    if not Path(path).exists():
        return []
    data = json.loads(Path(path).read_text())
    test = data.get("test results", {})
    p_value = test.get("p-value", "")
    return [
        {
            "test": "BUSTED",
            "lrt": test.get("LRT", ""),
            "p_value": p_value,
            "q_value": p_value,
            "evidence_type": "gene-wide episodic positive selection",
        }
    ]


def parse_absrel(path, p_threshold, significant_only=True):
    if not Path(path).exists():
        return []
    data = json.loads(Path(path).read_text())
    attrs = data.get("branch attributes", {}).get("0", {})
    rows = []
    for branch, values in attrs.items():
        corrected = values.get("Corrected P-value")
        uncorrected = values.get("Uncorrected P-value")
        if corrected is None:
            continue
        if significant_only and corrected > p_threshold:
            continue
        rows.append(
            {
                "branch": branch,
                "corrected_p_value": corrected,
                "uncorrected_p_value": uncorrected,
                "lrt": values.get("LRT", ""),
                "rate_classes": values.get("Rate classes", ""),
                "evidence_type": "episodic selection on branch",
            }
        )
    return rows


def parse_relax(path):
    if not Path(path).exists():
        return []
    data = json.loads(Path(path).read_text())
    test = data.get("test results", {})
    fits = data.get("fits", {})
    k_value = (
        test.get("relaxation or intensification parameter")
        or test.get("K")
        or test.get("Relaxation/intensification parameter")
        or ""
    )
    direction = ""
    k_numeric = numeric(k_value)
    if k_numeric is not None:
        if k_numeric > 1:
            direction = "intensified"
        elif k_numeric < 1:
            direction = "relaxed"
        else:
            direction = "unchanged"
    return [
        {
            "test": "RELAX",
            "comparison": "human-derived Test vs reservoir-derived Reference",
            "k": k_value,
            "direction": direction,
            "lrt": test.get("LRT", ""),
            "p_value": test.get("p-value", ""),
            "q_value": test.get("p-value", ""),
            "evidence_type": "selection intensity difference",
            "relax_json": str(path),
        }
    ]


def alignment_qc(alignment_path):
    records = read_fasta(alignment_path)
    rows = []
    codon_gap_counts = defaultdict(int)
    codon_ambig_counts = defaultdict(int)
    codon_stop_counts = defaultdict(int)

    for record in records:
        seq = record["seq"].upper().replace("U", "T")
        gaps = seq.count("-")
        ambiguous = sum(1 for char in seq if char not in "ACGT-")
        stops = 0
        for i in range(0, len(seq) - 2, 3):
            codon_index = (i // 3) + 1
            codon = seq[i : i + 3]
            if "-" in codon:
                codon_gap_counts[codon_index] += 1
            if any(base not in "ACGT-" for base in codon):
                codon_ambig_counts[codon_index] += 1
            if codon in STOP_CODONS:
                codon_stop_counts[codon_index] += 1
                stops += 1
        rows.append(
            {
                "id": record["id"],
                "length_nt": len(seq),
                "length_codons": len(seq) // 3 if len(seq) % 3 == 0 else "",
                "gap_fraction": gaps / len(seq) if seq else "",
                "ambiguous_fraction": ambiguous / len(seq) if seq else "",
                "internal_stop_codons": stops,
            }
        )

    n = len(records)
    site_rows = []
    max_codon = max(
        list(codon_gap_counts.keys()) + list(codon_ambig_counts.keys()) + list(codon_stop_counts.keys()) + [0]
    )
    for codon in range(1, max_codon + 1):
        site_rows.append(
            {
                "codon": codon,
                "gap_fraction": codon_gap_counts[codon] / n if n else 0,
                "ambiguous_fraction": codon_ambig_counts[codon] / n if n else 0,
                "stop_count": codon_stop_counts[codon],
            }
        )
    return rows, site_rows


def pairwise_identity_distribution(alignment_path):
    records = read_fasta(alignment_path)
    rows = []
    values = []
    for i in range(len(records)):
        seq_i = records[i]["seq"].upper()
        for j in range(i + 1, len(records)):
            seq_j = records[j]["seq"].upper()
            compared = 0
            matches = 0
            for a, b in zip(seq_i, seq_j):
                if a == "-" or b == "-":
                    continue
                if a not in "ACGT" or b not in "ACGT":
                    continue
                compared += 1
                if a == b:
                    matches += 1
            identity = matches / compared if compared else 0
            values.append(identity)
            rows.append({"seq1": records[i]["id"], "seq2": records[j]["id"], "identity": identity})
    return rows, values


def year_from_date(value):
    match = re.search(r"\d{4}", str(value or ""))
    return match.group(0) if match else "Unknown"


def metadata_summaries(metadata_path):
    rows = read_tsv(metadata_path)
    host = Counter((row.get("hostNameScientific") or "Unknown_host").strip() or "Unknown_host" for row in rows)
    country = Counter((row.get("geoLocCountry") or "Unknown").strip() or "Unknown" for row in rows)
    year = Counter(year_from_date(row.get("sampleCollectionDate")) for row in rows)
    return rows, host, country, year


def safe_to_metadata_id(seq_id):
    parts = seq_id.rsplit("_", 1)
    if len(parts) == 2 and parts[1].isdigit():
        return f"{parts[0]}.{parts[1]}"
    return seq_id


def selection_input_host_summary(hyphy_qc_path, metadata_path, duplicate_map_path):
    qc_rows = read_tsv(hyphy_qc_path)
    metadata_rows = read_tsv(metadata_path)
    metadata_by_id = {clean_metadata_id(row.get("accessionVersion", "")): row for row in metadata_rows}

    duplicate_members = defaultdict(set)
    for row in read_tsv(duplicate_map_path):
        representative = row.get("representative_id", "").strip()
        member = row.get("member_id", "").strip()
        if representative and member:
            duplicate_members[representative].add(member)

    summary = Counter()
    country_summary = Counter()
    year_summary = Counter()
    details = []
    for row in qc_rows:
        if row.get("status") != "kept":
            continue
        tip = row.get("hyphy_id") or row.get("id")
        members = duplicate_members.get(tip) or {tip}
        hosts = []
        metadata_ids = []
        countries = []
        years = []
        for member in sorted(members):
            candidates = [member, safe_to_metadata_id(member)]
            for candidate in candidates:
                metadata_row = metadata_by_id.get(clean_metadata_id(candidate))
                if metadata_row:
                    metadata_ids.append(metadata_row.get("accessionVersion", candidate))
                    host = (metadata_row.get("hostNameScientific") or "").strip()
                    if host:
                        hosts.append(host)
                    country = (metadata_row.get("geoLocCountry") or "").strip()
                    if country:
                        countries.append(country)
                    years.append(year_from_date(metadata_row.get("sampleCollectionDate")))
                    break
        unique_hosts = sorted(set(hosts))
        has_human = "Homo sapiens" in unique_hosts
        has_other = any(host != "Homo sapiens" for host in unique_hosts)
        if has_human and has_other:
            host_group = "ambiguous_mixed_human_nonhuman"
        elif has_human:
            host_group = "human"
        elif has_other:
            host_group = "known_nonhuman"
        else:
            host_group = "unknown"
        summary[host_group] += 1
        unique_countries = sorted(set(country for country in countries if country))
        unique_years = sorted(set(year for year in years if year))
        country_bucket = unique_countries[0] if len(unique_countries) == 1 else "Mixed" if unique_countries else "Unknown"
        year_bucket = unique_years[0] if len(unique_years) == 1 else "Mixed" if unique_years else "Unknown"
        country_summary[country_bucket] += 1
        year_summary[year_bucket] += 1
        details.append(
            {
                "tip": tip,
                "host_group": host_group,
                "hostNameScientific_values": ";".join(unique_hosts),
                "geoLocCountry_values": ";".join(unique_countries),
                "collection_year_values": ";".join(unique_years),
                "duplicate_member_count": len(members),
                "metadata_ids": ",".join(sorted(set(metadata_ids))),
            }
        )
    return summary, country_summary, year_summary, details


def main():
    parser = argparse.ArgumentParser(description="Build dashboard input files.")
    parser.add_argument("--outdir", default="results/dashboard")
    parser.add_argument(
        "--raw-m-fasta",
        default="data/andv_nuc-M_2026-05-11T1953.fasta/andv_nuc-M_2026-05-11T1953.fasta",
        help="Raw M-segment nucleotide FASTA used for overview counts.",
    )
    parser.add_argument(
        "--metadata",
        default="data/andv_metadata_2026-05-11T1953.tsv/andv_metadata_2026-05-11T1953.tsv",
        help="Pathoplexus metadata TSV used for host/country/year summaries.",
    )
    parser.add_argument("--fel-q", type=float, default=0.1, help="FDR q-value threshold for FEL site calls.")
    parser.add_argument("--meme-q", type=float, default=0.1, help="FDR q-value threshold for MEME site calls.")
    parser.add_argument("--absrel-p", type=float, default=0.05)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    paths = {
        "raw_m_fasta": args.raw_m_fasta,
        "metadata": args.metadata,
        "recovered_qc": "results/GPC_CDS.recovered.raw.qc.tsv",
        "recovered_fasta": "results/GPC_CDS.recovered.raw.fasta",
        "full_alignment": "results/GPC_CDS.mafft_codon_aligned.fasta",
        "dedup_alignment": "results/GPC_CDS.mafft_codon_aligned.deduplicated.fasta",
        "dedup_map": "results/GPC_CDS.mafft_codon_aligned.duplicates.tsv",
        "hyphy_ready_alignment": "results/hyphy_global_selection/ANDV_GPC_global.hyphy_ready.fasta",
        "hyphy_ready_qc": "results/hyphy_global_selection/ANDV_GPC_global.hyphy_ready.qc.tsv",
        "annotation": "data/reference/ANDV_GPC_reference_annotation.tsv",
        "reference_protein": "data/reference/NC_003467.2_GPC_protein.no_stop.fasta",
        "fel": "results/hyphy_global_selection/ANDV_GPC_global_FEL.json",
        "meme": "results/hyphy_global_selection/ANDV_GPC_global_MEME.json",
        "busted": "results/hyphy_global_selection/ANDV_GPC_global_BUSTED.json",
        "absrel": "results/hyphy_global_selection/ANDV_GPC_global_aBSREL.json",
        "relax": "results/hyphy_host_relax/ANDV_GPC_human_vs_reservoir_RELAX.json",
        "tree": "results/hyphy_global_selection/ANDV_GPC_global.hyphy_ready.treefile",
    }

    annotation_rows = read_reference_annotation(paths["annotation"])
    protein_seq = load_protein(paths["reference_protein"])
    fel_all_sites = parse_fel(paths["fel"], annotation_rows, protein_seq, args.fel_q, significant_only=False)
    meme_all_sites = parse_meme(paths["meme"], annotation_rows, protein_seq, args.meme_q, significant_only=False)
    selected_sites = [row for row in fel_all_sites + meme_all_sites if row["significant"] == "yes"]

    write_tsv(
        outdir / "ANDV_GPC_selection_annotated_sites.tsv",
        selected_sites,
        [
            "hyphy_site",
            "reference_codon",
            "reference_nt_start",
            "reference_nt_end",
            "aa_ref",
            "region",
            "test",
            "direction",
            "p_value",
            "minus_log10_p",
            "significant",
            "q_value",
            "omega",
            "effect_size",
            "evidence_type",
        ],
    )
    write_tsv(
        outdir / "ANDV_GPC_FEL.all_sites.tsv",
        fel_all_sites,
        [
            "hyphy_site",
            "reference_codon",
            "reference_nt_start",
            "reference_nt_end",
            "aa_ref",
            "region",
            "test",
            "direction",
            "p_value",
            "minus_log10_p",
            "significant",
            "q_value",
            "alpha",
            "beta",
            "omega",
            "effect_size",
            "evidence_type",
        ],
    )
    write_tsv(
        outdir / "ANDV_GPC_MEME.all_sites.tsv",
        meme_all_sites,
        [
            "hyphy_site",
            "reference_codon",
            "reference_nt_start",
            "reference_nt_end",
            "aa_ref",
            "region",
            "test",
            "direction",
            "p_value",
            "minus_log10_p",
            "significant",
            "q_value",
            "alpha",
            "beta_plus",
            "p_plus",
            "branches_under_selection",
            "omega",
            "effect_size",
            "evidence_type",
        ],
    )

    busted_summary = parse_busted(paths["busted"])
    write_tsv(
        outdir / "ANDV_GPC_BUSTED.summary.tsv",
        busted_summary,
        ["test", "lrt", "p_value", "q_value", "evidence_type"],
    )

    absrel_all_branches = parse_absrel(paths["absrel"], args.absrel_p, significant_only=False)
    absrel_branches = [row for row in absrel_all_branches if numeric(row["corrected_p_value"], 1) <= args.absrel_p]
    write_tsv(
        outdir / "ANDV_GPC_ABSREL.branches.tsv",
        absrel_branches,
        ["branch", "corrected_p_value", "uncorrected_p_value", "lrt", "rate_classes", "evidence_type"],
    )
    write_tsv(
        outdir / "ANDV_GPC_ABSREL.all_branches.tsv",
        absrel_all_branches,
        ["branch", "corrected_p_value", "uncorrected_p_value", "lrt", "rate_classes", "evidence_type"],
    )

    relax_summary = parse_relax(paths["relax"])
    write_tsv(
        outdir / "ANDV_GPC_RELAX.summary.tsv",
        relax_summary,
        ["test", "comparison", "k", "direction", "lrt", "p_value", "q_value", "evidence_type", "relax_json"],
    )

    align_rows, site_rows = alignment_qc(paths["dedup_alignment"])
    write_tsv(
        outdir / "alignment_sequence_qc.tsv",
        align_rows,
        ["id", "length_nt", "length_codons", "gap_fraction", "ambiguous_fraction", "internal_stop_codons"],
    )
    write_tsv(outdir / "alignment_site_qc.tsv", site_rows, ["codon", "gap_fraction", "ambiguous_fraction", "stop_count"])

    pairwise_rows, pairwise_values = pairwise_identity_distribution(paths["dedup_alignment"])
    write_tsv(outdir / "pairwise_identity.tsv", pairwise_rows, ["seq1", "seq2", "identity"])

    metadata_rows, host_counts, country_counts, year_counts = metadata_summaries(paths["metadata"])
    write_tsv(outdir / "metadata_host_counts.tsv", [{"host": k, "count": v} for k, v in host_counts.most_common()], ["host", "count"])
    write_tsv(outdir / "metadata_country_counts.tsv", [{"country": k, "count": v} for k, v in country_counts.most_common()], ["country", "count"])
    write_tsv(outdir / "metadata_year_counts.tsv", [{"year": k, "count": v} for k, v in sorted(year_counts.items())], ["year", "count"])

    selection_host_summary, selection_country_summary, selection_year_summary, selection_host_details = selection_input_host_summary(
        paths["hyphy_ready_qc"],
        paths["metadata"],
        paths["dedup_map"],
    )
    write_tsv(
        outdir / "selection_input_host_group_counts.tsv",
        [{"host_group": key, "count": value} for key, value in sorted(selection_host_summary.items())],
        ["host_group", "count"],
    )
    write_tsv(
        outdir / "selection_input_host_group_details.tsv",
        selection_host_details,
        [
            "tip",
            "host_group",
            "hostNameScientific_values",
            "geoLocCountry_values",
            "collection_year_values",
            "duplicate_member_count",
            "metadata_ids",
        ],
    )
    write_tsv(
        outdir / "selection_input_country_counts.tsv",
        [{"country": key, "count": value} for key, value in selection_country_summary.most_common()],
        ["country", "count"],
    )
    write_tsv(
        outdir / "selection_input_year_counts.tsv",
        [{"year": key, "count": value} for key, value in sorted(selection_year_summary.items())],
        ["year", "count"],
    )

    raw_records = read_fasta(paths["raw_m_fasta"])
    recovered_records = read_fasta(paths["recovered_fasta"])
    full_records = read_fasta(paths["full_alignment"])
    dedup_records = read_fasta(paths["dedup_alignment"])
    hyphy_records = read_fasta(paths["hyphy_ready_alignment"])
    recovered_qc = read_tsv(paths["recovered_qc"])
    hyphy_qc = read_tsv(paths["hyphy_ready_qc"])
    duplicate_map = read_tsv(paths["dedup_map"])
    duplicate_reps = {row["representative_id"] for row in duplicate_map}
    duplicate_removed = len(duplicate_map) - len(duplicate_reps) if duplicate_map else 0

    metrics = {
        "raw_m_segments": len(raw_records),
        "recovered_gpc_cds": len(recovered_records),
        "full_codon_alignment_sequences": len(full_records),
        "deduplicated_alignment_sequences": len(dedup_records),
        "hyphy_ready_sequences": len(hyphy_records),
        "hyphy_ready_human_sequences": selection_host_summary.get("human", 0),
        "hyphy_ready_known_nonhuman_sequences": selection_host_summary.get("known_nonhuman", 0),
        "hyphy_ready_unknown_host_sequences": selection_host_summary.get("unknown", 0),
        "hyphy_ready_ambiguous_host_sequences": selection_host_summary.get("ambiguous_mixed_human_nonhuman", 0),
        "duplicate_records_removed": duplicate_removed,
        "alignment_length_nt": len(dedup_records[0]["seq"]) if dedup_records else 0,
        "alignment_length_codons": len(dedup_records[0]["seq"]) // 3 if dedup_records else 0,
        "host_species_count": len(host_counts),
        "fel_sites": sum(1 for row in selected_sites if row["test"] == "FEL"),
        "meme_sites": sum(1 for row in selected_sites if row["test"] == "MEME"),
        "busted_p_value": busted_summary[0]["p_value"] if busted_summary else "",
        "absrel_positive_branches": len(absrel_branches),
        "relax_k": relax_summary[0]["k"] if relax_summary else "",
        "relax_p_value": relax_summary[0]["p_value"] if relax_summary else "",
        "relax_direction": relax_summary[0]["direction"] if relax_summary else "",
        "reference": "NC_003467.2 / CHI-7913",
        "reference_id": "PP_006W0E7_1",
        "files": paths,
        "recovery_status_counts": dict(Counter(row.get("status", "") for row in recovered_qc)),
        "hyphy_ready_status_counts": dict(Counter(row.get("status", "") for row in hyphy_qc)),
        "median_pairwise_identity": sorted(pairwise_values)[len(pairwise_values) // 2] if pairwise_values else "",
        "fel_q_threshold": args.fel_q,
        "meme_q_threshold": args.meme_q,
    }
    write_json(outdir / "overview_metrics.json", metrics)

    print(f"Saved dashboard data to: {outdir}")
    print(f"Annotated selected sites: {len(selected_sites)}")
    print(f"aBSREL positive branches: {len(absrel_branches)}")


if __name__ == "__main__":
    main()
