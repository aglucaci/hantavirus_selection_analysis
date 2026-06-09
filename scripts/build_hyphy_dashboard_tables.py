#!/usr/bin/env python3
"""Build normalized dashboard tables from the hantavirus HyPhy run."""

from __future__ import annotations

import csv
import json
import math
import os
import re
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_RESULTS = ROOT / "results" / "ANDV_trees_aln-hyphy"
RESULTS = Path(os.environ.get("HYPHY_RESULTS_DIR", DEFAULT_RESULTS))
if not RESULTS.is_absolute():
    RESULTS = ROOT / RESULTS
INPUTS = RESULTS / "inputs"
LOGS = RESULTS / "logs" if (RESULTS / "logs").exists() else ROOT / "results" / "logs"
OUT = RESULTS / "dashboard_tables"
METHODS = ("FEL", "MEME", "aBSREL", "RELAX", "MSS")
PLACEHOLDER_METHODS = ("BUSTED", "CFEL", "GARD")


def read_json(path: Path) -> dict:
    with path.open() as handle:
        return json.load(handle)


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def header_index(headers: list[list[str]], name: str) -> int | None:
    for idx, header in enumerate(headers):
        if header and header[0] == name:
            return idx
    return None


def normalize_header(value: str) -> str:
    text = re.sub(r"<[^>]+>", "", value)
    return re.sub(r"[^a-z0-9]+", "", text.lower())


def header_index_any(headers: list[list[str]], *names: str) -> int | None:
    wanted = {normalize_header(name) for name in names}
    for idx, header in enumerate(headers):
        if header and normalize_header(header[0]) in wanted:
            return idx
    return None


def as_float(value) -> float | None:
    try:
        if value in ("", None):
            return None
        return float(value)
    except (TypeError, ValueError):
        return None


def safe_neg_log10(value) -> float | None:
    number = as_float(value)
    if number is None or number <= 0:
        return None
    return -math.log10(number)


def bh_q_values(p_values: list[float | None]) -> list[float | None]:
    indexed = [(idx, value) for idx, value in enumerate(p_values) if value is not None]
    q_values: list[float | None] = [None] * len(p_values)
    if not indexed:
        return q_values
    ranked = sorted(indexed, key=lambda item: item[1], reverse=True)
    m = len(indexed)
    running = 1.0
    for reverse_rank, (idx, p_value) in enumerate(ranked, start=1):
        rank = m - reverse_rank + 1
        running = min(running, p_value * m / rank)
        q_values[idx] = min(running, 1.0)
    return q_values


def fdr_status(q_value: float | None, p_value: float | None) -> str:
    if q_value is not None and q_value <= 0.05:
        return "strict"
    if q_value is not None and q_value <= 0.1:
        return "candidate"
    if p_value is not None and p_value <= 0.1:
        return "exploratory"
    return "not_selected"


def fel_selection_call(direction: str) -> tuple[str, str]:
    if direction == "diversifying":
        return "positive_selection", "Positive selection (dN > dS)"
    if direction == "purifying":
        return "negative_selection", "Negative selection (dN < dS)"
    if direction == "neutral":
        return "neutral", "Neutral (dN = dS)"
    return "unknown", "Unknown"


def bayes_factor(prior: float | None, posterior: float | None) -> float | None:
    if prior is None or posterior is None:
        return None
    if prior >= 1:
        return 1.0 if posterior >= 1 else 0.0
    if prior <= 0:
        return 1.0
    if posterior >= 1:
        return math.inf
    if posterior <= 0:
        return 0.0
    return posterior * (1 - prior) / (1 - posterior) / prior


def parse_binary_model_key(value: str) -> list[int]:
    numbers = re.findall(r"[-+]?\d+", value)
    return [int(number) for number in numbers]


def first_model_metrics(value) -> list[float]:
    if isinstance(value, list) and value:
        first = value[0]
        if isinstance(first, list):
            return [number for number in first if isinstance(number, (int, float))]
        if all(isinstance(item, (int, float)) for item in value):
            return value
    return []


def median_number(values: list[int | float]) -> int | float | str:
    if not values:
        return ""
    sorted_values = sorted(values)
    midpoint = len(sorted_values) // 2
    if len(sorted_values) % 2:
        return sorted_values[midpoint]
    median = (sorted_values[midpoint - 1] + sorted_values[midpoint]) / 2
    return int(median) if median.is_integer() else median


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, list[str]] = {}
    name = ""
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            records[name] = []
        elif name:
            records[name].append(line.strip())
    return {key: "".join(value).upper() for key, value in records.items()}


def shannon_entropy(values: list[str]) -> float:
    if not values:
        return 0.0
    counts = Counter(values)
    total = len(values)
    return -sum((count / total) * math.log2(count / total) for count in counts.values())


def discovered_runs() -> list[tuple[str, str, str]]:
    runs = []
    for path in sorted(RESULTS.glob("*.json")):
        if path.stat().st_size == 0:
            continue
        parsed = parse_run_name(path)
        if parsed is not None:
            runs.append(parsed)
    return runs


def discovered_segments_label_sets(runs: list[tuple[str, str, str]]) -> tuple[list[str], list[str]]:
    segments = sorted({segment for segment, _, _ in runs})
    label_sets = sorted({label_set for _, label_set, _ in runs})
    for path in sorted(INPUTS.glob("*.qc.tsv")):
        name = path.name.replace(".hyphy_ready.qc.tsv", "")
        if "_" not in name:
            continue
        segment, label_set = name.split("_", 1)
        if segment not in segments:
            segments.append(segment)
        if label_set not in label_sets:
            label_sets.append(label_set)
    return sorted(segments), sorted(label_sets)


def alignment_quality_rows(segments: list[str], label_sets: list[str]) -> list[dict]:
    rows = []
    for segment in segments:
        fasta = INPUTS / f"{segment}.hyphy_ready.fasta"
        if not fasta.exists():
            continue
        records = read_fasta(fasta)
        if not records:
            continue
        sequences = list(records.values())
        codons = len(sequences[0]) // 3
        for label_set in label_sets:
            for codon_index in range(codons):
                codon_states = [seq[codon_index * 3 : codon_index * 3 + 3] for seq in sequences]
                gap_count = sum("-" in codon for codon in codon_states)
                ambiguous_count = sum(bool(re.search(r"[^ACGT-]", codon)) for codon in codon_states)
                non_gap = [codon for codon in codon_states if "-" not in codon]
                rows.append(
                    {
                        "segment": segment,
                        "label_set": label_set,
                        "codon": codon_index + 1,
                        "gap_fraction": gap_count / len(codon_states),
                        "ambiguous_fraction": ambiguous_count / len(codon_states),
                        "entropy": shannon_entropy(non_gap),
                        "non_gap_sequences": len(non_gap),
                    }
                )
    return rows


def parse_run_name(path: Path) -> tuple[str, str, str] | None:
    match = re.match(r"^(?P<segment>[^_]+)_(?P<label_set>.+)_(?P<method>FEL|MEME|aBSREL|RELAX|MSS)\.json$", path.name)
    if not match:
        return None
    return match.group("segment"), match.group("label_set"), match.group("method")


def parse_relax_branch_counts(log_path: Path) -> tuple[int | None, int | None]:
    if not log_path.exists():
        return None, None
    text = log_path.read_text(errors="replace")
    ref = re.search(r"Selected\s+(\d+)\s+branches as the _Reference_ set", text)
    test = re.search(r"Selected\s+(\d+)\s+branches as the _Test_ set", text)
    return (
        int(ref.group(1)) if ref else None,
        int(test.group(1)) if test else None,
    )


def duplicate_count(log_path: Path) -> int:
    if not log_path.exists():
        return 0
    match = re.search(r"contains\s+(\d+)\s+duplicate sequences", log_path.read_text(errors="replace"))
    return int(match.group(1)) if match else 0


def convergence_warning(log_path: Path) -> bool:
    if not log_path.exists():
        return False
    text = log_path.read_text(errors="replace").lower()
    needles = (
        "potential convergence issues",
        "multiple local maxima",
        "negative lrt",
        "unreliable k inference",
    )
    return any(needle in text for needle in needles)


def qc_rows() -> tuple[list[dict], list[dict]]:
    summaries = []
    dropped = []
    for path in sorted(INPUTS.glob("*.qc.tsv")):
        name = path.name.replace(".hyphy_ready.qc.tsv", "")
        segment, label_set = name.split("_", 1)
        rows = list(csv.DictReader(path.open(), delimiter="\t"))
        kept = [row for row in rows if row["status"] == "kept"]
        removed = [row for row in rows if row["status"] != "kept"]
        codons = ""
        if kept:
            codons = int(int(kept[0]["length_nt"]) / 3)
        summaries.append(
            {
                "segment": segment,
                "label_set": label_set,
                "input_records": len(rows),
                "kept_sequences": len(kept),
                "dropped_sequences": len(removed),
                "codons": codons,
                "status": "warn" if removed else "pass",
                "notes": "; ".join(sorted({row["reasons"] for row in removed if row["reasons"]})),
            }
        )
        for row in removed:
            dropped.append(
                {
                    "segment": segment,
                    "label_set": label_set,
                    "id": row["id"],
                    "status": row["status"],
                    "length_nt": row["length_nt"],
                    "reasons": row["reasons"],
                    "stop_codons": row["stop_codons"],
                }
            )
    return summaries, dropped


def parse_meme(path: Path, segment: str, label_set: str, method: str = "MEME") -> tuple[dict, list[dict], list[dict]]:
    data = read_json(path)
    headers = data["MLE"]["headers"]
    rows = data["MLE"]["content"]["0"]
    p_idx = header_index(headers, "p-value")
    beta_plus_idx = header_index(headers, "&beta;<sup>+</sup>")
    prop_idx = header_index(headers, "p<sup>+</sup>")
    branch_idx = header_index(headers, "# branches under selection")
    lrt_idx = header_index(headers, "LRT")
    branch_attributes = data.get("branch attributes", {}).get("0", {})
    substitutions = data.get("substitutions", {}).get("0", {})

    p_values = [as_float(row[p_idx]) if p_idx is not None else None for row in rows]
    q_values = bh_q_values(p_values)
    sites = []
    branch_rows = []
    for site_number, (row, p_value, q_value) in enumerate(zip(rows, p_values, q_values), start=1):
        if p_value is not None and p_value <= 0.1:
            site_fdr_status = fdr_status(q_value, p_value)
            sites.append(
                {
                    "segment": segment,
                    "label_set": label_set,
                    "method": method,
                    "codon": site_number,
                    "p_value": p_value,
                    "q_value": q_value,
                    "neg_log10_p": safe_neg_log10(p_value),
                    "neg_log10_q": safe_neg_log10(q_value),
                    "omega_plus": row[beta_plus_idx] if beta_plus_idx is not None else "",
                    "branch_fraction": row[prop_idx] if prop_idx is not None else "",
                    "branches_under_selection": row[branch_idx] if branch_idx is not None else "",
                    "lrt": row[lrt_idx] if lrt_idx is not None else "",
                    "fdr_status": site_fdr_status,
                }
            )
            prior = as_float(row[prop_idx]) if prop_idx is not None else None
            site_substitutions = substitutions.get(str(site_number - 1)) or {}
            for branch, attrs in branch_attributes.items():
                posterior_by_class = attrs.get("Posterior prob omega class by site", [])
                if not posterior_by_class:
                    continue
                positive_posteriors = posterior_by_class[-1]
                posterior = as_float(positive_posteriors[site_number - 1]) if site_number - 1 < len(positive_posteriors) else None
                if posterior is None:
                    continue
                codon_state = site_substitutions.get(branch, "")
                has_substitution = bool(codon_state and codon_state != "---")
                ebf = bayes_factor(prior, posterior)
                if posterior <= 0 and not has_substitution:
                    continue
                branch_rows.append(
                    {
                        "segment": segment,
                        "label_set": label_set,
                        "method": method,
                        "codon": site_number,
                        "p_value": p_value,
                        "q_value": q_value,
                        "fdr_status": site_fdr_status,
                        "branch": branch,
                        "posterior_positive_class": posterior,
                        "reconstructed_ebf": ebf,
                        "ebf_ge_100": bool(ebf is not None and ebf >= 100),
                        "branch_length": attrs.get("Global MG94xREV", ""),
                        "reconstructed_codon": codon_state,
                        "has_reconstructed_substitution": has_substitution,
                    }
                )
    selected_count = sum(1 for value in q_values if value is not None and value <= 0.1)
    summary = {
        "segment": segment,
        "label_set": label_set,
        "method": method,
        "status": "pass",
        "n_sequences": data.get("input", {}).get("number of sequences", ""),
        "codons": data.get("input", {}).get("number of sites", ""),
        "tested": len(rows),
        "significant_count": selected_count,
        "p_value": "",
        "neg_log10_p": "",
        "k": "",
        "interpretation": "foreground-branch episodic sites at BH FDR q <= 0.1",
    }
    return summary, sites, branch_rows


def parse_fel(path: Path, segment: str, label_set: str) -> tuple[dict, list[dict]]:
    data = read_json(path)
    headers = data["MLE"]["headers"]
    rows = data["MLE"]["content"]["0"]
    alpha_idx = header_index_any(headers, "alpha", "&alpha;")
    beta_idx = header_index_any(headers, "beta", "&beta;")
    p_idx = header_index_any(headers, "p-value", "p value", "p")
    lrt_idx = header_index_any(headers, "LRT", "likelihood ratio test")

    p_values = [as_float(row[p_idx]) if p_idx is not None else None for row in rows]
    q_values = bh_q_values(p_values)
    sites = []
    diversifying_count = 0
    for site_number, (row, p_value, q_value) in enumerate(zip(rows, p_values, q_values), start=1):
        if p_value is not None and p_value <= 0.1:
            alpha = as_float(row[alpha_idx]) if alpha_idx is not None else None
            beta = as_float(row[beta_idx]) if beta_idx is not None else None
            omega = beta / alpha if alpha and beta is not None else ""
            if alpha is None or beta is None:
                direction = "unknown"
            elif beta > alpha:
                direction = "diversifying"
                if q_value is not None and q_value <= 0.1:
                    diversifying_count += 1
            elif beta < alpha:
                direction = "purifying"
            else:
                direction = "neutral"
            selection_call, selection_interpretation = fel_selection_call(direction)
            sites.append(
                {
                    "segment": segment,
                    "label_set": label_set,
                    "codon": site_number,
                    "p_value": p_value,
                    "q_value": q_value,
                    "neg_log10_p": safe_neg_log10(p_value),
                    "neg_log10_q": safe_neg_log10(q_value),
                    "alpha": row[alpha_idx] if alpha_idx is not None else "",
                    "beta": row[beta_idx] if beta_idx is not None else "",
                    "omega": omega,
                    "lrt": row[lrt_idx] if lrt_idx is not None else "",
                    "direction": direction,
                    "selection_call": selection_call,
                    "selection_interpretation": selection_interpretation,
                    "status": "selected" if q_value is not None and q_value <= 0.1 else "exploratory",
                    "fdr_status": fdr_status(q_value, p_value),
                }
            )
    summary = {
        "segment": segment,
        "label_set": label_set,
        "method": "FEL",
        "status": "pass",
        "n_sequences": data.get("input", {}).get("number of sequences", ""),
        "codons": data.get("input", {}).get("number of sites", ""),
        "tested": len(rows),
        "significant_count": diversifying_count,
        "p_value": "",
        "neg_log10_p": "",
        "k": "",
        "interpretation": "foreground-branch pervasive diversifying sites at BH FDR q <= 0.1",
    }
    return summary, sites


def parse_absrel(path: Path, segment: str, label_set: str) -> tuple[dict, list[dict]]:
    data = read_json(path)
    test_results = data.get("test results", {})
    branch_records = []
    for branch, attrs in data.get("branch attributes", {}).get("0", {}).items():
        if not isinstance(attrs, dict):
            continue
        corrected = as_float(attrs.get("Corrected P-value"))
        uncorrected = as_float(attrs.get("Uncorrected P-value"))
        if corrected is None:
            continue
        branch_records.append(
            {
                "segment": segment,
                "label_set": label_set,
                "branch": branch,
                "corrected_p_value": corrected,
                "uncorrected_p_value": uncorrected,
                "q_value": "",
                "neg_log10_corrected_p": safe_neg_log10(corrected),
                "neg_log10_q": "",
                "lrt": attrs.get("LRT", ""),
                "branch_length": attrs.get("Full adaptive model", ""),
                "omega_classes": json.dumps(attrs.get("Rate Distributions", "")),
                "status": "borderline",
            }
        )
    q_values = bh_q_values([as_float(row.get("uncorrected_p_value")) for row in branch_records])
    for row, q_value in zip(branch_records, q_values):
        row["q_value"] = q_value
        row["neg_log10_q"] = safe_neg_log10(q_value)
        row["status"] = "selected" if q_value is not None and q_value <= 0.05 else "borderline"
    branches = [
        row for row in branch_records
        if row["corrected_p_value"] <= 0.1
        or row["uncorrected_p_value"] == 0
        or (row["q_value"] is not None and row["q_value"] <= 0.1)
    ]
    summary = {
        "segment": segment,
        "label_set": label_set,
        "method": "aBSREL",
        "status": "pass",
        "n_sequences": data.get("input", {}).get("number of sequences", ""),
        "codons": data.get("input", {}).get("number of sites", ""),
        "tested": test_results.get("tested", ""),
        "significant_count": len([row for row in branch_records if row["q_value"] is not None and row["q_value"] <= 0.05]),
        "p_value": "",
        "neg_log10_p": "",
        "k": "",
        "interpretation": "foreground branches selected after aBSREL correction",
    }
    return summary, branches


def parse_relax(path: Path, segment: str, label_set: str) -> tuple[dict, dict]:
    data = read_json(path)
    results = data.get("test results", {})
    p_value = as_float(results.get("p-value"))
    k_value = as_float(results.get("relaxation or intensification parameter"))
    ref_count, test_count = parse_relax_branch_counts(LOGS / f"{segment}_{label_set}_RELAX.log")
    if k_value is None:
        interpretation = "unknown"
    elif k_value > 1:
        interpretation = "intensified"
    elif k_value < 1:
        interpretation = "relaxed"
    else:
        interpretation = "neutral"
    row = {
        "segment": segment,
        "label_set": label_set,
        "test_branches": test_count,
        "reference_branches": ref_count,
        "k": k_value,
        "lrt": results.get("LRT", ""),
        "p_value": p_value,
        "q_value": "",
        "neg_log10_p": safe_neg_log10(p_value),
        "neg_log10_q": "",
        "significant": False,
        "interpretation": interpretation,
    }
    summary = {
        "segment": segment,
        "label_set": label_set,
        "method": "RELAX",
        "status": "pass",
        "n_sequences": data.get("input", {}).get("number of sequences", ""),
        "codons": data.get("input", {}).get("number of sites", ""),
        "tested": "",
        "significant_count": int(bool(p_value is not None and p_value <= 0.05)),
        "p_value": p_value,
        "q_value": "",
        "neg_log10_p": safe_neg_log10(p_value),
        "neg_log10_q": "",
        "k": k_value,
        "interpretation": interpretation,
    }
    return summary, row


def apply_relax_fdr(analysis_summary: list[dict], relax_results: list[dict]) -> None:
    q_values = bh_q_values([as_float(row.get("p_value")) for row in relax_results])
    for row, q_value in zip(relax_results, q_values):
        row["q_value"] = q_value
        row["neg_log10_q"] = safe_neg_log10(q_value)
        row["significant"] = bool(q_value is not None and q_value <= 0.05)

    relax_by_run = {
        (row["segment"], row["label_set"]): row
        for row in relax_results
    }
    for summary in analysis_summary:
        if summary.get("method") != "RELAX" or summary.get("status") != "pass":
            continue
        row = relax_by_run.get((summary.get("segment"), summary.get("label_set")))
        if not row:
            continue
        summary["q_value"] = row["q_value"]
        summary["neg_log10_q"] = row["neg_log10_q"]
        summary["significant_count"] = int(bool(row["significant"]))


def parse_mss(path: Path, segment: str, label_set: str) -> tuple[dict, dict, list[dict], list[dict]]:
    data = read_json(path)
    models = data.get("masterList", {})
    files = data.get("files", [])
    mapping = data.get("mapping", [])
    parameter_names = [
        str(item[0]).replace("mss_0.model_object.", "") if isinstance(item, list) and item else str(item)
        for item in mapping
    ]

    def first_number(value) -> float | None:
        if isinstance(value, (int, float)):
            return float(value)
        if isinstance(value, list):
            for item in value:
                number = first_number(item)
                if number is not None:
                    return number
        if isinstance(value, dict):
            for item in value.values():
                number = first_number(item)
                if number is not None:
                    return number
        return None

    scores = [number for number in (first_number(value) for value in models.values()) if number is not None]
    best_ic = min(scores) if scores else None
    model_rows = []
    parameter_counts = Counter()
    sorted_models = []
    if isinstance(models, dict):
        for model_key, model_value in models.items():
            assignments = parse_binary_model_key(model_key)
            metrics = first_model_metrics(model_value)
            ic = as_float(metrics[0]) if len(metrics) > 0 else None
            log_likelihood = as_float(metrics[1]) if len(metrics) > 1 else None
            sorted_models.append((ic if ic is not None else math.inf, model_key, assignments, metrics))
            for index, assignment in enumerate(assignments):
                if assignment and index < len(parameter_names):
                    parameter_counts[parameter_names[index]] += 1
    sorted_models.sort(key=lambda item: item[0])
    for rank, (ic, _model_key, assignments, metrics) in enumerate(sorted_models, start=1):
        active_parameters = sum(1 for assignment in assignments if assignment)
        class_rates = [number for number in metrics[2:-1] if isinstance(number, (int, float))]
        model_rows.append(
            {
                "segment": segment,
                "label_set": label_set,
                "rank": rank,
                "ic": "" if math.isinf(ic) else ic,
                "delta_ic": "" if best_ic is None or math.isinf(ic) else ic - best_ic,
                "log_likelihood": metrics[1] if len(metrics) > 1 else "",
                "active_parameters": active_parameters,
                "inactive_parameters": max(0, len(assignments) - active_parameters),
                "class_count": len(class_rates),
                "class_rates": ",".join(str(rate) for rate in class_rates),
            }
        )
    parameter_rows = [
        {
            "segment": segment,
            "label_set": label_set,
            "parameter": parameter,
            "active_models": count,
            "model_count": len(models) if isinstance(models, dict) else "",
            "active_fraction": count / len(models) if isinstance(models, dict) and models else "",
        }
        for parameter, count in sorted(parameter_counts.items(), key=lambda item: (-item[1], item[0]))
    ]
    active_counts = [row["active_parameters"] for row in model_rows]
    row = {
        "segment": segment,
        "label_set": label_set,
        "files": len(files) if isinstance(files, list) else "",
        "model_count": len(models) if isinstance(models, dict) else "",
        "best_ic": best_ic,
        "min_active_parameters": min(active_counts) if active_counts else "",
        "median_active_parameters": median_number(active_counts),
        "max_active_parameters": max(active_counts) if active_counts else "",
        "top_parameter": parameter_rows[0]["parameter"] if parameter_rows else "",
        "top_parameter_frequency": parameter_rows[0]["active_fraction"] if parameter_rows else "",
        "classes": "",
        "interpretation": "Model-inclusion frequencies for synonymous-rate heterogeneity; not a p/q-value significance test",
    }
    summary = {
        "segment": segment,
        "label_set": label_set,
        "method": "MSS",
        "status": "pass",
        "n_sequences": data.get("input", {}).get("number of sequences", ""),
        "codons": data.get("input", {}).get("number of sites", ""),
        "tested": len(models) if isinstance(models, dict) else "",
        "significant_count": "",
        "p_value": "",
        "neg_log10_p": "",
        "k": "",
        "interpretation": row["interpretation"],
    }
    return summary, row, model_rows, parameter_rows


def build_warnings(
    qc_summary: list[dict],
    dropped_rows: list[dict],
    summaries: list[dict],
    segments: list[str],
    label_sets: list[str],
) -> list[dict]:
    warnings = []
    for row in dropped_rows:
        warnings.append(
            {
                "segment": row["segment"],
                "label_set": row["label_set"],
                "method": "input_qc",
                "severity": "high" if "stop_codon" in row["reasons"] else "medium",
                "warning": f"{row['id']} dropped: {row['reasons']}",
                "suggested_action": "Inspect the original alignment record before interpreting S-segment results.",
            }
        )
    for segment in segments:
        for label_set in label_sets:
            for method in METHODS:
                log_path = LOGS / f"{segment}_{label_set}_{method}.log"
                dupes = duplicate_count(log_path)
                if dupes:
                    warnings.append(
                        {
                            "segment": segment,
                            "label_set": label_set,
                            "method": method,
                            "severity": "medium",
                            "warning": f"{dupes} duplicate sequences reported by HyPhy",
                            "suggested_action": "Consider a duplicate-aware sensitivity run if this signal becomes manuscript-critical.",
                        }
                    )
                if method == "RELAX" and convergence_warning(log_path):
                    warnings.append(
                        {
                            "segment": segment,
                            "label_set": label_set,
                            "method": method,
                            "severity": "high",
                            "warning": "RELAX convergence/local-maxima diagnostics were reported",
                            "suggested_action": "Treat K and p-values cautiously; inspect the RELAX log before biological interpretation.",
                        }
                    )
    for row in summaries:
        if row["status"] == "fail":
            warnings.append(
                {
                    "segment": row["segment"],
                    "label_set": row["label_set"],
                    "method": row["method"],
                    "severity": "high",
                    "warning": f"{row['method']} output could not be parsed",
                    "suggested_action": "Rerun this analysis; the JSON file is empty or incomplete.",
                }
            )
        if row["method"] in PLACEHOLDER_METHODS and row["status"] == "not_run":
            warnings.append(
                {
                    "segment": row["segment"],
                    "label_set": row["label_set"],
                    "method": row["method"],
                    "severity": "info",
                    "warning": f"{row['method']} output is not present in this run",
                    "suggested_action": "Add that HyPhy method to the workflow before using the corresponding dashboard tab scientifically.",
                }
            )
    return warnings


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    runs = discovered_runs()
    segments, label_sets = discovered_segments_label_sets(runs)
    analysis_summary = []
    fel_sites = []
    meme_sites = []
    meme_branch_ebf = []
    absrel_branches = []
    relax_results = []
    mss_results = []
    mss_models = []
    mss_parameter_frequency = []

    for path in sorted(RESULTS.glob("*.json")):
        if path.stat().st_size == 0:
            continue
        parsed = parse_run_name(path)
        if parsed is None:
            continue
        segment, label_set, method = parsed
        try:
            if method == "FEL":
                summary, rows = parse_fel(path, segment, label_set)
                analysis_summary.append(summary)
                fel_sites.extend(rows)
            elif method == "MEME":
                summary, rows, branch_rows = parse_meme(path, segment, label_set, method)
                analysis_summary.append(summary)
                meme_sites.extend(rows)
                meme_branch_ebf.extend(branch_rows)
            elif method == "aBSREL":
                summary, rows = parse_absrel(path, segment, label_set)
                analysis_summary.append(summary)
                absrel_branches.extend(rows)
            elif method == "RELAX":
                summary, row = parse_relax(path, segment, label_set)
                analysis_summary.append(summary)
                relax_results.append(row)
            elif method == "MSS":
                summary, row, model_rows, parameter_rows = parse_mss(path, segment, label_set)
                analysis_summary.append(summary)
                mss_results.append(row)
                mss_models.extend(model_rows)
                mss_parameter_frequency.extend(parameter_rows)
        except json.JSONDecodeError:
            analysis_summary.append(
                {
                    "segment": segment,
                    "label_set": label_set,
                    "method": method,
                    "status": "fail",
                    "n_sequences": "",
                    "codons": "",
                    "tested": "",
                    "significant_count": "",
                    "p_value": "",
                    "neg_log10_p": "",
                    "k": "",
                    "interpretation": f"could not parse {path.name}",
                }
            )

    present = {(row["segment"], row["label_set"], row["method"]) for row in analysis_summary}
    for segment in segments:
        for label_set in label_sets:
            for method in METHODS + PLACEHOLDER_METHODS:
                if (segment, label_set, method) not in present:
                    analysis_summary.append(
                        {
                            "segment": segment,
                            "label_set": label_set,
                            "method": method,
                            "status": "not_run",
                            "n_sequences": "",
                            "codons": "",
                            "tested": "",
                            "significant_count": "",
                            "p_value": "",
                            "q_value": "",
                            "neg_log10_p": "",
                            "neg_log10_q": "",
                            "k": "",
                            "interpretation": "not available in current hantavirus run",
                        }
                    )

    apply_relax_fdr(analysis_summary, relax_results)

    qc_summary, dropped_rows = qc_rows()
    quality_rows = alignment_quality_rows(segments, label_sets)
    warnings = build_warnings(qc_summary, dropped_rows, analysis_summary, segments, label_sets)

    write_tsv(
        OUT / "analysis_summary.tsv",
        sorted(analysis_summary, key=lambda row: (row["segment"], row["label_set"], row["method"])),
        ["segment", "label_set", "method", "status", "n_sequences", "codons", "tested", "significant_count", "p_value", "q_value", "neg_log10_p", "neg_log10_q", "k", "interpretation"],
    )
    write_tsv(
        OUT / "qc_summary.tsv",
        qc_summary,
        ["segment", "label_set", "input_records", "kept_sequences", "dropped_sequences", "codons", "status", "notes"],
    )
    write_tsv(
        OUT / "dropped_sequences.tsv",
        dropped_rows,
        ["segment", "label_set", "id", "status", "length_nt", "reasons", "stop_codons"],
    )
    write_tsv(
        OUT / "fel_sites.tsv",
        fel_sites,
        [
            "segment",
            "label_set",
            "codon",
            "p_value",
            "q_value",
            "neg_log10_p",
            "neg_log10_q",
            "alpha",
            "beta",
            "omega",
            "lrt",
            "direction",
            "selection_call",
            "selection_interpretation",
            "status",
            "fdr_status",
        ],
    )
    write_tsv(
        OUT / "meme_sites.tsv",
        meme_sites,
        ["segment", "label_set", "method", "codon", "p_value", "q_value", "neg_log10_p", "neg_log10_q", "omega_plus", "branch_fraction", "branches_under_selection", "lrt", "fdr_status"],
    )
    write_tsv(
        OUT / "meme_branch_ebf.tsv",
        sorted(
            meme_branch_ebf,
            key=lambda row: (
                row["segment"],
                row["label_set"],
                row["codon"],
                -(row["reconstructed_ebf"] if isinstance(row["reconstructed_ebf"], (int, float)) and math.isfinite(row["reconstructed_ebf"]) else 1e308),
                row["branch"],
            ),
        ),
        [
            "segment",
            "label_set",
            "method",
            "codon",
            "p_value",
            "q_value",
            "fdr_status",
            "branch",
            "posterior_positive_class",
            "reconstructed_ebf",
            "ebf_ge_100",
            "branch_length",
            "reconstructed_codon",
            "has_reconstructed_substitution",
        ],
    )
    write_tsv(
        OUT / "absrel_branches.tsv",
        absrel_branches,
        ["segment", "label_set", "branch", "corrected_p_value", "uncorrected_p_value", "q_value", "neg_log10_corrected_p", "neg_log10_q", "lrt", "branch_length", "omega_classes", "status"],
    )
    write_tsv(
        OUT / "relax_results.tsv",
        relax_results,
        ["segment", "label_set", "test_branches", "reference_branches", "k", "lrt", "p_value", "q_value", "neg_log10_p", "neg_log10_q", "significant", "interpretation"],
    )
    write_tsv(
        OUT / "mss_results.tsv",
        mss_results,
        [
            "segment",
            "label_set",
            "files",
            "model_count",
            "best_ic",
            "min_active_parameters",
            "median_active_parameters",
            "max_active_parameters",
            "top_parameter",
            "top_parameter_frequency",
            "classes",
            "interpretation",
        ],
    )
    write_tsv(
        OUT / "mss_models.tsv",
        mss_models,
        ["segment", "label_set", "rank", "ic", "delta_ic", "log_likelihood", "active_parameters", "inactive_parameters", "class_count", "class_rates"],
    )
    write_tsv(
        OUT / "mss_parameter_frequency.tsv",
        mss_parameter_frequency,
        ["segment", "label_set", "parameter", "active_models", "model_count", "active_fraction"],
    )
    write_tsv(
        OUT / "warnings.tsv",
        warnings,
        ["segment", "label_set", "method", "severity", "warning", "suggested_action"],
    )
    write_tsv(
        OUT / "alignment_quality_by_site.tsv",
        quality_rows,
        ["segment", "label_set", "codon", "gap_fraction", "ambiguous_fraction", "entropy", "non_gap_sequences"],
    )
    print(f"Wrote dashboard tables to {OUT}")


if __name__ == "__main__":
    main()
