import json
import math
import html
import re
import sys
from pathlib import Path

import altair as alt
import pandas as pd
import streamlit as st
import streamlit.components.v1 as components


ROOT = Path(__file__).resolve().parents[1]
DASHBOARD_DIR = ROOT / "results" / "dashboard"
RESULTS_DIR = ROOT / "results"
STRUCTURE_PDB = (
    ROOT
    / "results"
    / "colabfold"
    / "PP_006W0E7_1_97bc3"
    / "PP_006W0E7_1_97bc3_unrelaxed_rank_005_alphafold2_ptm_model_1_seed_000.pdb"
)
SCRIPTS_DIR = ROOT / "scripts"
GPC_X_DOMAIN = [1, 1138]
DOMAIN_COLOR_DOMAIN = ["Gn", "WAASA", "Gc", "TM", "Cytoplasmic tail"]
DOMAIN_COLOR_RANGE = ["#7fcdbb", "#f6c85f", "#bcb2d6", "#c49a4a", "#5b5b5b"]
FEL_COLOR_DOMAIN = ["purifying", "diversifying", "non-significant"]
FEL_COLOR_RANGE = ["#0047ff", "#ff1493", "#cfcfcf"]
FEL_LOLLIPOP_COLOR_RANGE = ["#00a6d6", "#7a00cc", "#b8b8b8"]
MEME_COLOR_DOMAIN = ["episodic diversifying", "non-significant"]
MEME_COLOR_RANGE = ["#39a800", "#e1e1e1"]
SIGNIFICANT_RING_COLOR = "#111111"


st.set_page_config(page_title="ANDV GPC Selection Dashboard", layout="wide")


def load_tsv(name):
    path = DASHBOARD_DIR / name
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t")


def load_json(name):
    path = DASHBOARD_DIR / name
    if not path.exists():
        return {}
    return json.loads(path.read_text())


def file_note(path):
    path = Path(path)
    if path.exists():
        st.caption(f"Loaded: `{path}`")
    else:
        st.warning(f"Missing file: `{path}`")


def read_fasta_lengths(path):
    path = Path(path)
    if not path.exists():
        return pd.DataFrame(columns=["id", "length_nt"])
    rows = []
    seq_id = None
    seq = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if seq_id is not None:
                rows.append({"id": seq_id, "length_nt": len("".join(seq))})
            seq_id = line[1:].split()[0]
            seq = []
        elif line.strip():
            seq.append(line.strip())
    if seq_id is not None:
        rows.append({"id": seq_id, "length_nt": len("".join(seq))})
    return pd.DataFrame(rows)


def read_pdb_ca(path):
    rows = []
    path = Path(path)
    if not path.exists():
        return pd.DataFrame(columns=["chain", "residue", "aa3", "x", "y", "z", "plddt"])
    seen = set()
    for line in path.read_text().splitlines():
        if not line.startswith("ATOM"):
            continue
        atom_name = line[12:16].strip()
        if atom_name != "CA":
            continue
        chain = line[21].strip() or "A"
        residue = int(line[22:26])
        key = (chain, residue)
        if key in seen:
            continue
        seen.add(key)
        rows.append(
            {
                "chain": chain,
                "residue": residue,
                "aa3": line[17:20].strip(),
                "x": float(line[30:38]),
                "y": float(line[38:46]),
                "z": float(line[46:54]),
                "plddt": float(line[60:66]),
            }
        )
    return pd.DataFrame(rows)


def metric_row(metrics, keys):
    cols = st.columns(len(keys))
    for col, (label, key) in zip(cols, keys):
        value = metrics.get(key, "NA")
        col.metric(label, value)


def domain_bands():
    path = ROOT / "data" / "reference" / "ANDV_GPC_reference_annotation.tsv"
    if not path.exists():
        return pd.DataFrame(columns=["feature", "display", "start_aa", "end_aa"])
    annotation = pd.read_csv(path, sep="\t")
    features = [
        "Gn",
        "WAASA_cleavage_motif",
        "Gc",
        "Gn_TM",
        "Gc_TM",
        "Gn_cytoplasmic_tail",
        "Gc_cytoplasmic_tail",
    ]
    keep = annotation[annotation["feature"].isin(features)].copy()
    if keep.empty:
        return pd.DataFrame(columns=["feature", "display", "start_aa", "end_aa"])
    keep["start_aa"] = pd.to_numeric(keep["start_aa"], errors="coerce")
    keep["end_aa"] = pd.to_numeric(keep["end_aa"], errors="coerce")
    keep = keep.dropna(subset=["start_aa", "end_aa"])
    keep["mid_aa"] = (keep["start_aa"] + keep["end_aa"]) / 2
    keep["display"] = keep["feature"].map(
        {
            "Gn": "Gn",
            "WAASA_cleavage_motif": "WAASA",
            "Gc": "Gc",
            "Gn_TM": "TM",
            "Gc_TM": "TM",
            "Gn_cytoplasmic_tail": "Cytoplasmic tail",
            "Gc_cytoplasmic_tail": "Cytoplasmic tail",
        }
    )
    domain_colors = dict(zip(DOMAIN_COLOR_DOMAIN, DOMAIN_COLOR_RANGE))
    keep["domain_color"] = keep["display"].map(domain_colors)
    return keep[["feature", "display", "start_aa", "end_aa", "mid_aa", "domain_color"]]


def protein_schematic(bands):
    if bands.empty:
        return
    chart = (
        alt.Chart(bands)
        .mark_rect(cornerRadius=3)
        .encode(
            x=alt.X("start_aa:Q", title=None, scale=alt.Scale(domain=GPC_X_DOMAIN)),
            x2="end_aa:Q",
            y=alt.value(36),
            y2=alt.value(72),
            color=alt.Color("domain_color:N", scale=None, legend=None),
            tooltip=["feature", "start_aa", "end_aa"],
        )
        .properties(height=92)
    )
    labels = (
        alt.Chart(bands)
        .mark_text(fontSize=12, fontWeight="bold", baseline="bottom")
        .encode(
            x=alt.X("mid_aa:Q", scale=alt.Scale(domain=GPC_X_DOMAIN)),
            y=alt.value(28),
            text="display:N",
            color=alt.value("#263238"),
        )
    )
    st.altair_chart(chart + labels, use_container_width=True)


def add_region_background(base_chart, bands):
    if bands.empty:
        return base_chart
    background = (
        alt.Chart(bands)
        .mark_rect(opacity=0.18)
        .encode(
            x=alt.X("start_aa:Q", scale=alt.Scale(domain=GPC_X_DOMAIN)),
            x2="end_aa:Q",
            color=alt.Color("domain_color:N", scale=None, legend=None),
        )
    )
    return background + base_chart


def fel_lollipop(fel, bands):
    if fel.empty:
        st.info("No FEL all-sites table found. Run scripts/build_dashboard_data.py after HyPhy completes.")
        return
    plot = fel.copy()
    plot["reference_codon"] = pd.to_numeric(plot["reference_codon"], errors="coerce")
    plot["effect_size"] = pd.to_numeric(plot["effect_size"], errors="coerce")
    plot["p_value"] = pd.to_numeric(plot["p_value"], errors="coerce")
    plot = plot.dropna(subset=["reference_codon", "effect_size"])
    plot["point_size"] = plot["significant"].map({"yes": 62, "no": 18}).fillna(18)
    tooltip = [
        "hyphy_site",
        "reference_codon",
        "aa_ref",
        "region",
        "direction",
        alt.Tooltip("p_value:Q", format=".3g"),
        alt.Tooltip("effect_size:Q", title="beta-alpha", format=".3g"),
        alt.Tooltip("omega:Q", format=".3g"),
    ]
    stems = (
        alt.Chart(plot)
        .mark_rule(opacity=0.45)
        .encode(
            x=alt.X("reference_codon:Q", title="Reference codon, NC_003467.2 / CHI-7913", scale=alt.Scale(domain=GPC_X_DOMAIN)),
            y=alt.Y("zero:Q", title="FEL beta-alpha"),
            y2="effect_size:Q",
            color=alt.Color("direction:N", scale=alt.Scale(domain=FEL_COLOR_DOMAIN, range=FEL_LOLLIPOP_COLOR_RANGE), legend=None),
            tooltip=tooltip,
        )
        .transform_calculate(zero="0")
    )
    points = (
        alt.Chart(plot)
        .mark_circle(size=48, opacity=0.9)
        .encode(
            x=alt.X("reference_codon:Q", scale=alt.Scale(domain=GPC_X_DOMAIN)),
            y="effect_size:Q",
            size=alt.Size("point_size:Q", legend=None),
            color=alt.Color(
                "direction:N",
                scale=alt.Scale(domain=FEL_COLOR_DOMAIN, range=FEL_LOLLIPOP_COLOR_RANGE),
                legend=alt.Legend(title="FEL call"),
            ),
            tooltip=tooltip,
        )
    )
    rings = (
        alt.Chart(plot[plot["significant"] == "yes"])
        .mark_circle(size=130, fillOpacity=0, stroke=SIGNIFICANT_RING_COLOR, strokeWidth=1.8)
        .encode(
            x=alt.X("reference_codon:Q", scale=alt.Scale(domain=GPC_X_DOMAIN)),
            y="effect_size:Q",
            tooltip=tooltip,
        )
    )
    zero = alt.Chart(pd.DataFrame({"y": [0]})).mark_rule(strokeDash=[4, 4], color="#555").encode(y="y:Q")
    chart = (stems + points + rings + zero).properties(height=360)
    st.altair_chart(add_region_background(chart, bands), use_container_width=True)


def meme_lollipop(meme, bands, threshold):
    if meme.empty:
        st.info("No MEME all-sites table found. Run scripts/build_dashboard_data.py after HyPhy completes.")
        return
    plot = meme.copy()
    plot["reference_codon"] = pd.to_numeric(plot["reference_codon"], errors="coerce")
    plot["minus_log10_p"] = pd.to_numeric(plot["minus_log10_p"], errors="coerce")
    plot["p_value"] = pd.to_numeric(plot["p_value"], errors="coerce")
    plot = plot.dropna(subset=["reference_codon", "minus_log10_p"])
    plot["point_size"] = plot["significant"].map({"yes": 62, "no": 18}).fillna(18)
    threshold_y = -math.log10(threshold)
    tooltip = [
        "hyphy_site",
        "reference_codon",
        "aa_ref",
        "region",
        "direction",
        alt.Tooltip("p_value:Q", format=".3g"),
        alt.Tooltip("minus_log10_p:Q", title="-log10 p", format=".3g"),
        alt.Tooltip("branches_under_selection:Q", format=".3g"),
    ]
    stems = (
        alt.Chart(plot)
        .mark_rule(opacity=0.45)
        .encode(
            x=alt.X("reference_codon:Q", title="Reference codon, NC_003467.2 / CHI-7913", scale=alt.Scale(domain=GPC_X_DOMAIN)),
            y=alt.Y("zero:Q", title="MEME -log10(p-value)"),
            y2="minus_log10_p:Q",
            color=alt.Color("direction:N", scale=alt.Scale(domain=MEME_COLOR_DOMAIN, range=MEME_COLOR_RANGE), legend=None),
            tooltip=tooltip,
        )
        .transform_calculate(zero="0")
    )
    points = (
        alt.Chart(plot)
        .mark_circle(size=48, opacity=0.9)
        .encode(
            x=alt.X("reference_codon:Q", scale=alt.Scale(domain=GPC_X_DOMAIN)),
            y="minus_log10_p:Q",
            size=alt.Size("point_size:Q", legend=None),
            color=alt.Color(
                "direction:N",
                scale=alt.Scale(domain=MEME_COLOR_DOMAIN, range=MEME_COLOR_RANGE),
                legend=alt.Legend(title="MEME call"),
            ),
            tooltip=tooltip,
        )
    )
    rings = (
        alt.Chart(plot[plot["significant"] == "yes"])
        .mark_circle(size=130, fillOpacity=0, stroke=SIGNIFICANT_RING_COLOR, strokeWidth=1.8)
        .encode(
            x=alt.X("reference_codon:Q", scale=alt.Scale(domain=GPC_X_DOMAIN)),
            y="minus_log10_p:Q",
            tooltip=tooltip,
        )
    )
    cutoff = (
        alt.Chart(pd.DataFrame({"threshold": [threshold_y]}))
        .mark_rule(strokeDash=[6, 4], color=FEL_LOLLIPOP_COLOR_RANGE[1])
        .encode(y="threshold:Q")
    )
    chart = (stems + points + rings + cutoff).properties(height=360)
    st.altair_chart(add_region_background(chart, bands), use_container_width=True)


def selection_overview_chart(fel, meme, bands):
    if fel.empty or meme.empty:
        st.info("Missing FEL or MEME all-sites table. Run scripts/build_dashboard_data.py after HyPhy completes.")
        return

    fel_plot = fel.copy()
    fel_plot["reference_codon"] = pd.to_numeric(fel_plot["reference_codon"], errors="coerce")
    fel_plot["effect_size"] = pd.to_numeric(fel_plot["effect_size"], errors="coerce")
    fel_plot["p_value"] = pd.to_numeric(fel_plot["p_value"], errors="coerce")
    fel_plot["q_value"] = pd.to_numeric(fel_plot["q_value"], errors="coerce")
    fel_plot["point_size"] = fel_plot["significant"].map({"yes": 62, "no": 18}).fillna(18)
    fel_plot = fel_plot.dropna(subset=["reference_codon", "effect_size"])

    meme_plot = meme.copy()
    meme_plot["reference_codon"] = pd.to_numeric(meme_plot["reference_codon"], errors="coerce")
    meme_plot["minus_log10_p"] = pd.to_numeric(meme_plot["minus_log10_p"], errors="coerce")
    meme_plot["p_value"] = pd.to_numeric(meme_plot["p_value"], errors="coerce")
    meme_plot["q_value"] = pd.to_numeric(meme_plot["q_value"], errors="coerce")
    meme_plot["point_size"] = meme_plot["significant"].map({"yes": 62, "no": 18}).fillna(18)
    meme_plot = meme_plot.dropna(subset=["reference_codon", "minus_log10_p"])

    x_axis = alt.X(
        "reference_codon:Q",
        title="Reference codon position mapped to NC_003467.2 / CHI-7913",
        scale=alt.Scale(domain=GPC_X_DOMAIN),
    )
    band_start_x = alt.X("start_aa:Q", title=None, scale=alt.Scale(domain=GPC_X_DOMAIN), axis=alt.Axis(labels=False, ticks=False))
    band_mid_x = alt.X("mid_aa:Q", title=None, scale=alt.Scale(domain=GPC_X_DOMAIN), axis=alt.Axis(labels=False, ticks=False))

    domain_track = (
        alt.Chart(bands)
        .mark_rect(cornerRadius=4)
        .encode(
            x=band_start_x,
            x2="end_aa:Q",
            y=alt.value(36),
            y2=alt.value(76),
            color=alt.Color("domain_color:N", scale=None, legend=None),
            tooltip=["feature", "start_aa", "end_aa"],
        )
        .properties(title="GPC Domain Schematic", height=105)
    )
    domain_labels = (
        alt.Chart(bands)
        .mark_text(fontSize=12, fontWeight="bold", baseline="bottom")
        .encode(x=band_mid_x, y=alt.value(28), text="display:N", color=alt.value("#263238"))
    )
    domain_boundaries = (
        alt.Chart(bands)
        .mark_rule(color="white", strokeWidth=1.5, opacity=0.85)
        .encode(x=band_start_x, y=alt.value(36), y2=alt.value(76))
    )

    fel_background = add_region_background(alt.Chart(fel_plot).mark_point(opacity=0).encode(x=x_axis), bands)
    fel_stems = (
        alt.Chart(fel_plot)
        .mark_rule(opacity=0.45)
        .encode(
            x=x_axis,
            y=alt.Y("zero:Q", title="FEL beta-alpha"),
            y2="effect_size:Q",
            color=alt.Color(
                "direction:N",
                scale=alt.Scale(domain=FEL_COLOR_DOMAIN, range=FEL_LOLLIPOP_COLOR_RANGE),
                legend=alt.Legend(title="FEL call", orient="top"),
            ),
        )
        .transform_calculate(zero="0")
    )
    fel_points = (
        alt.Chart(fel_plot)
        .mark_circle(opacity=0.9)
        .encode(
            x=x_axis,
            y="effect_size:Q",
            size=alt.Size("point_size:Q", legend=None),
            color=alt.Color(
                "direction:N",
                scale=alt.Scale(domain=FEL_COLOR_DOMAIN, range=FEL_LOLLIPOP_COLOR_RANGE),
                legend=None,
            ),
            tooltip=[
                "hyphy_site",
                "reference_codon",
                "aa_ref",
                "region",
                "direction",
                alt.Tooltip("p_value:Q", format=".3g"),
                alt.Tooltip("q_value:Q", title="FDR q", format=".3g"),
                alt.Tooltip("effect_size:Q", title="beta-alpha", format=".3g"),
                alt.Tooltip("omega:Q", format=".3g"),
            ],
        )
    )
    fel_rings = (
        alt.Chart(fel_plot[fel_plot["significant"] == "yes"])
        .mark_circle(size=140, fillOpacity=0, stroke=SIGNIFICANT_RING_COLOR, strokeWidth=1.8)
        .encode(
            x=x_axis,
            y="effect_size:Q",
            tooltip=[
                "hyphy_site",
                "reference_codon",
                "aa_ref",
                "region",
                "direction",
                alt.Tooltip("p_value:Q", format=".3g"),
                alt.Tooltip("q_value:Q", title="FDR q", format=".3g"),
                alt.Tooltip("effect_size:Q", title="beta-alpha", format=".3g"),
            ],
        )
    )
    fel_zero = alt.Chart(pd.DataFrame({"y": [0]})).mark_rule(strokeDash=[4, 4], color="#555").encode(y="y:Q")
    fel_track = (fel_background + fel_stems + fel_points + fel_rings + fel_zero).properties(
        title="FEL: Pervasive Site-Level Selection",
        height=260,
    )

    meme_background = add_region_background(alt.Chart(meme_plot).mark_point(opacity=0).encode(x=x_axis), bands)
    meme_stems = (
        alt.Chart(meme_plot)
        .mark_rule(opacity=0.45)
        .encode(
            x=x_axis,
            y=alt.Y("zero:Q", title="MEME -log10(p-value)"),
            y2="minus_log10_p:Q",
            color=alt.Color(
                "direction:N",
                scale=alt.Scale(domain=MEME_COLOR_DOMAIN, range=MEME_COLOR_RANGE),
                legend=alt.Legend(title="MEME call", orient="top"),
            ),
        )
        .transform_calculate(zero="0")
    )
    meme_points = (
        alt.Chart(meme_plot)
        .mark_circle(opacity=0.9)
        .encode(
            x=x_axis,
            y="minus_log10_p:Q",
            size=alt.Size("point_size:Q", legend=None),
            color=alt.Color(
                "direction:N",
                scale=alt.Scale(domain=MEME_COLOR_DOMAIN, range=MEME_COLOR_RANGE),
                legend=None,
            ),
            tooltip=[
                "hyphy_site",
                "reference_codon",
                "aa_ref",
                "region",
                "direction",
                alt.Tooltip("p_value:Q", format=".3g"),
                alt.Tooltip("q_value:Q", title="FDR q", format=".3g"),
                alt.Tooltip("minus_log10_p:Q", title="-log10 p", format=".3g"),
                alt.Tooltip("branches_under_selection:Q", format=".3g"),
            ],
        )
    )
    meme_rings = (
        alt.Chart(meme_plot[meme_plot["significant"] == "yes"])
        .mark_circle(size=140, fillOpacity=0, stroke=SIGNIFICANT_RING_COLOR, strokeWidth=1.8)
        .encode(
            x=x_axis,
            y="minus_log10_p:Q",
            tooltip=[
                "hyphy_site",
                "reference_codon",
                "aa_ref",
                "region",
                "direction",
                alt.Tooltip("p_value:Q", format=".3g"),
                alt.Tooltip("q_value:Q", title="FDR q", format=".3g"),
                alt.Tooltip("minus_log10_p:Q", title="-log10 p", format=".3g"),
            ],
        )
    )
    meme_track = (meme_background + meme_stems + meme_points + meme_rings).properties(
        title="MEME: Episodic Diversifying Selection",
        height=260,
    )

    chart = (
        alt.vconcat(domain_track + domain_boundaries + domain_labels, fel_track, meme_track, spacing=12)
        .resolve_scale(x="shared", color="independent")
        .properties(
            title=alt.TitleParams(
                "Andes Virus GPC Selection Overview",
                subtitle="FEL and MEME site-level results mapped to NC_003467.2 / CHI-7913",
            )
        )
        .configure_view(stroke=None)
    )
    st.altair_chart(chart, use_container_width=True)


def selection_takeaways(fel_sig, meme_sig):
    purifying = int((fel_sig["direction"] == "purifying").sum()) if not fel_sig.empty else 0
    diversifying = int((fel_sig["direction"] == "diversifying").sum()) if not fel_sig.empty else 0
    episodic = len(meme_sig) if not meme_sig.empty else 0
    meme_region = "exposed Gn/Gc regions"
    if not meme_sig.empty and "region" in meme_sig.columns:
        region_counts = meme_sig["region"].fillna("Unannotated").value_counts()
        if not region_counts.empty:
            meme_region = region_counts.index[0].replace(";", " / ")

    st.info(
        "\n".join(
            [
                "Key Takeaways",
                f"1. FEL identifies {purifying} pervasively purifying sites and {diversifying} pervasively diversifying sites.",
                f"2. MEME identifies {episodic} sites with episodic diversifying selection.",
                f"3. The most common MEME annotation in this run is {meme_region}.",
                "4. All sites are mapped to NC_003467.2 / CHI-7913 reference codons for region-aware interpretation.",
            ]
        )
    )


def page_overview(metrics):
    st.title("ANDV GPC Selection Dashboard")
    st.write(
        "This dashboard summarizes codon-level selection analyses of Andes virus GPC "
        "sequences mapped to the NC_003467.2 / CHI-7913 M-segment reference. Sites "
        "are annotated by Gn/Gc region, host metadata, and HyPhy selection test."
    )

    metric_row(
        metrics,
        [
            ("Raw M segments", "raw_m_segments"),
            ("Recovered GPC CDS", "recovered_gpc_cds"),
            ("Deduplicated seqs", "deduplicated_alignment_sequences"),
            ("HyPhy-ready seqs", "hyphy_ready_sequences"),
        ],
    )
    metric_row(
        metrics,
        [
            ("Alignment codons", "alignment_length_codons"),
            ("FEL sites", "fel_sites"),
            ("MEME sites", "meme_sites"),
            ("aBSREL branches", "absrel_positive_branches"),
        ],
    )
    metric_row(
        metrics,
        [
            ("Host species", "host_species_count"),
            ("Duplicate records removed", "duplicate_records_removed"),
            ("Median pairwise identity", "median_pairwise_identity"),
            ("Reference", "reference"),
        ],
    )

    busted_p = metrics.get("busted_p_value", "")
    if busted_p != "":
        st.info(f"BUSTED gene-wide p-value / FDR q-value: `{busted_p}`")

    st.subheader("Selection-Analysis Input Host Groups")
    host_input = load_tsv("selection_input_host_group_counts.tsv")
    if not host_input.empty:
        cols = st.columns(4)
        values = dict(zip(host_input["host_group"], host_input["count"]))
        cols[0].metric("Human", values.get("human", 0))
        cols[1].metric("Known non-human", values.get("known_nonhuman", 0))
        cols[2].metric("Unknown host", values.get("unknown", 0))
        cols[3].metric("Ambiguous", values.get("ambiguous_mixed_human_nonhuman", 0))
        st.bar_chart(host_input.set_index("host_group"))

    st.subheader("Pipeline Status")
    status = pd.DataFrame(
        [
            {"step": "Raw M segments", "count": metrics.get("raw_m_segments", 0)},
            {"step": "Recovered GPC CDS", "count": metrics.get("recovered_gpc_cds", 0)},
            {"step": "MAFFT codon alignment", "count": metrics.get("full_codon_alignment_sequences", 0)},
            {"step": "Deduplicated alignment", "count": metrics.get("deduplicated_alignment_sequences", 0)},
            {"step": "HyPhy-ready alignment", "count": metrics.get("hyphy_ready_sequences", 0)},
        ]
    )
    st.bar_chart(status.set_index("step"))


def page_qc(metrics):
    st.header("Dataset QC")

    lengths = read_fasta_lengths(ROOT / "results" / "GPC_CDS.recovered.raw.fasta")
    if not lengths.empty:
        st.subheader("CDS Length Histogram")
        st.caption("Expected full GPC CDS is near 3414 nt without terminal stop; 2880 nt cutoff shown in table context.")
        st.bar_chart(lengths["length_nt"].value_counts().sort_index())
        st.dataframe(lengths.describe(), use_container_width=True)

    host_counts = load_tsv("metadata_host_counts.tsv")
    country_counts = load_tsv("metadata_country_counts.tsv")
    year_counts = load_tsv("metadata_year_counts.tsv")

    left, right = st.columns(2)
    with left:
        st.subheader("Host Species Counts")
        if not host_counts.empty:
            st.bar_chart(host_counts.head(25).set_index("host"))
            st.dataframe(host_counts, use_container_width=True)
    with right:
        st.subheader("Country Counts")
        if not country_counts.empty:
            st.bar_chart(country_counts.head(25).set_index("country"))
            st.dataframe(country_counts, use_container_width=True)

    st.subheader("Collection Year Counts")
    if not year_counts.empty:
        st.bar_chart(year_counts.set_index("year"))

    st.subheader("Dropped / Filtered Sequences")
    hyphy_qc = pd.read_csv(ROOT / "results" / "hyphy_global_selection" / "ANDV_GPC_global.hyphy_ready.qc.tsv", sep="\t")
    dropped = hyphy_qc[hyphy_qc["status"] != "kept"]
    st.dataframe(dropped, use_container_width=True)

    st.subheader("Recovery Status Counts")
    st.json(metrics.get("recovery_status_counts", {}))


def page_alignment():
    st.header("Alignment QC")
    seq_qc = load_tsv("alignment_sequence_qc.tsv")
    site_qc = load_tsv("alignment_site_qc.tsv")
    pairwise = load_tsv("pairwise_identity.tsv")

    if not seq_qc.empty:
        metric_row(
            {
                "n": len(seq_qc),
                "codons": int(seq_qc["length_codons"].dropna().iloc[0]) if not seq_qc.empty else 0,
                "mean_gap": round(seq_qc["gap_fraction"].mean(), 6),
                "stops": int(seq_qc["internal_stop_codons"].sum()),
            },
            [
                ("Sequences", "n"),
                ("Codons", "codons"),
                ("Mean gap fraction", "mean_gap"),
                ("Internal stops", "stops"),
            ],
        )

        st.subheader("Gap Fraction by Sequence")
        st.bar_chart(seq_qc.set_index("id")["gap_fraction"])
        st.dataframe(seq_qc.sort_values("gap_fraction", ascending=False), use_container_width=True)

    if not pairwise.empty:
        st.subheader("Pairwise Identity Distribution")
        bins = pd.cut(pairwise["identity"], bins=30)
        hist = pairwise.groupby(bins, observed=False).size().reset_index(name="count")
        hist["identity_bin"] = hist["identity"].astype(str)
        st.bar_chart(hist.set_index("identity_bin")["count"])
        st.dataframe(pairwise["identity"].describe(), use_container_width=True)

    if not site_qc.empty:
        st.subheader("Gap Fraction by Codon")
        st.line_chart(site_qc.set_index("codon")[["gap_fraction", "ambiguous_fraction"]])
        st.dataframe(site_qc.sort_values("gap_fraction", ascending=False).head(50), use_container_width=True)


def page_selection():
    st.header("Selection Sites")
    sites = load_tsv("ANDV_GPC_selection_annotated_sites.tsv")
    fel_all = load_tsv("ANDV_GPC_FEL.all_sites.tsv")
    meme_all = load_tsv("ANDV_GPC_MEME.all_sites.tsv")
    busted = load_tsv("ANDV_GPC_BUSTED.summary.tsv")
    bands = domain_bands()

    if not busted.empty:
        st.subheader("BUSTED Summary")
        st.dataframe(busted, use_container_width=True)

    if sites.empty:
        st.info("No selected-site table found. Run scripts/build_dashboard_data.py after HyPhy completes.")
        return

    fel_sig = fel_all[fel_all["significant"] == "yes"] if not fel_all.empty else sites[sites["test"] == "FEL"]
    meme_sig = meme_all[meme_all["significant"] == "yes"] if not meme_all.empty else sites[sites["test"] == "MEME"]

    st.subheader("Andes Virus GPC Selection Overview")
    st.caption("FEL and MEME site-level results mapped to NC_003467.2 / CHI-7913 reference codons.")

    with st.container(border=True):
        metric_row(
            {
                "purifying": int((fel_sig["direction"] == "purifying").sum()) if not fel_sig.empty else 0,
                "diversifying": int((fel_sig["direction"] == "diversifying").sum()) if not fel_sig.empty else 0,
                "episodic": int((meme_sig["direction"] == "episodic diversifying").sum()) if not meme_sig.empty else 0,
                "threshold": f"FEL {metrics.get('fel_q_threshold', 0.1)} / MEME {metrics.get('meme_q_threshold', 0.1)}",
            },
            [
                ("Purifying sites", "purifying"),
                ("Diversifying sites", "diversifying"),
                ("MEME episodic sites", "episodic"),
                ("FDR q thresholds", "threshold"),
            ],
        )
        selection_overview_chart(fel_all, meme_all, bands)

        selection_takeaways(fel_sig, meme_sig)

    st.subheader("Top FEL Sites")
    if not fel_sig.empty:
        st.dataframe(
            fel_sig.sort_values("q_value").head(30)[
                ["reference_codon", "aa_ref", "region", "direction", "p_value", "q_value", "effect_size", "omega"]
            ],
            use_container_width=True,
        )

    st.subheader("Top MEME Sites")
    if not meme_sig.empty:
        st.dataframe(
            meme_sig.sort_values("q_value").head(30)[
                [
                    "reference_codon",
                    "aa_ref",
                    "region",
                    "p_value",
                    "q_value",
                    "minus_log10_p",
                    "branches_under_selection",
                    "omega",
                ]
            ],
            use_container_width=True,
        )

    st.subheader("Selected-Site Counts by Region")
    region_counts = sites.assign(region=sites["region"].fillna("Unannotated")).groupby(["region", "test"]).size()
    if not region_counts.empty:
        st.bar_chart(region_counts.unstack(fill_value=0))

    st.subheader("Selected-Site Heatmap Table")
    heat = sites.pivot_table(index="reference_codon", columns="test", values="q_value", aggfunc="min")
    st.dataframe(heat, use_container_width=True)

    st.subheader("All Significant Site-Test Rows")
    st.dataframe(sites.sort_values(["reference_codon", "test"]), use_container_width=True)


def extract_relax_label(name):
    match = re.search(r"\{([^}]+)\}", str(name))
    label = match.group(1) if match else ""
    clean_name = re.sub(r"\{[^}]+\}", "", str(name))
    return clean_name, label


def relax_tree_segments(tree_path, qc):
    if not tree_path.exists():
        return pd.DataFrame(), pd.DataFrame()
    if str(SCRIPTS_DIR) not in sys.path:
        sys.path.insert(0, str(SCRIPTS_DIR))
    from label_tree_mrca_by_metadata import NewickParser

    root = NewickParser(tree_path.read_text()).parse()
    host_by_tip = {}
    group_by_tip = {}
    if not qc.empty:
        for row in qc.to_dict("records"):
            tip = str(row.get("tip", ""))
            host_by_tip[tip] = row.get("hostNameScientific", "")
            group_by_tip[tip] = row.get("relax_group", "")

    y_positions = {}
    counter = {"y": 0}

    def assign_y(node):
        if node.is_leaf():
            clean_name, _label = extract_relax_label(node.name)
            y_positions[id(node)] = counter["y"]
            counter["y"] += 1
            return y_positions[id(node)]
        child_ys = [assign_y(child) for child in node.children]
        y_positions[id(node)] = sum(child_ys) / len(child_ys)
        return y_positions[id(node)]

    assign_y(root)

    x_positions = {id(root): 0.0}

    def assign_x(node):
        parent_x = x_positions[id(node)]
        for child in node.children:
            try:
                length = float(child.length) if child.length else 0.0
            except ValueError:
                length = 0.0
            x_positions[id(child)] = parent_x + length
            assign_x(child)

    assign_x(root)

    segments = []
    tips = []

    def group_for_node(node):
        clean_name, label = extract_relax_label(node.name)
        if label == "Test":
            return "test"
        if label == "Reference":
            return "reference"
        return group_by_tip.get(clean_name, "") or "unlabeled"

    def visit(node):
        node_x = x_positions[id(node)]
        node_y = y_positions[id(node)]
        if node.children:
            child_ys = [y_positions[id(child)] for child in node.children]
            segments.append(
                {
                    "x": node_x,
                    "x2": node_x,
                    "y": min(child_ys),
                    "y2": max(child_ys),
                    "group": "internal",
                    "tip": "",
                    "host": "",
                    "branch_length": "",
                    "segment": "vertical",
                }
            )
            for child in node.children:
                clean_name, _label = extract_relax_label(child.name)
                try:
                    branch_length = float(child.length) if child.length else 0.0
                except ValueError:
                    branch_length = 0.0
                group = group_for_node(child)
                segments.append(
                    {
                        "x": node_x,
                        "x2": x_positions[id(child)],
                        "y": y_positions[id(child)],
                        "y2": y_positions[id(child)],
                        "group": group,
                        "tip": clean_name,
                        "host": host_by_tip.get(clean_name, ""),
                        "branch_length": branch_length,
                        "segment": "branch",
                    }
                )
                visit(child)
        else:
            clean_name, _label = extract_relax_label(node.name)
            tips.append(
                {
                    "x": node_x,
                    "y": node_y,
                    "tip": clean_name,
                    "host": host_by_tip.get(clean_name, ""),
                    "group": group_for_node(node),
                }
            )

    visit(root)
    return pd.DataFrame(segments), pd.DataFrame(tips)


def host_details_lookup():
    details = load_tsv("selection_input_host_group_details.tsv")
    host_by_tip = {}
    group_by_tip = {}
    year_by_tip = {}
    if details.empty:
        return host_by_tip, group_by_tip, year_by_tip
    for row in details.to_dict("records"):
        tip = str(row.get("tip", ""))
        host_by_tip[tip] = row.get("hostNameScientific_values", "")
        group_by_tip[tip] = row.get("host_group", "")
        year_by_tip[tip] = row.get("collection_year_values", "")
    return host_by_tip, group_by_tip, year_by_tip


def generic_tree_segments(tree_path, order_tips=False, show_species=True, show_year=True):
    if not tree_path.exists():
        return pd.DataFrame(), pd.DataFrame()
    if str(SCRIPTS_DIR) not in sys.path:
        sys.path.insert(0, str(SCRIPTS_DIR))
    from label_tree_mrca_by_metadata import NewickParser

    root = NewickParser(tree_path.read_text()).parse()
    host_by_tip, group_by_tip, year_by_tip = host_details_lookup()
    y_positions = {}
    counter = {"y": 0}

    def tip_sort_key(node):
        clean_name, _label = extract_relax_label(node.name)
        host = str(host_by_tip.get(clean_name, "") or "Unknown")
        year = str(year_by_tip.get(clean_name, "") or "Unknown")
        year_key = 999999 if year in {"", "Unknown", "Mixed"} else int(year) if year.isdigit() else 999998
        return (host.lower(), year_key, clean_name.lower())

    def subtree_sort_key(node):
        if node.is_leaf():
            return tip_sort_key(node)
        return min(subtree_sort_key(child) for child in node.children)

    if order_tips:
        def sort_children(node):
            if node.children:
                for child in node.children:
                    sort_children(child)
                node.children.sort(key=subtree_sort_key)

        sort_children(root)

    def assign_y(node):
        if node.is_leaf():
            clean_name, _label = extract_relax_label(node.name)
            y_positions[id(node)] = counter["y"]
            counter["y"] += 1
            return y_positions[id(node)]
        child_ys = [assign_y(child) for child in node.children]
        y_positions[id(node)] = sum(child_ys) / len(child_ys)
        return y_positions[id(node)]

    assign_y(root)
    x_positions = {id(root): 0.0}
    unnamed_counter = {"n": 0}

    def assign_x(node):
        parent_x = x_positions[id(node)]
        for child in node.children:
            try:
                length = float(child.length) if child.length else 0.0
            except ValueError:
                length = 0.0
            x_positions[id(child)] = parent_x + length
            assign_x(child)

    assign_x(root)
    segments = []
    tips = []

    def branch_name(node):
        clean_name, _label = extract_relax_label(node.name)
        if clean_name:
            return clean_name
        unnamed_counter["n"] += 1
        return f"Node{unnamed_counter['n']}"

    def visit(node, parent_name=""):
        node_x = x_positions[id(node)]
        if node.children:
            child_ys = [y_positions[id(child)] for child in node.children]
            segments.append(
                {
                    "branch": branch_name(node),
                    "x": node_x,
                    "x2": node_x,
                    "y": min(child_ys),
                    "y2": max(child_ys),
                    "tip": "",
                    "host": "",
                    "host_group": "",
                    "branch_length": "",
                    "segment": "vertical",
                }
            )
            for child in node.children:
                clean_name, _label = extract_relax_label(child.name)
                child_branch = branch_name(child)
                try:
                    branch_length = float(child.length) if child.length else 0.0
                except ValueError:
                    branch_length = 0.0
                segments.append(
                    {
                        "branch": child_branch,
                        "x": node_x,
                        "x2": x_positions[id(child)],
                        "y": y_positions[id(child)],
                        "y2": y_positions[id(child)],
                    "tip": clean_name,
                    "host": host_by_tip.get(clean_name, ""),
                    "host_group": group_by_tip.get(clean_name, ""),
                    "collection_year": year_by_tip.get(clean_name, ""),
                    "branch_length": branch_length,
                    "segment": "branch",
                }
                )
                visit(child, child_branch)
        else:
            clean_name, _label = extract_relax_label(node.name)
            tips.append(
                {
                    "x": node_x,
                    "y": y_positions[id(node)],
                    "tip": clean_name,
                    "host": host_by_tip.get(clean_name, ""),
                    "host_group": group_by_tip.get(clean_name, ""),
                    "collection_year": year_by_tip.get(clean_name, ""),
                }
            )

    visit(root)
    tip_frame = pd.DataFrame(tips)
    if not tip_frame.empty:
        labels = []
        for row in tip_frame.to_dict("records"):
            parts = []
            if show_species:
                parts.append(str(row.get("host", "") or "Unknown host"))
            if show_year:
                parts.append(str(row.get("collection_year", "") or "Unknown year"))
            labels.append(" | ".join(parts))
        tip_frame["tip_label"] = labels
    return pd.DataFrame(segments), tip_frame


def absrel_tree_chart(order_tips=False, show_species=True, show_year=True):
    tree_path = ROOT / "results" / "hyphy_global_selection" / "ANDV_GPC_global.hyphy_ready.treefile"
    absrel_all = load_tsv("ANDV_GPC_ABSREL.all_branches.tsv")
    absrel_sig = load_tsv("ANDV_GPC_ABSREL.branches.tsv")
    segments, tips = generic_tree_segments(tree_path, order_tips=order_tips, show_species=show_species, show_year=show_year)
    if segments.empty:
        st.info("HyPhy-ready tree not found.")
        return

    if absrel_all.empty:
        absrel_all = absrel_sig
    if not absrel_all.empty:
        absrel_all = absrel_all.copy()
        absrel_all["corrected_p_value"] = pd.to_numeric(absrel_all["corrected_p_value"], errors="coerce")
        absrel_all["uncorrected_p_value"] = pd.to_numeric(absrel_all["uncorrected_p_value"], errors="coerce")
        absrel_all["absrel_call"] = "not significant"
        absrel_all.loc[absrel_all["uncorrected_p_value"] < 0.05, "absrel_call"] = "nominal"
        absrel_all.loc[absrel_all["corrected_p_value"] < 0.05, "absrel_call"] = "corrected significant"
        segments = segments.merge(
            absrel_all[["branch", "corrected_p_value", "uncorrected_p_value", "lrt", "absrel_call"]],
            on="branch",
            how="left",
        )
    else:
        segments["corrected_p_value"] = pd.NA
        segments["uncorrected_p_value"] = pd.NA
        segments["lrt"] = pd.NA
        segments["absrel_call"] = "not significant"
    segments["absrel_call"] = segments["absrel_call"].fillna("not significant")
    segments["stroke_width"] = segments["absrel_call"].map(
        {
            "corrected significant": 3.2,
            "nominal": 2.0,
            "not significant": 1.0,
        }
    ).fillna(1.0)

    calls = segments[segments["segment"] == "branch"]["absrel_call"].value_counts().to_dict()
    cols = st.columns(3)
    cols[0].metric("Corrected significant branches", calls.get("corrected significant", 0))
    cols[1].metric("Nominal branches", calls.get("nominal", 0))
    cols[2].metric("Tested branches shown", int((segments["segment"] == "branch").sum()))

    branch_colors = alt.Scale(
        domain=["corrected significant", "nominal", "not significant"],
        range=["#d7301f", "#fdae6b", "#cfcfcf"],
    )
    tree = (
        alt.Chart(segments)
        .mark_rule()
        .encode(
            x=alt.X("x:Q", title="Branch length from root"),
            x2="x2:Q",
            y=alt.Y("y:Q", axis=None),
            y2="y2:Q",
            color=alt.Color("absrel_call:N", scale=branch_colors, legend=alt.Legend(title="aBSREL branch result")),
            strokeWidth=alt.StrokeWidth("stroke_width:Q", legend=None),
            tooltip=[
                "branch",
                "tip",
                "host",
                "host_group",
                "collection_year",
                "absrel_call",
                alt.Tooltip("corrected_p_value:Q", format=".3g"),
                alt.Tooltip("uncorrected_p_value:Q", format=".3g"),
                alt.Tooltip("branch_length:Q", format=".4g"),
            ],
        )
    )
    tip_colors = alt.Scale(
        domain=["human", "known_nonhuman", "unknown", "ambiguous_mixed_human_nonhuman"],
        range=["#d95f02", "#1b9e77", "#969696", "#7570b3"],
    )
    tip_points = (
        alt.Chart(tips)
        .mark_circle(size=24, opacity=0.9)
        .encode(
            x="x:Q",
            y="y:Q",
            color=alt.Color("host_group:N", scale=tip_colors, legend=alt.Legend(title="Tip host group")),
            tooltip=["tip", "host", "host_group", "collection_year"],
        )
    )
    layers = tree + tip_points
    if not tips.empty and "tip_label" in tips.columns and any(tips["tip_label"].fillna("").astype(str)):
        tip_labels = (
            alt.Chart(tips)
            .mark_text(align="left", baseline="middle", dx=8, fontSize=11, color="#263238")
            .encode(
                x="x:Q",
                y="y:Q",
                text="tip_label:N",
                tooltip=["tip", "host", "host_group", "collection_year"],
            )
        )
        layers = layers + tip_labels
    st.altair_chart(
        layers
        .properties(
            height=720,
            title=alt.TitleParams(
                "aBSREL Branch-Level Selection",
                subtitle="Red/orange branches indicate episodic diversifying selection support; tips are colored by host group.",
            ),
        )
        .configure_view(stroke=None),
        use_container_width=True,
    )
    if not absrel_sig.empty:
        st.subheader("aBSREL Significant Branch Table")
        st.dataframe(absrel_sig.sort_values("corrected_p_value"), use_container_width=True)
    else:
        st.info("aBSREL detected no branches with significant episodic diversifying selection after correction.")


def relax_tree_chart(tree_path, qc, relax):
    segments, tips = relax_tree_segments(tree_path, qc)
    if segments.empty:
        st.info("RELAX-labeled tree not found yet.")
        return

    title = "RELAX Tree View: Foreground / Background Branch Labels"
    subtitle = "Foreground Test branches are human-derived; Reference branches are known non-human hosts."
    if not relax.empty:
        row = relax.iloc[0]
        k = row.get("k", "")
        p = row.get("p_value", "")
        direction = row.get("direction", "")
        subtitle = f"K = {k}, p = {p}; {direction or 'no direction parsed'}"

    branch_colors = alt.Scale(
        domain=["test", "reference", "unlabeled", "unknown", "internal"],
        range=["#d95f02", "#4d4d4d", "#d9d9d9", "#d9d9d9", "#bdbdbd"],
    )
    branch_width = alt.condition(alt.datum.group == "test", alt.value(2.8), alt.value(1.2))
    tree = (
        alt.Chart(segments)
        .mark_rule()
        .encode(
            x=alt.X("x:Q", title="Branch length from root"),
            x2="x2:Q",
            y=alt.Y("y:Q", axis=None),
            y2="y2:Q",
            color=alt.Color("group:N", scale=branch_colors, legend=alt.Legend(title="RELAX branch group")),
            strokeWidth=branch_width,
            tooltip=[
                "tip",
                "host",
                "group",
                alt.Tooltip("branch_length:Q", format=".4g"),
                "segment",
            ],
        )
    )
    tip_points = (
        alt.Chart(tips)
        .mark_circle(size=28, opacity=0.9)
        .encode(
            x="x:Q",
            y="y:Q",
            color=alt.Color("group:N", scale=branch_colors, legend=None),
            tooltip=["tip", "host", "group"],
        )
    )
    tip_labels = (
        alt.Chart(tips)
        .transform_filter(alt.datum.group == "test")
        .mark_text(align="left", dx=4, fontSize=10, color="#333")
        .encode(x="x:Q", y="y:Q", text="tip:N")
    )
    st.altair_chart(
        (tree + tip_points + tip_labels)
        .properties(height=720, title=alt.TitleParams(title, subtitle=subtitle))
        .configure_view(stroke=None),
        use_container_width=True,
    )


def relax_visualization(relax):
    qc_path = ROOT / "results" / "hyphy_host_relax" / "ANDV_GPC_human_vs_reservoir_RELAX.qc.tsv"
    tree_path = ROOT / "results" / "hyphy_host_relax" / "ANDV_GPC_human_vs_reservoir_RELAX.treefile"
    qc = pd.read_csv(qc_path, sep="\t") if qc_path.exists() else pd.DataFrame()

    st.subheader("RELAX Selection Intensity")
    st.caption("RELAX tests whether selection is relaxed or intensified in foreground branches; it is not a site-level positive-selection test.")

    if not qc.empty:
        counts = (
            qc.assign(relax_group=qc["relax_group"].fillna("").replace("", "unlabeled"))
            .groupby("relax_group")
            .size()
            .reset_index(name="branches")
        )
        colors = alt.Scale(
            domain=["test", "reference", "unknown", "unlabeled"],
            range=[FEL_COLOR_RANGE[0], MEME_COLOR_RANGE[0], "#bdbdbd", "#d9d9d9"],
        )
        branch_chart = (
            alt.Chart(counts)
            .mark_bar(cornerRadiusTopLeft=3, cornerRadiusTopRight=3)
            .encode(
                x=alt.X("relax_group:N", title="RELAX branch group"),
                y=alt.Y("branches:Q", title="Terminal branches"),
                color=alt.Color("relax_group:N", scale=colors, legend=None),
                tooltip=["relax_group", "branches"],
            )
            .properties(height=220, title="Branch Labels Used For RELAX")
        )

    if relax.empty:
        st.info("RELAX result JSON not found yet. Run `bash scripts/run_relax_by_host_group.sh`, then rebuild dashboard data.")
        if not qc.empty:
            st.altair_chart(branch_chart, use_container_width=True)
            st.dataframe(qc, use_container_width=True)
        if tree_path.exists():
            st.subheader("RELAX Tree View")
            relax_tree_chart(tree_path, qc, relax)
            with st.expander("RELAX-labeled Newick"):
                st.text_area("Newick with Test/Reference labels", tree_path.read_text(), height=180)
        return

    row = relax.iloc[0]
    k_value = pd.to_numeric(row.get("k", ""), errors="coerce")
    p_value = pd.to_numeric(row.get("p_value", ""), errors="coerce")
    q_value = pd.to_numeric(row.get("q_value", ""), errors="coerce")
    direction = row.get("direction", "") or "unknown"
    significant = (not pd.isna(q_value) and q_value <= 0.05) or (pd.isna(q_value) and not pd.isna(p_value) and p_value <= 0.05)
    result_color = "#808080"
    if significant and not pd.isna(k_value):
        result_color = "#1f9e89" if k_value < 1 else "#d95f02"
    interpretation = "no significant difference from K = 1"
    if significant and direction == "relaxed":
        interpretation = "relaxed selection in human-derived lineages"
    elif significant and direction == "intensified":
        interpretation = "intensified selection in human-derived lineages"

    with st.container(border=True):
        st.markdown("**Comparison: Human-derived vs reservoir-derived**")
        metric_row(
            {
                "foreground": "Homo sapiens",
                "background": "Known non-human hosts",
                "k": "" if pd.isna(k_value) else round(float(k_value), 4),
                "p": "" if pd.isna(p_value) else f"{p_value:.3g}",
            },
            [
                ("Foreground", "foreground"),
                ("Background", "background"),
                ("RELAX K", "k"),
                ("p-value", "p"),
            ],
        )
        st.markdown(f"**Interpretation:** {interpretation}")

    if not pd.isna(k_value):
        forest = pd.DataFrame(
            [
                {
                    "comparison": "Human-derived vs reservoir-derived",
                    "k": float(k_value),
                    "direction": direction,
                    "significant": "significant" if significant else "not significant",
                    "plot_color": result_color,
                    "p_value": p_value,
                    "q_value": q_value,
                }
            ]
        )
        domain_max = max(2, float(k_value) * 1.25)
        forest_point = (
            alt.Chart(forest)
            .mark_circle(size=240)
            .encode(
                x=alt.X("k:Q", title="RELAX K selection-intensity parameter", scale=alt.Scale(domain=[0, domain_max])),
                y=alt.Y("comparison:N", title=None),
                color=alt.Color("plot_color:N", scale=None, legend=None),
                tooltip=[
                    "comparison",
                    alt.Tooltip("k:Q", format=".4g"),
                    "direction",
                    "significant",
                    alt.Tooltip("p_value:Q", format=".3g"),
                    alt.Tooltip("q_value:Q", format=".3g"),
                ],
            )
        )
        guide = (
            alt.Chart(pd.DataFrame({"x": [1]}))
            .mark_rule(strokeDash=[5, 4], color="#111", strokeWidth=2)
            .encode(x="x:Q")
        )
        relaxed_band = (
            alt.Chart(pd.DataFrame({"start": [0], "end": [1]}))
            .mark_rect(opacity=0.12, color="#1f9e89")
            .encode(x="start:Q", x2="end:Q")
        )
        intensified_band = (
            alt.Chart(pd.DataFrame({"start": [1], "end": [domain_max]}))
            .mark_rect(opacity=0.12, color="#d95f02")
            .encode(x="start:Q", x2="end:Q")
        )
        st.altair_chart(
            (relaxed_band + intensified_band + guide + forest_point).properties(
                height=150,
                title="RELAX K Forest Plot: <1 Relaxed, =1 Neutral, >1 Intensified",
            ),
            use_container_width=True,
        )

    if not qc.empty:
        st.altair_chart(branch_chart, use_container_width=True)

    if tree_path.exists():
        st.subheader("Foreground / Background Tree Context")
        st.caption("Inspect whether Test branches are spread across the tree or concentrated in a single clade.")
        relax_tree_chart(tree_path, qc, relax)
        with st.expander("RELAX-labeled Newick"):
            st.text_area("RELAX-labeled Newick tree", tree_path.read_text(), height=180)

    st.warning(
        "Interpret RELAX cautiously: if foreground branches are concentrated in one clade, "
        "the result may reflect clade structure as much as host-derived selection intensity."
    )

    st.dataframe(relax, use_container_width=True)


def page_host_metadata():
    st.header("Host / Metadata Association")
    st.warning(
        "Host effect may be confounded with viral clade if each host corresponds to one lineage. "
        "Use terms like human-derived rather than human-adapted unless supported by explicit tests."
    )

    host_counts = load_tsv("metadata_host_counts.tsv")
    country_counts = load_tsv("metadata_country_counts.tsv")
    year_counts = load_tsv("metadata_year_counts.tsv")

    cols = st.columns(3)
    with cols[0]:
        st.subheader("Host Composition")
        if not host_counts.empty:
            st.bar_chart(host_counts.head(20).set_index("host"))
    with cols[1]:
        st.subheader("Country Composition")
        if not country_counts.empty:
            st.bar_chart(country_counts.head(20).set_index("country"))
    with cols[2]:
        st.subheader("Year Composition")
        if not year_counts.empty:
            st.bar_chart(year_counts.set_index("year"))

    st.subheader("Selection-Analysis FASTA/Tree Host Composition")
    input_host_counts = load_tsv("selection_input_host_group_counts.tsv")
    input_host_details = load_tsv("selection_input_host_group_details.tsv")
    if not input_host_counts.empty:
        st.bar_chart(input_host_counts.set_index("host_group"))
        st.dataframe(input_host_details, use_container_width=True)
    else:
        st.info("Selection input host summary not found. Run `python scripts/build_dashboard_data.py`.")

    relax = load_tsv("ANDV_GPC_RELAX.summary.tsv")
    relax_visualization(relax)

    st.subheader("Other Host-Labeled Tests")
    st.write("Recommended follow-up analyses after RELAX:")
    st.markdown(
        """
        - BUSTED foreground: selected host groups
        - CFEL: reservoir vs human-derived if enough independent sequences
        """
    )


def page_hyphy_hosts():
    st.header("HyPhy-Ready Host Split")
    st.caption("Host composition for the deduplicated alignment/tree used by the global HyPhy scan.")

    counts = load_tsv("selection_input_host_group_counts.tsv")
    details = load_tsv("selection_input_host_group_details.tsv")
    country_counts = load_tsv("selection_input_country_counts.tsv")
    year_counts = load_tsv("selection_input_year_counts.tsv")
    if counts.empty:
        st.info("HyPhy-ready host summary not found. Rebuild dashboard data first.")
        return

    values = dict(zip(counts["host_group"], counts["count"]))
    human = int(values.get("human", 0))
    other = int(values.get("known_nonhuman", 0))
    unknown = int(values.get("unknown", 0))
    ambiguous = int(values.get("ambiguous_mixed_human_nonhuman", 0))
    known_total = human + other

    cols = st.columns(4)
    cols[0].metric("Human tips", human)
    cols[1].metric("Other known-host tips", other)
    cols[2].metric("Known-host total", known_total)
    cols[3].metric("Unknown / ambiguous", unknown + ambiguous)

    if known_total:
        human_pct = round(100 * human / known_total, 1)
        other_pct = round(100 * other / known_total, 1)
        st.info(
            f"Among HyPhy-ready tips with resolved host groups: "
            f"{human_pct}% human and {other_pct}% other known hosts."
        )
    else:
        st.warning("No resolved human-versus-other host split is present in the current dashboard summary.")

    display = pd.DataFrame(
        [
            {"group": "Human", "count": human},
            {"group": "Other known host", "count": other},
            {"group": "Unknown", "count": unknown},
            {"group": "Ambiguous mixed", "count": ambiguous},
        ]
    )
    chart = (
        alt.Chart(display)
        .mark_bar(cornerRadiusTopLeft=4, cornerRadiusTopRight=4)
        .encode(
            x=alt.X(
                "group:N",
                title=None,
                sort=["Human", "Other known host", "Unknown", "Ambiguous mixed"],
            ),
            y=alt.Y("count:Q", title="HyPhy-ready tips"),
            color=alt.Color(
                "group:N",
                scale=alt.Scale(
                    domain=["Human", "Other known host", "Unknown", "Ambiguous mixed"],
                    range=["#d1495b", "#247ba0", "#8d99ae", "#6a4c93"],
                ),
                legend=None,
            ),
            tooltip=["group", "count"],
        )
        .properties(height=320)
    )
    st.altair_chart(chart, use_container_width=True)

    trend_cols = st.columns(2)
    with trend_cols[0]:
        st.subheader("Collection Country")
        if not country_counts.empty:
            country_chart = (
                alt.Chart(country_counts)
                .mark_bar(cornerRadiusTopLeft=4, cornerRadiusTopRight=4, color="#247ba0")
                .encode(
                    x=alt.X("country:N", title=None, sort="-y"),
                    y=alt.Y("count:Q", title="HyPhy-ready tips"),
                    tooltip=["country", "count"],
                )
                .properties(height=300)
            )
            st.altair_chart(country_chart, use_container_width=True)
        else:
            st.info("No HyPhy-ready country summary found.")
    with trend_cols[1]:
        st.subheader("Collection Year")
        if not year_counts.empty:
            year_chart = (
                alt.Chart(year_counts)
                .mark_bar(cornerRadiusTopLeft=4, cornerRadiusTopRight=4, color="#d1495b")
                .encode(
                    x=alt.X("year:N", title=None, sort=None),
                    y=alt.Y("count:Q", title="HyPhy-ready tips"),
                    tooltip=["year", "count"],
                )
                .properties(height=300)
            )
            st.altair_chart(year_chart, use_container_width=True)
        else:
            st.info("No HyPhy-ready collection-year summary found.")

    st.subheader("Tip-Level Host Assignment")
    if not details.empty:
        show_cols = [
            col
            for col in [
                "tip",
                "host_group",
                "hostNameScientific_values",
                "geoLocCountry_values",
                "collection_year_values",
                "duplicate_member_count",
                "metadata_ids",
            ]
            if col in details.columns
        ]
        st.dataframe(details[show_cols], use_container_width=True, hide_index=True)


def page_tree():
    st.header("Tree / Clade")
    mode = st.radio("Tree overlay", ["aBSREL branch selection", "Raw HyPhy-ready Newick"], horizontal=True)
    if mode == "aBSREL branch selection":
        control_cols = st.columns(3)
        with control_cols[0]:
            order_tips = st.checkbox("Order tips by species / year", value=True)
        with control_cols[1]:
            show_species = st.checkbox("Show species labels", value=True)
        with control_cols[2]:
            show_year = st.checkbox("Show collection year", value=True)
        absrel_tree_chart(order_tips=order_tips, show_species=show_species, show_year=show_year)
        return

    tree_path = ROOT / "results" / "hyphy_global_selection" / "ANDV_GPC_global.hyphy_ready.treefile"
    file_note(tree_path)
    if tree_path.exists():
        tree_text = tree_path.read_text()
        st.text_area("HyPhy-Ready Newick Tree", tree_text, height=220)

    absrel = load_tsv("ANDV_GPC_ABSREL.branches.tsv")
    st.subheader("aBSREL-Positive Branches")
    if absrel.empty:
        st.info("No aBSREL-positive branch table found, or no positive branches at the threshold.")
    else:
        st.dataframe(absrel, use_container_width=True)

    labeled_tree = ROOT / "results" / "iqtree2" / "GPC_CDS.mafft_codon_aligned.deduplicated.host_mrca_labeled.treefile"
    st.subheader("Optional Host MRCA-Labeled Tree")
    if labeled_tree.exists():
        file_note(labeled_tree)
        st.text_area("Host-Labeled Newick", labeled_tree.read_text(), height=220)
    else:
        st.caption(
            "Host MRCA-labeled tree has not been generated yet. The tree above is the "
            "HyPhy-ready tree used for the global FEL/MEME/BUSTED/aBSREL scan."
        )


def structure_selection_table():
    structure = read_pdb_ca(STRUCTURE_PDB)
    sites = load_tsv("ANDV_GPC_selection_annotated_sites.tsv")
    if structure.empty:
        return structure

    structure["reference_codon"] = structure["residue"]
    structure["selection_class"] = "not significant"
    structure["selection_summary"] = ""
    structure["min_p_value"] = ""
    structure["min_q_value"] = ""
    structure["region"] = ""

    if sites.empty:
        return structure

    sites = sites.copy()
    sites["reference_codon"] = pd.to_numeric(sites["reference_codon"], errors="coerce")
    sites["p_value"] = pd.to_numeric(sites["p_value"], errors="coerce")
    sites["q_value"] = pd.to_numeric(sites["q_value"], errors="coerce")

    summaries = []
    for codon, group in sites.dropna(subset=["reference_codon"]).groupby("reference_codon"):
        labels = []
        if ((group["test"] == "MEME") & (group["direction"] == "episodic diversifying")).any():
            labels.append("MEME episodic")
        if ((group["test"] == "FEL") & (group["direction"] == "diversifying")).any():
            labels.append("FEL diversifying")
        if ((group["test"] == "FEL") & (group["direction"] == "purifying")).any():
            labels.append("FEL purifying")
        selection_class = labels[0] if labels else "selected"
        summaries.append(
            {
                "reference_codon": int(codon),
                "selection_class": selection_class,
                "selection_summary": "; ".join(labels),
                "min_p_value": group["p_value"].min(),
                "min_q_value": group["q_value"].min(),
                "region": "; ".join(sorted(set(str(value) for value in group["region"].fillna("") if str(value)))),
            }
        )

    summary = pd.DataFrame(summaries)
    merged = structure.drop(columns=["selection_class", "selection_summary", "min_p_value", "min_q_value", "region"]).merge(
        summary,
        on="reference_codon",
        how="left",
    )
    merged["selection_class"] = merged["selection_class"].fillna("not significant")
    merged["selection_summary"] = merged["selection_summary"].fillna("")
    merged["region"] = merged["region"].fillna("")
    merged["min_q_value"] = merged["min_q_value"].fillna("")
    return merged


def structure_projection_plot(structure):
    if structure.empty:
        st.warning(f"Missing PDB file: `{STRUCTURE_PDB}`")
        return
    color_domain = ["not significant", "FEL purifying", "FEL diversifying", "MEME episodic"]
    color_range = ["#d9d9d9", FEL_COLOR_RANGE[0], FEL_COLOR_RANGE[1], MEME_COLOR_RANGE[0]]
    trace = (
        alt.Chart(structure)
        .mark_line(color="#bdbdbd", opacity=0.6)
        .encode(
            x=alt.X("x:Q", title="PDB CA x coordinate"),
            y=alt.Y("y:Q", title="PDB CA y coordinate"),
            order="residue:Q",
        )
    )
    points = (
        alt.Chart(structure)
        .mark_circle(opacity=0.9)
        .encode(
            x="x:Q",
            y="y:Q",
            size=alt.condition(alt.datum.selection_class == "not significant", alt.value(16), alt.value(90)),
            color=alt.Color(
                "selection_class:N",
                scale=alt.Scale(domain=color_domain, range=color_range),
                legend=alt.Legend(title="Selection result"),
            ),
            tooltip=[
                "residue",
                "aa3",
                "region",
                "selection_summary",
                alt.Tooltip("min_p_value:Q", format=".3g"),
                alt.Tooltip("min_q_value:Q", title="FDR q", format=".3g"),
                alt.Tooltip("plddt:Q", title="AlphaFold pLDDT", format=".1f"),
                "chain",
            ],
        )
    )
    st.altair_chart((trace + points).properties(height=620), use_container_width=True)


def pdb_selection_viewer(structure, show_fel=True, show_meme=True):
    if structure.empty or not STRUCTURE_PDB.exists():
        return
    selected = structure[structure["selection_class"] != "not significant"]
    meme = selected[selected["selection_class"] == "MEME episodic"]["residue"].astype(str).tolist() if show_meme else []
    fel_div = selected[selected["selection_class"] == "FEL diversifying"]["residue"].astype(str).tolist() if show_fel else []
    fel_pur = selected[selected["selection_class"] == "FEL purifying"]["residue"].astype(str).tolist() if show_fel else []
    pdb_text = STRUCTURE_PDB.read_text()
    bands = domain_bands()

    def residue_ranges(values):
        residues = sorted({int(value) for value in values if str(value).isdigit()})
        if not residues:
            return []
        ranges = []
        start = previous = residues[0]
        for residue in residues[1:]:
            if residue == previous + 1:
                previous = residue
                continue
            ranges.append(f"{start}-{previous}" if start != previous else str(start))
            start = previous = residue
        ranges.append(f"{start}-{previous}" if start != previous else str(start))
        return ranges

    def add_sphere_style_lines(ranges, color, radius):
        return " ".join(
            f'viewer.addStyle({{chain: "A", resi: "{residue_range}"}}, '
            f'{{sphere: {{radius: {radius}, color: "{color}", opacity: 0.95}}}});'
            for residue_range in ranges
        )

    def set_cartoon_style_lines(ranges, color):
        return " ".join(
            f'viewer.setStyle({{chain: "A", resi: "{residue_range}"}}, '
            f'{{cartoon: {{color: "{color}", opacity: 1.0}}}});'
            for residue_range in ranges
        )

    fel_pur_ranges = residue_ranges(fel_pur)
    fel_div_ranges = residue_ranges(fel_div)
    meme_ranges = residue_ranges(meme)
    fel_pur_style = set_cartoon_style_lines(fel_pur_ranges, FEL_COLOR_RANGE[0])
    fel_div_style = set_cartoon_style_lines(fel_div_ranges, FEL_COLOR_RANGE[1])
    meme_style = add_sphere_style_lines(meme_ranges, MEME_COLOR_RANGE[0], 1.2)

    domain_styles = []
    domain_labels = []
    domain_colors = {
        "Gn": DOMAIN_COLOR_RANGE[0],
        "WAASA": DOMAIN_COLOR_RANGE[1],
        "Gc": DOMAIN_COLOR_RANGE[2],
        "TM": DOMAIN_COLOR_RANGE[3],
        "Cytoplasmic tail": DOMAIN_COLOR_RANGE[4],
    }
    if not bands.empty:
        for row in bands.to_dict("records"):
            start = int(row["start_aa"])
            end = int(row["end_aa"])
            label = row["display"]
            color = domain_colors.get(label, "#d9d9d9")
            domain_styles.append(
                f'viewer.setStyle({{resi: "{start}-{end}"}}, {{cartoon: {{color: "{color}", opacity: 0.92}}}});'
            )
            middle = (start + end) // 2
            match = structure[structure["residue"] == middle]
            if match.empty:
                match = structure[(structure["residue"] >= start) & (structure["residue"] <= end)].head(1)
            if not match.empty:
                point = match.iloc[0]
                domain_labels.append(
                    "viewer.addLabel("
                    f'"{label}", '
                    "{"
                    f"position: {{x:{point['x']:.3f}, y:{point['y']:.3f}, z:{point['z']:.3f}}}, "
                    'backgroundColor: "white", backgroundOpacity: 0.75, '
                    f'fontColor: "{color}", fontSize: 13, inFront: true'
                    "}"
                    ");"
                )

    escaped_pdb = html.escape(pdb_text).replace("`", "\\`")
    viewer_html = f"""
    <div style="display:flex; gap:18px; flex-wrap:wrap; align-items:center; font-family:Arial, sans-serif; font-size:13px; margin:0 0 8px 0;">
      <span><b>Domains:</b></span>
      <span><span style="background:{DOMAIN_COLOR_RANGE[0]};display:inline-block;width:14px;height:14px;margin-right:4px;"></span>Gn</span>
      <span><span style="background:{DOMAIN_COLOR_RANGE[1]};display:inline-block;width:14px;height:14px;margin-right:4px;"></span>WAASA</span>
      <span><span style="background:{DOMAIN_COLOR_RANGE[2]};display:inline-block;width:14px;height:14px;margin-right:4px;"></span>Gc</span>
      <span><span style="background:{DOMAIN_COLOR_RANGE[3]};display:inline-block;width:14px;height:14px;margin-right:4px;"></span>TM</span>
      <span><span style="background:{DOMAIN_COLOR_RANGE[4]};display:inline-block;width:14px;height:14px;margin-right:4px;"></span>Cytoplasmic tail</span>
      <span style="margin-left:10px;"><b>Selection overlays:</b></span>
      <span style="opacity:{1 if show_fel else 0.35};"><span style="background:{FEL_COLOR_RANGE[0]};display:inline-block;width:22px;height:8px;margin-right:4px;"></span>FEL purifying recolored residues ({len(fel_pur)})</span>
      <span style="opacity:{1 if show_fel else 0.35};"><span style="background:{FEL_COLOR_RANGE[1]};display:inline-block;width:22px;height:8px;margin-right:4px;"></span>FEL diversifying recolored residues ({len(fel_div)})</span>
      <span style="opacity:{1 if show_meme else 0.35};"><span style="background:{MEME_COLOR_RANGE[0]};border:2px solid #111;border-radius:50%;display:inline-block;width:14px;height:14px;margin-right:4px;"></span>MEME episodic ({len(meme)})</span>
    </div>
    <div id="andv_structure_viewer" style="width: 100%; height: 620px;"></div>
    <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    <script>
      const pdb = `{escaped_pdb}`;
      const viewer = $3Dmol.createViewer("andv_structure_viewer", {{backgroundColor: "white"}});
      viewer.addModel(pdb, "pdb");
      viewer.setStyle({{}}, {{cartoon: {{color: "lightgray", opacity: 0.8}}}});
      {" ".join(domain_styles)}
      {fel_pur_style}
      {fel_div_style}
      {meme_style}
      {" ".join(domain_labels)}
      viewer.zoomTo();
      viewer.render();
    </script>
    """
    components.html(viewer_html, height=650)


def page_structure():
    st.header("PDB Structure")
    st.caption("AlphaFold model PP_006W0E7_1_97bc3 annotated with FEL and MEME selection results.")
    file_note(STRUCTURE_PDB)

    structure = structure_selection_table()
    if structure.empty:
        st.warning("No structure coordinates were loaded.")
        return

    selected = structure[structure["selection_class"] != "not significant"]
    metric_row(
        {
            "residues": len(structure),
            "selected": len(selected),
            "meme": int((structure["selection_class"] == "MEME episodic").sum()),
            "plddt": round(structure["plddt"].median(), 1),
        },
        [
            ("PDB residues", "residues"),
            ("Selected residues", "selected"),
            ("MEME-highlighted", "meme"),
            ("Median pLDDT", "plddt"),
        ],
    )

    st.subheader("Structure Projection With Selection Overlay")
    structure_projection_plot(structure)

    st.subheader("Interactive 3D Structure")
    st.caption(
        "Cartoon is colored by GPC domain architecture; FDR-significant FEL residues recolor the cartoon directly. "
        "MEME hits, if present, are shown as optional sphere overlays. "
        "Requires browser access to the 3Dmol.js library. If this panel is blank, use the projection plot above."
    )
    overlay_cols = st.columns(2)
    with overlay_cols[0]:
        show_fel = st.checkbox("Highlight FDR-significant FEL sites", value=True)
    with overlay_cols[1]:
        show_meme = st.checkbox("Highlight FDR-significant MEME sites", value=True)
    pdb_selection_viewer(structure, show_fel=show_fel, show_meme=show_meme)

    st.subheader("Selected Residues On Structure")
    st.dataframe(
        selected.sort_values(["selection_class", "min_q_value"])[
            ["residue", "aa3", "region", "selection_class", "selection_summary", "min_p_value", "min_q_value", "plddt"]
        ],
        use_container_width=True,
    )


def page_downloads():
    st.header("Downloads")
    files = [
        DASHBOARD_DIR / "ANDV_GPC_selection_annotated_sites.tsv",
        DASHBOARD_DIR / "ANDV_GPC_FEL.all_sites.tsv",
        DASHBOARD_DIR / "ANDV_GPC_MEME.all_sites.tsv",
        DASHBOARD_DIR / "ANDV_GPC_BUSTED.summary.tsv",
        DASHBOARD_DIR / "ANDV_GPC_ABSREL.branches.tsv",
        DASHBOARD_DIR / "ANDV_GPC_ABSREL.all_branches.tsv",
        DASHBOARD_DIR / "ANDV_GPC_RELAX.summary.tsv",
        DASHBOARD_DIR / "selection_input_host_group_counts.tsv",
        DASHBOARD_DIR / "selection_input_host_group_details.tsv",
        DASHBOARD_DIR / "alignment_sequence_qc.tsv",
        DASHBOARD_DIR / "pairwise_identity.tsv",
        RESULTS_DIR / "GPC_CDS.mafft_codon_aligned.deduplicated.fasta",
        RESULTS_DIR / "GPC_CDS.mafft_codon_aligned.duplicates.tsv",
        STRUCTURE_PDB,
        ROOT / "data" / "reference" / "ANDV_GPC_reference_annotation.tsv",
    ]
    for path in files:
        if path.exists():
            with open(path, "rb") as handle:
                st.download_button(path.name, handle, file_name=path.name)
        else:
            st.caption(f"Missing: `{path}`")


metrics = load_json("overview_metrics.json")

page = st.sidebar.radio(
    "Page",
    ["Overview", "QC", "Alignment", "Selection Sites", "Structure", "HyPhy Hosts", "Host/Metadata", "Tree", "Downloads"],
)

st.sidebar.caption("Reference: NC_003467.2 / CHI-7913")
st.sidebar.caption("Safe reference ID: PP_006W0E7_1")

if page == "Overview":
    page_overview(metrics)
elif page == "QC":
    page_qc(metrics)
elif page == "Alignment":
    page_alignment()
elif page == "Selection Sites":
    page_selection()
elif page == "Structure":
    page_structure()
elif page == "HyPhy Hosts":
    page_hyphy_hosts()
elif page == "Host/Metadata":
    page_host_metadata()
elif page == "Tree":
    page_tree()
elif page == "Downloads":
    page_downloads()
