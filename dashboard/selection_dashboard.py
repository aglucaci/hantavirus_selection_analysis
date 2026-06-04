#!/usr/bin/env python3
"""SelectionScope dashboard for the current hantavirus HyPhy run."""

from __future__ import annotations

import os
from pathlib import Path

import pandas as pd
import streamlit as st

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_RESULTS = ROOT / "results" / "ANDV_trees_aln-hyphy"
RESULTS = Path(os.environ.get("HYPHY_RESULTS_DIR", DEFAULT_RESULTS))
if not RESULTS.is_absolute():
    RESULTS = ROOT / RESULTS
TABLES = RESULTS / "dashboard_tables"


def load_table(name: str) -> pd.DataFrame:
    path = TABLES / name
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t")


def metric(label: str, value) -> None:
    st.metric(label, "0" if pd.isna(value) else value)


def severity_order(value: str) -> int:
    return {"high": 0, "medium": 1, "info": 2}.get(str(value), 3)


st.set_page_config(page_title="SelectionScope", layout="wide")
st.title("SelectionScope: Hantavirus HyPhy Results")
st.caption("Interactive command center for the current FEL, MEME, aBSREL, RELAX, and MSS run.")

analysis = load_table("analysis_summary.tsv")
qc = load_table("qc_summary.tsv")
dropped = load_table("dropped_sequences.tsv")
fel = load_table("fel_sites.tsv")
meme = load_table("meme_sites.tsv")
meme_branch_ebf = load_table("meme_branch_ebf.tsv")
absrel = load_table("absrel_branches.tsv")
relax = load_table("relax_results.tsv")
mss = load_table("mss_results.tsv")
warnings = load_table("warnings.tsv")

if analysis.empty:
    st.error("Dashboard tables are missing. Run `python scripts/build_hyphy_dashboard_tables.py` first.")
    st.stop()

segments = sorted(analysis["segment"].dropna().unique())
label_sets = sorted(analysis["label_set"].dropna().unique())
methods_present = analysis.query("status == 'pass'")["method"].nunique()
tabs = st.tabs(["Overview", "QC", "FEL", "MEME", "aBSREL", "RELAX", "MSS", "Warnings", "Export"])

with tabs[0]:
    st.subheader("Run Overview")
    cols = st.columns(8)
    with cols[0]:
        metric("Segments", len(segments))
    with cols[1]:
        metric("Tree label sets", len(label_sets))
    with cols[2]:
        metric("Completed methods", methods_present)
    with cols[3]:
        if not fel.empty and "q_value" in fel:
            fel_diversifying = int(((fel["direction"] == "diversifying") & (pd.to_numeric(fel["q_value"], errors="coerce") <= 0.1)).sum())
        else:
            fel_diversifying = int((fel["direction"] == "diversifying").sum()) if "direction" in fel else len(fel)
        metric("Foreground FEL diversifying sites q <= 0.1", fel_diversifying)
    with cols[4]:
        meme_count = int((pd.to_numeric(meme["q_value"], errors="coerce") <= 0.1).sum()) if not meme.empty and "q_value" in meme else len(meme)
        metric("Foreground MEME sites q <= 0.1", meme_count)
    with cols[5]:
        metric("Foreground aBSREL selected branches", int((absrel["status"] == "selected").sum()) if not absrel.empty else 0)
    with cols[6]:
        mss_models = int(pd.to_numeric(mss["model_count"], errors="coerce").fillna(0).sum()) if not mss.empty and "model_count" in mss else 0
        metric("MSS models", mss_models)
    with cols[7]:
        metric("High warnings", int((warnings["severity"] == "high").sum()) if not warnings.empty else 0)

    st.subheader("Method Status")
    status = analysis.pivot_table(
        index=["segment", "label_set"],
        columns="method",
        values="status",
        aggfunc="first",
        fill_value="missing",
    )
    st.dataframe(status, width="stretch")

    st.subheader("Selection Burden")
    burden = analysis.query("method in ['FEL', 'MEME', 'aBSREL', 'RELAX', 'MSS'] and status == 'pass'").copy()
    burden["run"] = burden["segment"] + " " + burden["label_set"] + " " + burden["method"]
    burden["significant_count"] = pd.to_numeric(burden["significant_count"], errors="coerce").fillna(0)
    st.bar_chart(burden.set_index("run")["significant_count"])

with tabs[1]:
    st.subheader("Input QC")
    st.dataframe(qc, width="stretch")
    if not dropped.empty:
        st.subheader("Dropped Sequences")
        st.dataframe(dropped, width="stretch")
    st.info("Trees are not reinferred here. Invalid sequences are removed from the FASTA and pruned from the prepared labeled trees before HyPhy runs.")

with tabs[2]:
    st.subheader("FEL Foreground-Branch Pervasive Site-Level Selection")
    left, right = st.columns(2)
    segment_filter = left.multiselect("Segment", segments, default=segments, key="fel_segment")
    label_filter = right.multiselect("Tree label set", label_sets, default=label_sets, key="fel_label")
    view = fel[fel["segment"].isin(segment_filter) & fel["label_set"].isin(label_filter)] if not fel.empty else fel
    if "direction" in view and "q_value" in view:
        fdr_view = view[pd.to_numeric(view["q_value"], errors="coerce") <= 0.1]
        st.metric("Displayed foreground FEL diversifying sites q <= 0.1", int((fdr_view["direction"] == "diversifying").sum()))
    elif "direction" in view:
        fdr_view = view
        st.metric("Displayed FEL diversifying sites", int((view["direction"] == "diversifying").sum()))
    else:
        fdr_view = view
        st.metric("Displayed FEL sites", len(view))
    if not fdr_view.empty:
        chart = fdr_view.groupby(["segment", "label_set"]).size().reset_index(name="sites")
        chart["run"] = chart["segment"] + " " + chart["label_set"]
        st.bar_chart(chart.set_index("run")["sites"])
    st.dataframe(view, width="stretch")

with tabs[3]:
    st.subheader("MEME Foreground-Branch Episodic Site-Level Selection")
    left, right = st.columns(2)
    segment_filter = left.multiselect("Segment", segments, default=segments, key="meme_segment")
    label_filter = right.multiselect("Tree label set", label_sets, default=label_sets, key="meme_label")
    view = meme[meme["segment"].isin(segment_filter) & meme["label_set"].isin(label_filter)] if not meme.empty else meme
    fdr_view = view[pd.to_numeric(view["q_value"], errors="coerce") <= 0.1] if not view.empty and "q_value" in view else view
    st.metric("Displayed foreground MEME sites q <= 0.1", len(fdr_view))
    if not fdr_view.empty:
        chart = fdr_view.groupby(["segment", "label_set"]).size().reset_index(name="sites")
        chart["run"] = chart["segment"] + " " + chart["label_set"]
        st.bar_chart(chart.set_index("run")["sites"])
    st.dataframe(view, width="stretch")
    if not meme_branch_ebf.empty:
        st.subheader("MEME Foreground Branch EBF Table")
        branch_view = meme_branch_ebf[
            meme_branch_ebf["segment"].isin(segment_filter) & meme_branch_ebf["label_set"].isin(label_filter)
        ].copy()
        branch_view["reconstructed_ebf"] = pd.to_numeric(branch_view["reconstructed_ebf"], errors="coerce")
        branch_view = branch_view.sort_values(["reconstructed_ebf", "posterior_positive_class"], ascending=[False, False]).head(200)
        st.caption("Reconstructed from MEME foreground-branch posterior annotations for raw p <= 0.10 sites; FDR status marks BH-corrected site support.")
        st.dataframe(branch_view, width="stretch")

with tabs[4]:
    st.subheader("aBSREL Foreground Branch-Level Selection")
    if absrel.empty:
        st.info("No foreground aBSREL branches passed the dashboard display threshold.")
    else:
        st.dataframe(absrel, width="stretch")
        chart = absrel.groupby(["segment", "label_set", "status"]).size().reset_index(name="branches")
        chart["run"] = chart["segment"] + " " + chart["label_set"] + " " + chart["status"]
        st.bar_chart(chart.set_index("run")["branches"])

with tabs[5]:
    st.subheader("RELAX Selection Intensification or Relaxation")
    st.dataframe(relax, width="stretch")
    if not relax.empty:
        chart = relax.copy()
        chart["run"] = chart["segment"] + " " + chart["label_set"]
        chart["k"] = pd.to_numeric(chart["k"], errors="coerce")
        st.bar_chart(chart.set_index("run")["k"])
    st.warning("RELAX rows with convergence/local-maxima warnings should be treated cautiously, especially when interpreting K.")

with tabs[6]:
    st.subheader("MSS-GA Synonymous-Rate Class Search")
    st.dataframe(mss, width="stretch")
    if not mss.empty:
        chart = mss.copy()
        chart["run"] = chart["segment"] + " " + chart["label_set"]
        chart["model_count"] = pd.to_numeric(chart["model_count"], errors="coerce")
        st.bar_chart(chart.set_index("run")["model_count"])

with tabs[7]:
    st.subheader("Robustness Warnings")
    if warnings.empty:
        st.success("No dashboard warnings were generated.")
    else:
        display = warnings.copy()
        display["_order"] = display["severity"].map(severity_order)
        display = display.sort_values(["_order", "segment", "label_set", "method"]).drop(columns="_order")
        st.dataframe(display, width="stretch")

with tabs[8]:
    st.subheader("Export Tables")
    for filename in [
        "analysis_summary.tsv",
        "qc_summary.tsv",
        "dropped_sequences.tsv",
        "fel_sites.tsv",
        "meme_sites.tsv",
        "meme_branch_ebf.tsv",
        "absrel_branches.tsv",
        "relax_results.tsv",
        "mss_results.tsv",
        "warnings.tsv",
    ]:
        path = TABLES / filename
        if path.exists():
            st.download_button(
                label=f"Download {filename}",
                data=path.read_text(),
                file_name=filename,
                mime="text/tab-separated-values",
            )

    if not fel.empty and "q_value" in fel:
        fel_sites = int(((fel["direction"] == "diversifying") & (pd.to_numeric(fel["q_value"], errors="coerce") <= 0.1)).sum())
    else:
        fel_sites = int((fel["direction"] == "diversifying").sum()) if "direction" in fel else len(fel)
    meme_sites = int((pd.to_numeric(meme["q_value"], errors="coerce") <= 0.1).sum()) if not meme.empty and "q_value" in meme else len(meme)
    selected_branches = int((absrel["status"] == "selected").sum()) if not absrel.empty else 0
    relax_sig = int(relax["significant"].astype(str).str.lower().eq("true").sum()) if not relax.empty else 0
    mss_models = int(pd.to_numeric(mss["model_count"], errors="coerce").fillna(0).sum()) if not mss.empty and "model_count" in mss else 0
    paragraph = (
        f"Across the current hantavirus analyses, foreground-branch FEL identified {fel_sites} pervasive diversifying candidate "
        f"sites at BH FDR q <= 0.1, MEME identified {meme_sites} episodic candidate sites at BH FDR q <= 0.1, "
        f"aBSREL identified {selected_branches} selected "
        f"branches after correction, RELAX found {relax_sig} significant branch-set shifts, and "
        f"MSS-GA evaluated {mss_models} synonymous-rate class models. "
        "QC warnings include duplicate-sequence warnings from HyPhy and the S-segment "
        "KP202360.1 stop-codon removal."
    )
    st.text_area("Draft results text", paragraph, height=130)
