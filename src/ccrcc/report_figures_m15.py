"""Plotly figures for the M15 analyses (treatment response, subtypes, ADC suitability, immune classes, replication).

Same rules as report_figures.py: built from data/processed only, transparent background, colours from the shared palette.
"""
from __future__ import annotations

import pandas as pd
import plotly.graph_objects as go

from .paths import PROCESSED
from .report_figures import CAT, GRAY, TRACK_COLOR, _html


def _read(name):
    p = PROCESSED / name
    return pd.read_parquet(p) if p.exists() else None


def treatment_arms(highlight: set[str]) -> str | None:
    """Per-gene log hazard ratio (per SD of expression) on sunitinib against atezolizumab + bevacizumab."""
    d = _read("layer_treatment.parquet")
    if d is None:
        return None
    a = d[d["arm"] == "atezo_bev"].set_index("gene")["cox_log_hr"]
    s = d[d["arm"] == "sunitinib"].set_index("gene")["cox_log_hr"]
    g = a.index.intersection(s.index)
    fig = go.Figure()
    fig.add_trace(go.Scattergl(x=s[g], y=a[g], mode="markers", name="all genes", text=g,
                               marker=dict(size=3, color="rgba(137,135,129,.35)"),
                               hovertemplate="%{text}<br>sunitinib log HR %{x:.2f}<br>atezolizumab + bevacizumab log HR %{y:.2f}<extra></extra>"))
    h = [x for x in g if x in highlight]
    fig.add_trace(go.Scatter(x=s[h], y=a[h], mode="markers+text", text=h, textposition="top center", name="report candidates",
                             marker=dict(size=9, color=CAT[0], symbol="diamond"),
                             hovertemplate="%{text}<br>sunitinib %{x:.2f}<br>atezolizumab + bevacizumab %{y:.2f}<extra></extra>"))
    lim = 1.5
    fig.add_trace(go.Scatter(x=[-lim, lim], y=[-lim, lim], mode="lines", line=dict(color=GRAY, dash="dot"), showlegend=False,
                             hoverinfo="skip"))
    fig.update_xaxes(title="log hazard ratio of progression on sunitinib (per SD)", range=[-lim, lim])
    fig.update_yaxes(title="log hazard ratio on atezolizumab + bevacizumab", range=[-lim, lim])
    return _html(fig, 440, "Which genes separate the two treatments? (IMmotion150)")


def treatment_controls(summary: dict) -> str | None:
    c = summary.get("signature_controls")
    if not c:
        return None
    sigs = ["T_effector", "angiogenesis", "myeloid_inflammation"]
    names = {"T_effector": "T-effector", "angiogenesis": "Angiogenesis", "myeloid_inflammation": "Myeloid inflammation"}
    arms = [("sunitinib", "Sunitinib", CAT[1]), ("atezo_bev", "Atezolizumab + bevacizumab", CAT[0]),
            ("atezo", "Atezolizumab", CAT[2])]
    fig = go.Figure()
    for key, label, color in arms:
        fig.add_trace(go.Bar(x=[names[s] for s in sigs], y=[c[key][s]["hr_per_sd"] for s in sigs], name=label, marker_color=color,
                             text=[f"p {c[key][s]['pvalue']:.2g}" for s in sigs], textposition="outside",
                             hovertemplate="%{x}<br>HR %{y:.2f}<extra>" + label + "</extra>"))
    fig.add_hline(y=1, line_color=GRAY, line_dash="dot")
    fig.update_yaxes(title="hazard ratio of progression per SD", type="log", range=[-0.45, 0.4])
    fig.update_layout(barmode="group")
    return _html(fig, 380, "Published signatures reproduce in this download")


def subtype_counts(summary: dict) -> str | None:
    dr = summary.get("drivers")
    if not dr:
        return None
    layers = [("cptac_rna", "CPTAC RNA", CAT[0]), ("cptac_protein", "CPTAC protein", CAT[1]), ("tcga_rna", "TCGA-KIRC RNA", CAT[2])]
    fig = go.Figure()
    for key, label, color in layers:
        fig.add_trace(go.Bar(x=list(dr), y=[dr[g]["significant_by_layer"].get(key, 0) for g in dr], name=label, marker_color=color))
    fig.update_yaxes(title="genes that differ, mutant vs wild type (q < 0.1)")
    fig.update_layout(barmode="group")
    return _html(fig, 360, "How much each driver changes the tumour")


def ccab_figure() -> str | None:
    d = _read("subtype_ccab_calls.parquet")
    if d is None:
        return None
    fig = go.Figure()
    for coh, color in (("CPTAC", CAT[0]), ("TCGA-KIRC", CAT[1])):
        fig.add_trace(go.Violin(y=d[d["cohort"] == coh]["ccab_score"], name=coh, box_visible=True, meanline_visible=True,
                                line_color=color, fillcolor="rgba(0,0,0,0)", points=False))
    fig.add_hline(y=0, line_color=GRAY, line_dash="dot", annotation_text="above 0: ccA, below 0: ccB",
                  annotation_position="top left")
    fig.update_yaxes(title="ccA score minus ccB score")
    return _html(fig, 360, "ccA / ccB score of every tumour")


def adc_figure() -> str | None:
    d = _read("adc_suitability.parquet")
    if d is None:
        return None
    s = d[(d["c_surface"] > 0) & d["c_homogeneity"].notna()].copy()
    s["x"] = s["protein_log2fc"].fillna(s["rna_log2fc"] / 2)
    lab = pd.concat([s.nsmallest(14, "rank"), s[s["control"]]]).drop_duplicates("gene")
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=s["x"], y=s["rcc_pct_donor_p25"], mode="markers", name="surface proteins", text=s["gene"],
                             marker=dict(size=7, color=s["ihc_vital_max_level"].fillna(1.5), cmin=0, cmax=3, opacity=.8,
                                         colorscale=[[0, "#cde2fb"], [1, "#e34948"]],
                                         colorbar=dict(title="staining in<br>vital organs", tickvals=[0, 1, 2, 3],
                                                       ticktext=["none", "low", "medium", "high"], len=.7)),
                             hovertemplate="%{text}<br>tumour log2FC %{x:.2f}<br>tumour cells expressing %{y:.0f}%<extra></extra>"))
    fig.add_trace(go.Scatter(x=lab["x"], y=lab["rcc_pct_donor_p25"], mode="text", text=lab["gene"], textposition="top center",
                             showlegend=False, hoverinfo="skip"))
    fig.update_xaxes(title="change in tumour (protein log2 fold change, or RNA / 2)")
    fig.update_yaxes(title="tumour cells expressing it, weakest quarter of donors (%)")
    return _html(fig, 470, "Surface proteins: how much, how evenly, and where else")


def immune_classes_figure() -> str | None:
    d = _read("immune_classes.parquet")
    if d is None:
        return None
    fig = go.Figure()
    for col, label, color in (("cd8", "CD8 T cells", CAT[0]), ("treg", "Regulatory T cells", CAT[3]), ("m2", "M2 macrophages", CAT[1])):
        fig.add_trace(go.Box(x=d["immune_class"], y=d[col], name=label, marker_color=color, boxpoints=False))
    fig.update_layout(boxmode="group")
    fig.update_xaxes(categoryorder="array", categoryarray=["cold", "intermediate", "hot"],
                     title="immune class (xCell ImmuneScore tertile)")
    fig.update_yaxes(title="CIBERSORT absolute score")
    return _html(fig, 360, "Immune cells in cold, intermediate and hot tumours")


def sc_replication_figure() -> str | None:
    d = _read("sc_replication.parquet")
    disc = _read("discovery_candidates.parquet")
    if d is None or disc is None:
        return None
    m = disc[["gene", "track"]].merge(d, on="gene")
    fig = go.Figure()
    for track, color in TRACK_COLOR.items():
        t = m[m["track"] == track]
        fig.add_trace(go.Scatter(x=t["li_rcc_share"], y=t["zhang_rcc_share"], mode="markers", name=track, text=t["gene"],
                                 marker=dict(size=6, color=color, opacity=.7),
                                 hovertemplate="%{text}<br>Li 2022 %{x:.2f}<br>Zhang 2021 %{y:.2f}<extra></extra>"))
    fig.add_trace(go.Scatter(x=[0, 1], y=[0, 1], mode="lines", line=dict(color=GRAY, dash="dot"), showlegend=False,
                             hoverinfo="skip"))
    fig.update_xaxes(title="share of tumour transcripts from tumour cells, Li 2022")
    fig.update_yaxes(title="same share, Zhang 2021")
    return _html(fig, 420, "Does a second dataset give the same cell attribution?")


def confirmatory_scatter(kind: str) -> str | None:
    """Discovery vs confirmatory log2 fold change for proteins or phosphosites measured in both cohorts."""
    d = _read("confirmatory_protein.parquet" if kind == "protein" else "confirmatory_phospho.parquet")
    if d is None:
        return None
    label = "gene" if kind == "protein" else "site_id"
    d = d.dropna(subset=["disc_log2fc", "conf_log2fc"])
    rep = d[d["replicates"]]
    non = d[~d["replicates"]]
    fig = go.Figure()
    fig.add_trace(go.Scattergl(x=non["disc_log2fc"], y=non["conf_log2fc"], mode="markers", name="does not replicate", text=non[label],
                               marker=dict(size=3, color="rgba(137,135,129,.4)"),
                               hovertemplate="%{text}<br>discovery %{x:.2f}<br>confirmatory %{y:.2f}<extra></extra>"))
    fig.add_trace(go.Scattergl(x=rep["disc_log2fc"], y=rep["conf_log2fc"], mode="markers", name="replicates", text=rep[label],
                               marker=dict(size=3, color=CAT[0], opacity=.55),
                               hovertemplate="%{text}<br>discovery %{x:.2f}<br>confirmatory %{y:.2f}<extra></extra>"))
    if kind == "protein":
        hl = d[d["gene"].isin(["CA9", "NDUFA4L2", "ENPP3", "NNMT", "SCARB1", "SLC16A3", "HK2", "PSAT1", "PCK1", "ALDOB", "EGFR", "MET"])]
        fig.add_trace(go.Scatter(x=hl["disc_log2fc"], y=hl["conf_log2fc"], mode="markers+text", text=hl["gene"], textposition="top center",
                                 name="candidates", marker=dict(size=8, color=CAT[1], symbol="diamond"),
                                 hovertemplate="%{text}<br>discovery %{x:.2f}<br>confirmatory %{y:.2f}<extra></extra>"))
    lo = float(min(d["disc_log2fc"].min(), d["conf_log2fc"].min()))
    hi = float(max(d["disc_log2fc"].max(), d["conf_log2fc"].max()))
    fig.add_trace(go.Scatter(x=[lo, hi], y=[lo, hi], mode="lines", line=dict(color=GRAY, dash="dot"), showlegend=False, hoverinfo="skip"))
    fig.update_xaxes(title="change in tumour, discovery cohort (log2 fold change)")
    fig.update_yaxes(title="change in tumour, confirmatory cohort (log2 fold change)")
    return _html(fig, 440, "Proteins in two cohorts" if kind == "protein" else "Phosphosites in two cohorts")
