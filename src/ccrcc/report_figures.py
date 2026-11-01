"""Plotly figures for the report (deterministic; every figure is built from data/processed, never from LLM text).

Colour rules (dataviz): categorical hues in fixed slot order, never cycled; one-hue sequential ramp for magnitude;
blue <-> red diverging pair with a gray midpoint for polarity; scatter charts carry <= 3 categories (the first three
slots validate under the all-pairs rule); text uses ink tokens, never series colours. Figures are transparent and the
page re-themes fonts and grids from the CSS tokens (see report.html.j2: applyTheme).
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
import plotly.graph_objects as go

from .paths import LITERATURE, PROCESSED

CAT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
SEQ = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
DIV = [[0, "#184f95"], [0.25, "#6da7ec"], [0.5, "#f0efec"], [0.75, "#ee8f8e"], [1, "#c03130"]]
GRAY = "#898781"
TRACK_COLOR = {"tumor": CAT[0], "stroma": CAT[1], "immune": CAT[2]}
STREAM_COLOR = {"omics": CAT[0], "cnv": CAT[1], "literature": CAT[2], "genetics": CAT[3], "tractability": CAT[4],
                "dependency": CAT[5], "replication": CAT[6]}
COMPONENT_COLOR = {"omics": CAT[0], "cell_specificity": CAT[1], "understudied": CAT[2], "undrugged": CAT[3],
                   "druggable": CAT[4], "safety": CAT[5], "aggressiveness": CAT[6],
                   # M14 independent layers: one hue (last slot) in light-to-dark steps, read as one family
                   "replication": "#f1a3a2", "tumor_intrinsic": "#e34948", "prognosis": "#a92e2d", "dependency": "#6e1c1b"}
# cell groups displayed in stacked bars: fixed colour per entity, the rest folds into "other"
SHARE_GROUPS = [("RCC", "ccRCC tumor cell", CAT[0]), ("Endothelial", "Endothelial", CAT[1]),
                ("Macrophage", "Macrophage", CAT[2]), ("T", "T cell", CAT[3]), ("Fibroblast", "Fibroblast", CAT[5]),
                ("NK", "NK", CAT[4]), ("other", "Other", GRAY)]
T_GROUPS = ["CD8_T", "CD4_T", "Treg", "T_other"]
FONT = dict(family="system-ui, -apple-system, Segoe UI, sans-serif", size=12)
CFG = {"displaylogo": False, "responsive": True, "modeBarButtonsToRemove": ["lasso2d", "select2d"]}


def _html(fig: go.Figure, height: int, title: str | None = None, **layout) -> str:
    fig.update_layout(height=height, template="plotly_white", font=FONT, paper_bgcolor="rgba(0,0,0,0)",
                      plot_bgcolor="rgba(0,0,0,0)", margin=dict(l=60, r=20, t=44 if title else 16, b=48),
                      title=dict(text=title, x=0, font=dict(size=14)) if title else None,
                      legend=dict(orientation="h", y=-0.18, font=dict(size=11)), hoverlabel=dict(font=dict(size=12)),
                      **layout)
    fig.update_xaxes(gridcolor="rgba(137,135,129,.22)", zerolinecolor="rgba(137,135,129,.45)", automargin=True)
    fig.update_yaxes(gridcolor="rgba(137,135,129,.22)", zerolinecolor="rgba(137,135,129,.45)", automargin=True)
    return fig.to_html(full_html=False, include_plotlyjs=False, config=CFG)


# ------------------------------------------------------------------------------------------ omics
def volcano(layer: str, highlight: set[str]) -> str:
    d = pd.read_parquet(PROCESSED / "de_results.parquet")
    d = d[d["layer"] == layer].dropna(subset=["hgnc_symbol"]).copy()
    d["mlog"] = -np.log10(d["padj"].clip(lower=1e-300))
    bg = d[~d["hgnc_symbol"].isin(highlight)]
    hi = d[d["hgnc_symbol"].isin(highlight)]
    fig = go.Figure()
    fig.add_trace(go.Scattergl(x=bg["log2fc"], y=bg["mlog"], mode="markers", name="all genes",
                               marker=dict(size=3, color="rgba(137,135,129,.35)"), text=bg["hgnc_symbol"],
                               hovertemplate="%{text}<br>log2FC %{x:.2f} · -log10 padj %{y:.1f}<extra></extra>"))
    fig.add_trace(go.Scatter(x=hi["log2fc"], y=hi["mlog"], mode="markers+text", name="top-25 convergence",
                             text=hi["hgnc_symbol"], textposition="top center", textfont=dict(size=10),
                             marker=dict(size=8, color=CAT[0], line=dict(width=1.5, color="#fff")),
                             hovertemplate="%{text}<br>log2FC %{x:.2f} · -log10 padj %{y:.1f}<extra></extra>"))
    fig.update_xaxes(title="log2 fold change (tumor − normal)")
    fig.update_yaxes(title="−log10 adjusted p")
    return _html(fig, 430, f"{layer.upper()}: tumor vs normal (paired)")


def rna_vs_protein(highlight: set[str]) -> str:
    de = pd.read_parquet(PROCESSED / "de_results.parquet").dropna(subset=["hgnc_symbol"])
    de = de.sort_values("padj").drop_duplicates(["layer", "hgnc_symbol"])
    w = de.pivot(index="hgnc_symbol", columns="layer", values="log2fc").dropna(subset=["rna", "protein"])
    hi = w[w.index.isin(highlight)]
    bg = w[~w.index.isin(highlight)]
    fig = go.Figure()
    fig.add_trace(go.Scattergl(x=bg["rna"], y=bg["protein"], mode="markers", name="genes measured in both",
                               marker=dict(size=3, color="rgba(137,135,129,.35)"), text=bg.index,
                               hovertemplate="%{text}<br>RNA %{x:.2f} · protein %{y:.2f}<extra></extra>"))
    fig.add_trace(go.Scatter(x=hi["rna"], y=hi["protein"], mode="markers+text", text=hi.index, textposition="top center",
                             textfont=dict(size=10), name="top-25 convergence",
                             marker=dict(size=8, color=CAT[0], line=dict(width=1.5, color="#fff")),
                             hovertemplate="%{text}<br>RNA %{x:.2f} · protein %{y:.2f}<extra></extra>"))
    lim = [-6, 8]
    fig.add_trace(go.Scatter(x=lim, y=lim, mode="lines", line=dict(color=GRAY, dash="dot", width=1), name="RNA = protein",
                             hoverinfo="skip"))
    fig.update_xaxes(title="RNA log2FC")
    fig.update_yaxes(title="protein log2FC")
    return _html(fig, 430, "Do RNA and protein agree? (tumor vs normal)")


def arm_loss() -> str:
    c = pd.read_parquet(PROCESSED / "cnv_frequencies.parquet").dropna(subset=["arm"])
    g = c.groupby("arm").agg(loss=("loss_freq", "median"), gain=("gain_freq", "median"), n=("loss_freq", "size"))
    g = g[g.n >= 30].sort_values("loss", ascending=False).head(14)
    colors = [CAT[0] if a.startswith("3p") else GRAY for a in g.index]
    fig = go.Figure(go.Bar(x=g.index, y=g["loss"], marker_color=colors, marker_cornerradius=4,
                           customdata=g[["gain", "n"]], hovertemplate="%{x}: median loss freq %{y:.2f}<br>gain %{customdata[0]:.2f} · "
                                                                       "%{customdata[1]} genes<extra></extra>"))
    fig.update_yaxes(title="median loss frequency", range=[0, 1])
    return _html(fig, 330, "Chromosome-arm loss in CPTAC tumors (3p highlighted)", showlegend=False)


# ------------------------------------------------------------------------------------------ convergence
def stream_stack(ranking: list[dict], weights: dict) -> str:
    rows = ranking[:25][::-1]
    fig = go.Figure()
    for k, w in weights.items():
        fig.add_trace(go.Bar(y=[r["gene"] for r in rows], x=[w * r["components"][k] for r in rows], name=k, orientation="h",
                             marker_color=STREAM_COLOR[k], marker_line=dict(width=1.5, color="rgba(255,255,255,.9)"),
                             hovertemplate="%{y} · " + k + ": %{x:.2f}<extra></extra>"))
    fig.update_layout(barmode="stack")
    fig.update_xaxes(title="weighted contribution to the convergence score")
    return _html(fig, 640, "What drives each candidate's score? (top 25)")


def component_heatmap(ranking: list[dict], keys: list[str] | None = None) -> str:
    rows = ranking[:25]
    keys = keys or list(rows[0]["components"])   # convergence streams only (rows also carry profile components)
    z = [[r["components"][k] for k in keys] for r in rows]
    fig = go.Figure(go.Heatmap(z=z, x=keys, y=[r["gene"] for r in rows], colorscale=[[i / 6, c] for i, c in enumerate(SEQ)],
                               zmin=0, zmax=1, colorbar=dict(title="component", thickness=12),
                               hovertemplate="%{y} · %{x}: %{z:.2f}<extra></extra>"))
    fig.update_yaxes(autorange="reversed")
    return _html(fig, 640, "Evidence streams per candidate (0 – 1)")


# ------------------------------------------------------------------------------------------ cells & space
def _order_groups() -> list[str]:
    return ["RCC", "PT", "Epi_other", "Endothelial", "Fibroblast", "Macrophage", "Monocyte", "cDC", "pDC", "Mast",
            "CD8_T", "CD4_T", "Treg", "T_other", "NK", "B", "Plasma"]


def celltype_dotplot(genes: list[str], labels: dict) -> str:
    ct = pd.read_parquet(PROCESSED / "sc_celltype_expression.parquet")
    d = ct[(ct["tissue"] == "tumor") & ct["hgnc_symbol"].isin(genes)]
    order = _order_groups()
    fig = go.Figure()
    mx = float(d["mean_lognorm"].max())
    fig.add_trace(go.Scatter(
        x=[labels.get(g, g) for g in d["cell_group"]], y=d["hgnc_symbol"], mode="markers",
        marker=dict(size=np.sqrt(d["pct_expressing"].clip(lower=0.5)) * 3.2, color=d["mean_lognorm"], cmin=0, cmax=mx,
                    colorscale=[[i / 6, c] for i, c in enumerate(SEQ)], colorbar=dict(title="mean log-norm", thickness=12),
                    line=dict(width=0.6, color="rgba(137,135,129,.6)")),
        customdata=np.c_[d["pct_expressing"], d["share_of_tissue_transcripts"], d["n_cells"]],
        hovertemplate="%{y} in %{x}<br>%{customdata[0]:.1f} % of %{customdata[2]:,.0f} cells · "
                      "%{customdata[1]:.1%} of tumour transcripts<extra></extra>"))
    fig.update_xaxes(categoryorder="array", categoryarray=[labels.get(g, g) for g in order], tickangle=-40)
    fig.update_yaxes(categoryorder="array", categoryarray=genes[::-1])
    return _html(fig, max(420, 20 * len(genes) + 140), "Who expresses each candidate? Dot size = % cells, colour = mean expression",
                 showlegend=False)


def share_stack(genes: list[str]) -> str:
    ct = pd.read_parquet(PROCESSED / "sc_celltype_expression.parquet")
    d = ct[(ct["tissue"] == "tumor") & ct["hgnc_symbol"].isin(genes)]
    piv = d.pivot_table(index="hgnc_symbol", columns="cell_group", values="share_of_tissue_transcripts", aggfunc="sum").fillna(0)
    piv = piv.reindex(genes).dropna(how="all")
    series = {"RCC": piv.get("RCC"), "Endothelial": piv.get("Endothelial"), "Macrophage": piv.get("Macrophage"),
              "T": piv[[c for c in T_GROUPS if c in piv]].sum(axis=1), "Fibroblast": piv.get("Fibroblast"), "NK": piv.get("NK")}
    series["other"] = (1 - sum(v for v in series.values())).clip(lower=0)
    order = series["RCC"].sort_values(ascending=False).index.tolist()
    fig = go.Figure()
    for key, label, color in SHARE_GROUPS:
        s = series[key].reindex(order)
        fig.add_trace(go.Bar(y=order, x=s.values, name=label, orientation="h", marker_color=color,
                             marker_line=dict(width=1.5, color="rgba(255,255,255,.9)"),
                             hovertemplate="%{y}: %{x:.1%} in " + label + "<extra></extra>"))
    fig.update_layout(barmode="stack")
    fig.update_xaxes(tickformat=".0%", range=[0, 1])
    fig.update_yaxes(autorange="reversed")
    return _html(fig, max(420, 19 * len(order) + 120), "Share of each gene's tumour transcripts by cell type (scRNA-seq)")


def neighbor_heatmap(genes: list[str]) -> str:
    nb = pd.read_parquet(PROCESSED / "spatial_neighbors.parquet")
    nb = nb[nb["gene"].isin(genes)]
    z = nb.pivot(index="gene", columns="cell_group", values="log2_enrichment").reindex(genes).dropna(how="all")
    p = nb.pivot(index="gene", columns="cell_group", values="padj").reindex(index=z.index, columns=z.columns)
    text = z.round(2).astype(str).where(p < 0.05, "")
    fig = go.Figure(go.Heatmap(z=z.values, x=list(z.columns), y=list(z.index), zmid=0, colorscale=DIV, zmin=-1, zmax=1,
                               text=text.values, texttemplate="%{text}", textfont=dict(size=9),
                               colorbar=dict(title="log2 obs/exp", thickness=12),
                               hovertemplate="%{y}-high spots · neighbours %{x}<br>log2 obs/exp %{z:.2f}<extra></extra>"))
    fig.update_yaxes(autorange="reversed")
    return _html(fig, max(420, 20 * len(z) + 130), "Cell groups around target-high spots (labelled: padj < 0.05)")


SPATIAL_VARS = [("RCC", "ccRCC tumor cell share"), ("Endothelial", "Endothelial share"), ("Fibroblast", "Fibroblast share"),
                ("Macrophage", "Macrophage share"), ("T_cell", "T cell share")]


def spatial_map(region: str) -> str | None:
    path = PROCESSED / "spatial_section_maps.parquet"
    if not path.exists():
        return None
    d = pd.read_parquet(path)
    d = d[d["region"] == region]
    if d.empty:
        return None
    genes = [c[5:] for c in d.columns if c.startswith("gene_")]
    variables = [(c, lab, False) for c, lab in SPATIAL_VARS if c in d] + [(f"gene_{g}", f"{g} expression (CPM)", True) for g in genes]
    fig = go.Figure()
    buttons = []
    for i, (col, label, is_gene) in enumerate(variables):
        v = d[col].to_numpy(float)
        hi = np.quantile(v, 0.99) if is_gene else max(float(np.quantile(v, 0.99)), 0.05)
        fig.add_trace(go.Scattergl(x=d["x"], y=d["y"], mode="markers", visible=i == 0, name=label,
                                   marker=dict(size=5.5, color=v, cmin=0, cmax=max(hi, 1e-6), colorscale=[[j / 6, c] for j, c in enumerate(SEQ)],
                                               colorbar=dict(title=label.split(" (")[0], thickness=12)),
                                   hovertemplate=label + ": %{marker.color:.2f}<extra></extra>"))
        buttons.append(dict(label=label.replace(" share", "").replace(" expression (CPM)", " (gene)"), method="update",
                            args=[{"visible": [j == i for j in range(len(variables))]}]))
    fig.update_layout(updatemenus=[dict(buttons=buttons, direction="down", x=0, xanchor="left", y=1.14, yanchor="top",
                                        bgcolor="rgba(0,0,0,0)")])
    fig.update_yaxes(autorange="reversed", scaleanchor="x", showgrid=False, showticklabels=False, zeroline=False)
    fig.update_xaxes(showgrid=False, showticklabels=False, zeroline=False)
    title = "Tumour core section" if region == "core" else "Tumour–normal interface section"
    return _html(fig, 470, title, showlegend=False)


# ------------------------------------------------------------------------------------------ structure & MGD
def _pocket_frame(st: dict) -> pd.DataFrame:
    d = pd.read_parquet(PROCESSED / "structure_pockets.parquet")
    return d[d["rank"] == 1].copy()


def mgd_bubble(report_genes: list[str], discovery_genes: list[str]) -> str:
    d = _pocket_frame(None)
    if "control" not in d:
        d["control"] = False

    def kind(r):
        if r.control:
            return "known glue substrate (control)"
        return "discovery shortlist" if r.gene in discovery_genes else "report candidate"
    d["kind"] = [kind(r) for r in d.itertuples()]
    colors = {"report candidate": CAT[0], "discovery shortlist": CAT[1], "known glue substrate (control)": CAT[2]}
    fig = go.Figure()
    for k, c in colors.items():
        s = d[d["kind"] == k]
        fig.add_trace(go.Scatter(x=s["volume_a3"], y=s["mgd_score"], mode="markers+text", text=s["gene"], name=k,
                                 textposition="top center", textfont=dict(size=9),
                                 marker=dict(size=6 + 2.2 * s["n_accessible_lysines"], color=c, opacity=0.85,
                                             line=dict(width=1.5, color="#fff")),
                                 customdata=np.c_[s["max_depth_a"], s["mean_plddt"], s["n_accessible_lysines"], s["exposed_loop_run"]],
                                 hovertemplate="%{text}<br>%{x:.0f} Å³ · depth %{customdata[0]:.1f} Å · pLDDT %{customdata[1]:.0f}<br>"
                                               "lysines %{customdata[2]} · loop run %{customdata[3]} · MGD %{y:.2f}<extra></extra>"))
    fig.update_xaxes(title="pocket volume (Å³, log)", type="log")
    fig.update_yaxes(title="MGD site score", range=[0.2, 1])
    return _html(fig, 470, "Best pocket per protein: size, MGD score (bubble = accessible lysines)")


def mgd_components(weights: dict, n: int = 22) -> str:
    d = _pocket_frame(None)
    for k, w in weights.items():
        d[f"w_{k}"] = w * d[f"comp_{k}"]
    d = d.sort_values("mgd_score", ascending=False).head(n)[::-1]
    fig = go.Figure()
    for i, k in enumerate(weights):
        fig.add_trace(go.Bar(y=d["gene"], x=d[f"w_{k}"], name=k, orientation="h", marker_color=CAT[i],
                             marker_line=dict(width=1.5, color="rgba(255,255,255,.9)"),
                             hovertemplate="%{y} · " + k + ": %{x:.3f}<extra></extra>"))
    fig.update_layout(barmode="stack")
    fig.update_xaxes(title="weighted component (sums to the MGD score)")
    return _html(fig, max(420, 22 * len(d) + 110), "Why a site scores high: confidence, geometry, lysines, degron")


def ligand_validation(st: dict) -> str:
    rows = [(r["gene"], r["ligand_site_validation"]) for r in st["results"]
            if r["modeled"] and r["ligand_site_validation"].get("checked") and not r.get("control")]
    rows.sort(key=lambda x: x[1]["recall"])
    genes = [g for g, _ in rows]
    det = [v["recall"] for _, v in rows]
    top = [v["recall_in_top_pockets"] for _, v in rows]
    fig = go.Figure()
    for g, a, b in zip(genes, det, top):
        fig.add_trace(go.Scatter(x=[b, a], y=[g, g], mode="lines", line=dict(color="rgba(137,135,129,.55)", width=2),
                                 showlegend=False, hoverinfo="skip"))
    fig.add_trace(go.Scatter(x=det, y=genes, mode="markers", name="all detected cavities (detector)",
                             marker=dict(size=9, color=CAT[0], line=dict(width=1.5, color="#fff")),
                             hovertemplate="%{y}: detector covers %{x:.0%} of ligand residues<extra></extra>"))
    fig.add_trace(go.Scatter(x=top, y=genes, mode="markers", name="top 3 pockets by MGD score",
                             marker=dict(size=9, color=CAT[1], line=dict(width=1.5, color="#fff")),
                             hovertemplate="%{y}: top-3 covers %{x:.0%}<extra></extra>"))
    fig.update_xaxes(title="fraction of experimentally ligand-bound residues covered", tickformat=".0%", range=[-0.02, 1.02])
    return _html(fig, max(420, 20 * len(genes) + 120), "Validation against PDBe-KB ligand sites")


def control_ranking(controls: list[str]) -> str | None:
    d = _pocket_frame(None)
    if "control" not in d or not d["control"].any():
        return None
    d = d.sort_values("mgd_score", ascending=False).reset_index(drop=True)
    d["pct"] = 1 - d.index / max(len(d) - 1, 1)
    fig = go.Figure()
    fig.add_trace(go.Bar(x=d["gene"], y=d["mgd_score"], marker_color=[CAT[2] if c else "rgba(137,135,129,.45)" for c in d["control"]],
                         marker_cornerradius=3, hovertemplate="%{x}: %{y:.2f}<extra></extra>"))
    fig.update_xaxes(tickangle=-60, tickfont=dict(size=9))
    fig.update_yaxes(title="MGD site score", range=[0, 1])
    return _html(fig, 360, "Where do known glue substrates rank? (green = control)", showlegend=False)


# ------------------------------------------------------------------------------------------ discovery
def discovery_scatter() -> str:
    d = pd.read_parquet(PROCESSED / "discovery_candidates.parquet")
    fig = go.Figure()
    for track, color in TRACK_COLOR.items():
        s = d[d["track"] == track]
        fic = s["status"] == "first-in-class opportunity"
        for flag, sym, name in ((fic, "star", "first-in-class"), (~fic, "circle", "other")):
            t = s[flag]
            if t.empty:
                continue
            fig.add_trace(go.Scatter(
                x=t["ccrcc_papers"] + 1, y=t["protein_log2fc"], mode="markers", name=f"{track} · {name}",
                marker=dict(symbol=sym, size=5 + 14 * t["source_share"].clip(upper=1), color=color,
                            opacity=0.9 if sym == "star" else 0.45, line=dict(width=0.8, color="#fff")),
                text=t["gene"], customdata=np.c_[t["novelty_score"], t["total_papers"], t["source_group"], t["status"]],
                hovertemplate="<b>%{text}</b> (%{customdata[2]})<br>%{customdata[3]}<br>novelty %{customdata[0]:.2f} · "
                              "ccRCC papers %{x:.0f} · all papers %{customdata[1]}<br>protein log2FC %{y:.2f}<extra></extra>"))
    top = d[d["status"] == "first-in-class opportunity"].groupby("track").head(5)
    for r in top.itertuples():
        fig.add_annotation(x=np.log10(r.ccrcc_papers + 1), y=r.protein_log2fc, text=r.gene, showarrow=False, yshift=11,
                           font=dict(size=9), xref="x", yref="y")
    fig.update_xaxes(title="ccRCC papers + 1 (Europe PMC, log)", type="log")
    fig.update_yaxes(title="tumor protein log2FC (CPTAC)")
    fig.add_vline(x=16, line=dict(color=GRAY, dash="dot", width=1))
    return _html(fig, 520, "Discovery map: strong tumor-up genes with little ccRCC literature (marker size = cell-type specificity)")


def discovery_bars(weights: dict, n: int = 16) -> str:
    d = pd.read_parquet(PROCESSED / "discovery_candidates.parquet")
    fic = d[d["status"] == "first-in-class opportunity"]
    # all tumour and stroma first-in-class genes, then immune ones to fill: the immune track would fill the chart alone
    d = pd.concat([fic[fic.track != "immune"], fic[fic.track == "immune"]]).head(n)[::-1]
    weights = {k: v for k, v in weights.items() if f"c_{k}" in d.columns}
    tot = sum(weights.values())
    fig = go.Figure()
    for k, w in weights.items():
        fig.add_trace(go.Bar(y=d["gene"], x=w * d[f"c_{k}"] / tot, name=k.replace("_", " "), orientation="h",
                             marker_color=COMPONENT_COLOR[k], marker_line=dict(width=1.5, color="rgba(255,255,255,.9)"),
                             hovertemplate="%{y} · " + k + ": %{x:.3f}<extra></extra>"))
    fig.update_layout(barmode="stack")
    fig.update_xaxes(title="weighted component (weight × component ÷ total weight)")
    return _html(fig, max(420, 24 * len(d) + 120), "What makes a first-in-class candidate score high (tumour and stroma first)")


def safety_scatter() -> str:
    d = pd.read_parquet(PROCESSED / "discovery_candidates.parquet")
    d = d[d["novel"] & d["vital_max_ntpm"].notna()]
    fig = go.Figure()
    for track, color in TRACK_COLOR.items():
        s = d[d["track"] == track]
        fig.add_trace(go.Scatter(x=s["vital_max_ntpm"] + 0.1, y=s["novelty_score"], mode="markers", name=track,
                                 marker=dict(size=7, color=color, opacity=0.75, line=dict(width=0.8, color="#fff")), text=s["gene"],
                                 customdata=s["vital_max_tissue"],
                                 hovertemplate="<b>%{text}</b><br>%{customdata}: %{x:.0f} nTPM · novelty %{y:.2f}<extra></extra>"))
    fig.add_vline(x=50, line=dict(color=GRAY, dash="dot", width=1))
    fig.update_xaxes(title="highest expression in a vital organ (HPA nTPM, log)", type="log")
    fig.update_yaxes(title="novelty score")
    return _html(fig, 400, "Off-tumor risk of novel candidates (right of the line = expressed in vital organs)")


def grade_correlation(genes: list[str]) -> str:
    c = pd.read_parquet(PROCESSED / "clinical_associations.parquet")
    c = c[(c["variable"] == "grade") & c["hgnc_symbol"].isin(genes)]
    w = c.pivot_table(index="hgnc_symbol", columns="layer", values="rho").reindex(genes).dropna(how="all")
    fig = go.Figure()
    for layer, color in (("rna", CAT[0]), ("protein", CAT[1])):
        if layer in w:
            fig.add_trace(go.Bar(y=w.index, x=w[layer], name=layer, orientation="h", marker_color=color,
                                 hovertemplate="%{y} · " + layer + " ρ %{x:.2f}<extra></extra>"))
    fig.update_layout(barmode="group")
    fig.update_xaxes(title="Spearman ρ with histologic grade (tumors)", range=[-0.7, 0.7])
    fig.update_yaxes(autorange="reversed")
    return _html(fig, max(380, 26 * len(w) + 110), "Does expression rise with grade? (CPTAC, ~110 tumors)")


# ------------------------------------------------------------------------------------------ literature
def literature_charts() -> tuple[str, str]:
    rel = pd.read_parquet(LITERATURE / "relations.parquet")
    c = rel["relation"].value_counts().sort_values()
    f1 = go.Figure(go.Bar(y=c.index, x=c.values, orientation="h", marker_color=CAT[0], marker_cornerradius=4,
                          hovertemplate="%{y}: %{x} relations<extra></extra>"))
    f1.update_xaxes(title="validated relations")
    pubs = LITERATURE / "publications.parquet"
    if pubs.exists():
        p = pd.read_parquet(pubs)
        y = pd.to_numeric(p["year"], errors="coerce").dropna().astype(int).value_counts().sort_index()
        f2 = go.Figure(go.Bar(x=y.index, y=y.values, marker_color=CAT[2], marker_cornerradius=3,
                              hovertemplate="%{x}: %{y} publications<extra></extra>"))
        f2.update_yaxes(title="publications in the local corpus")
    else:
        f2 = go.Figure()
    return (_html(f1, 330, "Relation types extracted and verbatim-validated", showlegend=False),
            _html(f2, 330, "Publication year of the local corpus", showlegend=False))


def drug_stage(ranking: list[dict]) -> str:
    stages = ["PHASE_1", "PHASE_2", "PHASE_3", "PHASE_4", "APPROVAL"]
    names = {"PHASE_1": "Phase 1", "PHASE_2": "Phase 2", "PHASE_3": "Phase 3", "PHASE_4": "Phase 4", "APPROVAL": "Approved"}
    rows = [(r["gene"], r["drugs"]) for r in ranking[:25] if r["drugs"]]
    fig = go.Figure()
    if rows:
        genes = [g for g, _ in rows]
        z = [[sum(1 for d in drugs if d["max_stage"] == s) for s in stages] for _, drugs in rows]
        fig = go.Figure(go.Heatmap(z=z, x=[names[s] for s in stages], y=genes, colorscale=[[i / 6, c] for i, c in enumerate(SEQ)],
                                   zmin=0, colorbar=dict(title="drugs", thickness=12),
                                   hovertemplate="%{y} · %{x}: %{z} drugs<extra></extra>"))
        fig.update_yaxes(autorange="reversed")
    return _html(fig, max(320, 22 * len(rows) + 110), "Known drugs by maximum clinical stage (top-25 candidates)")


# ------------------------------------------------------------------------------------------ M13: validation layers
def rank_intervals(name: str, n: int = 25) -> str | None:
    p = PROCESSED / f"sensitivity_{name}.parquet"
    if not p.exists():
        return None
    d = pd.read_parquet(p).head(n).iloc[::-1]
    fig = go.Figure(go.Scatter(
        x=d["rank_median"], y=d["gene"], mode="markers", marker=dict(color=CAT[0], size=8),
        error_x=dict(type="data", symmetric=False, array=d["rank_p95"] - d["rank_median"],
                     arrayminus=d["rank_median"] - d["rank_p05"], color=CAT[0], thickness=1.4, width=3),
        customdata=np.c_[d["p_top10"], d["p_top25"], d["base_rank"]],
        hovertemplate="%{y}: median rank %{x}<br>P(top 10) %{customdata[0]:.2f} · P(top 25) %{customdata[1]:.2f}"
                      "<br>rank with configured weights %{customdata[2]}<extra></extra>"))
    fig.update_xaxes(title="rank under 1,000 random re-weightings (median, 5th–95th percentile)", type="log",
                     tickvals=[1, 2, 5, 10, 20, 50, 100, 200, 500, 1000], ticktext=["1", "2", "5", "10", "20", "50", "100",
                                                                                   "200", "500", "1000"])
    return _html(fig, max(420, 19 * len(d) + 120), f"How stable is the {name} ranking?", showlegend=False)


def mutation_recurrence(n_label: int = 10) -> str | None:
    p = PROCESSED / "layer_mutations.parquet"
    if not p.exists():
        return None
    d = pd.read_parquet(p)
    d = d[d["mutated"] >= 3].copy()
    d["mlog"] = -np.log10(d["padj"].clip(lower=1e-300))
    sig = d["padj"] < 0.05
    fig = go.Figure()
    fig.add_trace(go.Scattergl(x=d.loc[~sig, "obs_exp"], y=d.loc[~sig, "mlog"], mode="markers", name="not recurrent",
                               marker=dict(color=GRAY, size=4, opacity=.45), text=d.loc[~sig, "gene"],
                               hovertemplate="%{text}: obs/exp %{x:.2f}<extra></extra>"))
    s = d[sig]
    fig.add_trace(go.Scatter(x=s["obs_exp"], y=s["mlog"], mode="markers+text", name="recurrent (BH < 0.05)",
                             marker=dict(color=CAT[1], size=9), text=s["gene"].where(s.index.isin(s.head(n_label).index), ""),
                             textposition="top center", customdata=np.c_[s["mutated"], s["patients"]],
                             hovertemplate="%{text}: %{customdata[0]}/%{customdata[1]} tumours, obs/exp %{x:.1f}<extra></extra>"))
    fig.update_xaxes(title="observed / expected mutated tumours (coding-length background)", type="log",
                     tickvals=[0.1, 0.3, 1, 3, 10, 30, 100, 300], ticktext=["0.1", "0.3", "1", "3", "10", "30", "100", "300"])
    fig.update_yaxes(title="−log10 q (capped at 300)")
    return _html(fig, 430, "Somatic recurrence, CPTAC + TCGA-KIRC (511 tumours)")


def dependency_scatter(highlight: set[str]) -> str | None:
    p = PROCESSED / "layer_dependency.parquet"
    if not p.exists():
        return None
    d = pd.read_parquet(p).dropna(subset=["mean_effect_ccrcc", "mean_effect_other"])
    d = d[d["n_ccrcc_lines"] >= 10]
    sel = (d["padj"] < 0.1) & (d["selectivity"] > 0.3) & (d["frac_ccrcc_dependent"] >= 0.25) & ~d["common_essential"]
    lab = d[sel].nlargest(8, "selectivity")
    fig = go.Figure()
    fig.add_trace(go.Scattergl(x=d["mean_effect_other"], y=d["mean_effect_ccrcc"], mode="markers", name="genes",
                               marker=dict(color=GRAY, size=3, opacity=.35), text=d["gene"],
                               hovertemplate="%{text}: other %{x:.2f}, ccRCC %{y:.2f}<extra></extra>"))
    fig.add_trace(go.Scatter(x=d.loc[sel, "mean_effect_other"], y=d.loc[sel, "mean_effect_ccrcc"], mode="markers",
                             name="selective ccRCC dependency", marker=dict(color=CAT[1], size=7), text=d.loc[sel, "gene"],
                             hovertemplate="%{text}<extra></extra>"))
    h = d[d["gene"].isin(highlight)]
    fig.add_trace(go.Scatter(x=h["mean_effect_other"], y=h["mean_effect_ccrcc"], mode="markers", name="report candidates (hover)",
                             marker=dict(color=CAT[0], size=8, symbol="diamond"), text=h["gene"],
                             hovertemplate="%{text}: other %{x:.2f}, ccRCC %{y:.2f}<extra></extra>"))
    fig.add_trace(go.Scatter(x=lab["mean_effect_other"], y=lab["mean_effect_ccrcc"], mode="text", text=lab["gene"],
                             textposition="bottom left", showlegend=False, hoverinfo="skip"))
    lo = float(min(d["mean_effect_other"].min(), d["mean_effect_ccrcc"].min()))
    fig.add_trace(go.Scatter(x=[lo, 1], y=[lo, 1], mode="lines", line=dict(color=GRAY, dash="dot"), showlegend=False,
                             hoverinfo="skip"))
    fig.update_xaxes(title="mean Chronos effect, other lines")
    fig.update_yaxes(title="mean Chronos effect, ccRCC lines (lower = more essential)")
    return _html(fig, 470, "DepMap 26Q1: what do ccRCC lines need more than other lines?")


def tcga_replication(highlight: set[str]) -> str | None:
    p = PROCESSED / "layer_tcga.parquet"
    if not p.exists():
        return None
    t = pd.read_parquet(p).set_index("gene")
    de = pd.read_parquet(PROCESSED / "de_results.parquet")
    c = de[(de.layer == "rna")].dropna(subset=["hgnc_symbol"]).drop_duplicates("hgnc_symbol").set_index("hgnc_symbol")
    j = c[["log2fc", "padj"]].join(t[["log2fc", "padj"]], rsuffix="_tcga", how="inner").dropna()
    rho = float(j["log2fc"].corr(j["log2fc_tcga"], method="spearman"))
    sample = j.sample(min(len(j), 6000), random_state=1)
    fig = go.Figure(go.Scattergl(x=sample["log2fc"], y=sample["log2fc_tcga"], mode="markers", name="genes",
                                 marker=dict(color=GRAY, size=3, opacity=.35), text=sample.index,
                                 hovertemplate="%{text}: CPTAC %{x:.2f}, TCGA %{y:.2f}<extra></extra>"))
    h = j[j.index.isin(highlight)]
    fig.add_trace(go.Scatter(x=h["log2fc"], y=h["log2fc_tcga"], mode="markers+text", name="report candidates",
                             marker=dict(color=CAT[0], size=8), text=h.index, textposition="top center",
                             hovertemplate="%{text}<extra></extra>"))
    fig.update_xaxes(title="CPTAC RNA log2FC (tumour vs normal, 75 pairs)")
    fig.update_yaxes(title="TCGA-KIRC RNA log2FC (72 pairs)")
    return _html(fig, 450, f"External replication in TCGA-KIRC (Spearman ρ = {rho:.2f} over {len(j):,} genes)")


def purity_shift(genes: list[str]) -> str | None:
    p = PROCESSED / "layer_purity.parquet"
    if not p.exists():
        return None
    pu = pd.read_parquet(p).dropna(subset=["hgnc_symbol"]).drop_duplicates("hgnc_symbol").set_index("hgnc_symbol")
    de = pd.read_parquet(PROCESSED / "de_results.parquet")
    rna = de[de.layer == "rna"].dropna(subset=["hgnc_symbol"]).drop_duplicates("hgnc_symbol").set_index("hgnc_symbol")["log2fc"]
    rows = [(g, rna.get(g), pu.loc[g, "rna_lfc_at_purity1"], pu.loc[g, "rna_lfc_at_purity1_se"]) for g in genes
            if g in pu.index and pd.notna(rna.get(g)) and pd.notna(pu.loc[g, "rna_lfc_at_purity1"])]
    if not rows:
        return None
    d = pd.DataFrame(rows, columns=["gene", "bulk", "pure", "se"]).sort_values("pure")
    fig = go.Figure()
    for r in d.itertuples():
        fig.add_trace(go.Scatter(x=[r.bulk, r.pure], y=[r.gene, r.gene], mode="lines", line=dict(color=GRAY, width=1.5),
                                 showlegend=False, hoverinfo="skip"))
    fig.add_trace(go.Scatter(x=d["bulk"], y=d["gene"], mode="markers", name="bulk tumour vs normal",
                             marker=dict(color=GRAY, size=8), hovertemplate="%{y}: bulk %{x:.2f}<extra></extra>"))
    fig.add_trace(go.Scatter(x=d["pure"], y=d["gene"], mode="markers", name="extrapolated to purity 1",
                             marker=dict(color=CAT[0], size=9), error_x=dict(type="data", array=1.96 * d["se"], color=CAT[0]),
                             hovertemplate="%{y}: pure tumour %{x:.2f}<extra></extra>"))
    fig.update_xaxes(title="RNA log2FC (tumour − normal)")
    return _html(fig, max(380, 22 * len(d) + 120), "Infiltration or tumour cells? Bulk log2FC vs log2FC at purity 1 (ESTIMATE)")


def kinase_bars(n: int = 10) -> str | None:
    p = PROCESSED / "layer_kinase_activity.parquet"
    if not p.exists():
        return None
    k = pd.read_parquet(p)
    d = pd.concat([k.head(n), k.tail(n)]).drop_duplicates("gene").sort_values("z")
    col = [CAT[1] if (z > 0 and q < 0.1) else CAT[0] if (z < 0 and q < 0.1) else GRAY for z, q in zip(d["z"], d["padj"])]
    fig = go.Figure(go.Bar(x=d["z"], y=d["gene"], orientation="h", marker_color=col,
                           customdata=np.c_[d["n_substrate_sites"], d["padj"]],
                           hovertemplate="%{y}: z %{x:.2f}, %{customdata[0]} sites, q %{customdata[1]:.3g}<extra></extra>"))
    fig.update_xaxes(title="KSEA z (protein-adjusted site log2FC; coloured: BH < 0.1)")
    return _html(fig, max(380, 20 * len(d) + 110), "Kinase activity in tumours vs normal (CPTAC phosphosites)", showlegend=False)


def outcome_auroc(outcomes: dict | None) -> str | None:
    if not outcomes:
        return None
    fig = go.Figure()
    sets = [("clinical_targets", "RCC clinical targets"), ("temporal_holdout", "Holdout (first trial ≥ 2016)"),
            ("clearcode34", "ClearCode34 (prognostic)")]
    # ablation rows exist only for some gold sets; plot the labels every set has
    labels = [lab for lab in outcomes["clinical_targets"] if all(lab in outcomes[s[0]] for s in sets)]
    for i, lab in enumerate(labels):
        x = [s[1] for s in sets]
        y = [outcomes[s[0]][lab]["auroc"] for s in sets]
        lo = [outcomes[s[0]][lab]["auroc"] - outcomes[s[0]][lab]["auroc_ci95"][0] for s in sets]
        hi = [outcomes[s[0]][lab]["auroc_ci95"][1] - outcomes[s[0]][lab]["auroc"] for s in sets]
        fig.add_trace(go.Bar(x=x, y=y, name=lab, marker_color=CAT[i % len(CAT)],
                             error_y=dict(type="data", symmetric=False, array=hi, arrayminus=lo, color="#555")))
    fig.add_hline(y=0.5, line=dict(color=GRAY, dash="dot"))
    fig.update_yaxes(title="AUROC over the gene universe (95 % CI)", range=[0.3, 1])
    return _html(fig, 430, "Do the scores rank genes the clinic later used? (0.5 = random)", barmode="group")
