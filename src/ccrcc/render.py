"""M7: deterministic HTML report from runs/<date>/findings.json (Jinja2 + inline Plotly + inline Cytoscape.js).

The LLM only wrote the narrative paragraphs inside findings.json; layout, tables, figures, citations and
provenance are produced here. The output is a single self-contained file (no server, no network needed).
"""
from __future__ import annotations

import html
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from jinja2 import Environment, FileSystemLoader, select_autoescape
from plotly.offline import get_plotlyjs

from .graph import GraphStore
from .paths import PROCESSED, REPORTS

ASSETS = Path(__file__).parent / "report_assets"


def anchor(ev_id: str) -> str:
    return "ev-" + "".join(c if c.isalnum() else "-" for c in ev_id)


def collect_ids(f: dict) -> list[str]:
    ids: list[str] = []
    for p in f["narrative"].get("summary", []):
        ids += p["evidence_ids"]
    for paras in f["narrative"].get("sections", {}).values():
        for p in paras:
            ids += p["evidence_ids"]
    for t in f.get("targets", []):
        ids += t["evidence_ids"]
    for block in f.get("agent_findings", {}).values():
        for x in block.get("findings", []):
            ids += x["evidence_ids"]
    for r in f["convergence"]["ranking"][:25]:
        for v in r["evidence_ids"].values():
            ids += v
    for block in (f.get("cell_context") or {}).values():
        if isinstance(block, list):
            for row in block:
                ids += row.get("evidence_ids", [])
    return list(dict.fromkeys(ids))


def evidence_view(e: dict) -> dict:
    """Human-readable rendering of one evidence record."""
    t = e["type"]
    link, label, detail = None, e["source"], ""
    if t == "omics_de":
        detail = (f"{e['gene']} · {e['layer']}: log2FC {e['log2fc']}, padj {e['padj']}, {e['method']} "
                  f"(pairs {e['n_pairs']}; tumor {e['n_tumor']}, normal {e['n_normal']}), {e['unit']}")
    elif t == "omics_cnv":
        detail = f"{e['gene']} ({e.get('arm') or 'NA'}): loss {e['loss_freq']}, gain {e['gain_freq']} of {e['n_tumor']} tumors"
    elif t == "omics_rna_protein_correlation":
        detail = f"{e['gene']}: Spearman ρ {e['spearman_rho']} (padj {e['padj']}, n={e['n']})"
    elif t == "literature_relation":
        link = f"https://pubmed.ncbi.nlm.nih.gov/{e['pmid']}/"
        label = f"PMID {e['pmid']}"
        detail = f"“{e['sentence']}” — {e['subject_text']} → {e['object_text']} (confidence {e['confidence']})"
    elif t == "provenance":
        link = f"https://pubmed.ncbi.nlm.nih.gov/{e['pmid']}/"
        label = f"PMID {e['pmid']} (sentence)"
        detail = e.get("text", "")
    elif t == "open_targets_association":
        link, detail = e.get("url"), f"Open Targets association score {e['score']} ({e['disease']})"
    elif t == "open_targets_drug":
        link, detail = e.get("url"), f"{e['drug']} → {e['gene']} (max clinical stage {e['max_clinical_stage']})"
    elif t == "open_targets_tractability":
        detail = (f"{e['gene']}: modalities {', '.join(e['tractable_modalities']) or 'none'}; "
                  f"approved/clinical drug: {e['has_approved_or_clinical_drug']}")
    elif t == "uniprot_entry":
        link, detail = e.get("url"), f"UniProt {e['accession']}"
    elif t == "string_interaction":
        detail = f"{e['a']}–{e['b']} STRING score {e['score']}"
    elif t == "pathway_membership":
        link = f"https://reactome.org/content/detail/{e['reactome']}"
        detail = f"{e['gene']} in {e['reactome']}"
    elif t == "sc_celltype_expression":
        detail = (f"{e['gene']} in {e['cell_group']} ({e['cell_ontology']}, {e['tissue']}): {e['pct_expressing']} % of "
                  f"{e['n_cells']} cells expressing, mean log-norm {e['mean_lognorm']}, "
                  f"{e['share_of_tumor_transcripts']} of the gene's tumour transcripts ({e['n_donors']} donors)")
    elif t == "sc_tumor_vs_pt":
        detail = (f"{e['gene']}: tumour cells vs proximal tubule log2FC {e['log2fc']} (padj {e['padj']}, "
                  f"{e['n_donors']} paired donors; log2 CPM {e['mean_log2cpm_tumor']} vs {e['mean_log2cpm_pt']})")
    elif t == "spatial_expresser":
        detail = (f"{e['gene']} in situ: {e['cell_group']} explains {e['in_situ_contribution']} of spot expression "
                  f"(β {e['in_situ_beta_cpm']} CPM per unit share, Spearman {e['spearman_with_share']}; "
                  f"composition-corrected share {e['composition_corrected_share']})")
    elif t == "spatial_neighborhood":
        detail = (f"{e['gene']}-high spots are surrounded by {e['cell_group']}: observed "
                  f"{e['observed_neighbor_share']} vs expected {e['expected_neighbor_share']} "
                  f"(log2 {e['log2_enrichment']}, padj {e['padj']}, {e['n_high_spots']} spots, "
                  f"{e['n_sections']} sections)")
    elif t == "alphafold_model":
        link = e.get("url")
        detail = (f"{e['gene']} ({e['accession']}) AlphaFold {e['model']} v{e['version']}: mean pLDDT "
                  f"{e['mean_plddt']}, {e['fraction_confident_residues']} of {e['length']} residues confident")
    elif t == "mgd_site":
        link = f"https://alphafold.ebi.ac.uk/entry/{e['accession']}"
        c = e["components"]
        detail = (f"{e['gene']} pocket {e['rank']} ({e['cavity']}): {e['volume_a3']} Å³, depth {e['max_depth_a']} Å, "
                  f"mean pLDDT {e['mean_plddt']}, {len(e['accessible_lysines'])} accessible lysines, exposed loop run "
                  f"{e['exposed_loop_run']} → MGD score {e['mgd_score']} (confidence {c['confidence']}, geometry "
                  f"{c['geometry']}, lysines {c['lysines']}, degron {c['degron']}); residues "
                  f"{', '.join(map(str, e['residues'][:10]))}")
    elif t == "ligand_site_validation":
        link = e.get("url")
        detail = (f"{e['gene']} ({e['accession']}): the {e['n_cavities']} detected cavities cover {e['recall']} of "
                  f"{e['n_known_residues']} experimentally ligand-bound residues from {e['n_ligands']} ligands "
                  f"(top pockets alone: {e['recall_in_top_pockets']}); best-matching cavity ranked "
                  f"{e['best_cavity_rank']} of {e['n_cavities']} by the MGD score")
    return {"id": e["evidence_id"], "anchor": anchor(e["evidence_id"]), "type": t, "label": label,
            "link": link, "detail": detail}


def volcano(layer: str, highlight: set[str]) -> str:
    de = pd.read_parquet(PROCESSED / "de_results.parquet")
    d = de[de["layer"] == layer].dropna(subset=["hgnc_symbol"]).copy()
    d["mlog"] = -np.log10(d["padj"].clip(lower=1e-300))
    fig = go.Figure()
    bg = d[~d["hgnc_symbol"].isin(highlight)]
    fig.add_trace(go.Scattergl(x=bg["log2fc"], y=bg["mlog"], mode="markers", name="all genes",
                               marker=dict(size=3, color="rgba(140,140,150,0.35)"),
                               text=bg["hgnc_symbol"], hovertemplate="%{text}<br>log2FC %{x:.2f}<br>-log10 padj %{y:.1f}"))
    hi = d[d["hgnc_symbol"].isin(highlight)]
    fig.add_trace(go.Scatter(x=hi["log2fc"], y=hi["mlog"], mode="markers+text", name="top-25 convergence",
                             text=hi["hgnc_symbol"], textposition="top center", textfont=dict(size=10),
                             marker=dict(size=8, color="#d1495b"),
                             hovertemplate="%{text}<br>log2FC %{x:.2f}<br>-log10 padj %{y:.1f}"))
    fig.update_layout(title=f"{layer.upper()}: tumor vs normal (paired)", xaxis_title="log2 fold change",
                      yaxis_title="-log10 adjusted p", height=460, margin=dict(l=50, r=20, t=50, b=50),
                      legend=dict(orientation="h", y=-0.2), template="plotly_white")
    return fig.to_html(full_html=False, include_plotlyjs=False, config={"displaylogo": False})


def cells_figure(cell: dict) -> str:
    """Stacked bars: which cell groups hold each candidate's transcripts in the tumour compartment."""
    ct = pd.read_parquet(PROCESSED / "sc_celltype_expression.parquet")
    genes = [r["gene"] for r in cell["cell_rows"]]
    d = ct[(ct["tissue"] == "tumor") & ct["hgnc_symbol"].isin(genes)]
    order = d[d["cell_group"] == "RCC"].set_index("hgnc_symbol")["share_of_tissue_transcripts"].reindex(genes).fillna(0)
    genes = list(order.sort_values(ascending=False).index)
    fig = go.Figure()
    for grp, label in cell["cell_groups"].items():
        s = d[d["cell_group"] == grp].set_index("hgnc_symbol")["share_of_tissue_transcripts"].reindex(genes).fillna(0)
        if s.max() < 0.02:
            continue
        fig.add_trace(go.Bar(y=genes, x=s.values, name=label, orientation="h",
                             hovertemplate="%{y}: %{x:.1%} in " + label))
    fig.update_layout(barmode="stack", title="Share of each gene's tumour transcripts by cell type (scRNA-seq)",
                      xaxis=dict(tickformat=".0%", range=[0, 1]), height=760, template="plotly_white",
                      yaxis=dict(autorange="reversed"), margin=dict(l=80, r=20, t=50, b=40),
                      legend=dict(font=dict(size=9)))
    return fig.to_html(full_html=False, include_plotlyjs=False, config={"displaylogo": False})


def spatial_figure() -> str:
    """Heatmap: abundance-weighted enrichment of cell groups around each gene's high-expression spots."""
    nb = pd.read_parquet(PROCESSED / "spatial_neighbors.parquet")
    z = nb.pivot(index="gene", columns="cell_group", values="log2_enrichment")
    p = nb.pivot(index="gene", columns="cell_group", values="padj").reindex(index=z.index, columns=z.columns)
    text = z.round(2).astype(str).where(p < 0.05, "")
    fig = go.Figure(go.Heatmap(z=z.values, x=list(z.columns), y=list(z.index), zmid=0, colorscale="RdBu_r",
                               zmin=-1, zmax=1, text=text.values, texttemplate="%{text}", textfont=dict(size=8),
                               colorbar=dict(title="log2 obs/exp"),
                               hovertemplate="%{y} high spots · neighbours %{x}<br>log2 obs/exp %{z:.2f}"))
    fig.update_layout(title="Cell groups around target-high spots (labelled: padj < 0.05)", height=760,
                      template="plotly_white", margin=dict(l=80, r=20, t=50, b=60))
    return fig.to_html(full_html=False, include_plotlyjs=False, config={"displaylogo": False})


def network(top: list[str], store: GraphStore) -> dict:
    nodes, edges, ids = {}, [], {}
    for g in top:
        n = store.resolve(g)
        if n:
            ids[n["id"]] = g
            nodes[n["id"]] = {"id": n["id"], "label": g, "kind": "gene"}
    nodes["MONDO:0005005"] = {"id": "MONDO:0005005", "label": "ccRCC", "kind": "disease"}
    for nid in list(ids):
        for nb in store.neighbors(nid, limit=500):
            if nb["predicate"] == "targets" and sum(1 for e in edges if e["target"] == nid and e["kind"] == "drug") < 3:
                nodes.setdefault(nb["id"], {"id": nb["id"], "label": nb["name"].title(), "kind": "drug"})
                edges.append({"source": nb["id"], "target": nid, "kind": "drug", "ev": nb["evidence_ids"][0]})
            elif nb["predicate"] == "interacts_with" and nb["id"] in ids and nid < nb["id"]:
                edges.append({"source": nid, "target": nb["id"], "kind": "string", "ev": nb["evidence_ids"][0]})
            elif nb["predicate"] in ("upregulated_in", "downregulated_in") and nb["id"] == "MONDO:0005005":
                edges.append({"source": nid, "target": "MONDO:0005005", "kind": nb["predicate"],
                              "ev": nb["evidence_ids"][0]})
    seen, uniq = set(), []
    for e in edges:
        k = (e["source"], e["target"], e["kind"])
        if k not in seen:
            seen.add(k)
            uniq.append(e)
    return {"nodes": [{"data": n} for n in nodes.values()],
            "edges": [{"data": {**e, "id": f"e{i}"}} for i, e in enumerate(uniq)]}


def render(findings_path: Path, out: Path | None = None) -> Path:
    f = json.loads(findings_path.read_text(encoding="utf-8"))
    store = GraphStore()
    ids = collect_ids(f)
    records = {e["evidence_id"]: e for e in store.get_evidence(ids)}
    # chunk evidence: attach sentence text for display
    for eid, e in records.items():
        if e["type"] == "provenance":
            row = store.conn.execute("SELECT n.props FROM claims c JOIN nodes n ON n.id = c.subject "
                                     "WHERE c.evidence_id = ? LIMIT 1", (eid,)).fetchone()
            if row:
                e["text"] = json.loads(row[0]).get("text", "")
    ev_views = [evidence_view(records[i]) for i in ids if i in records]
    number = {v["id"]: k + 1 for k, v in enumerate(ev_views)}
    top25 = f["convergence"]["ranking"][:25]
    top_genes = [r["gene"] for r in top25]

    env = Environment(loader=FileSystemLoader(ASSETS), autoescape=select_autoescape(["html"]))

    def cite(ev_ids):
        parts = [f'<a class="cite" href="#{anchor(i)}" title="{html.escape(i)}">{number[i]}</a>'
                 for i in ev_ids if i in number]
        return "<sup>[" + ",".join(parts) + "]</sup>" if parts else ""

    env.filters["cite"] = cite
    env.filters["evattr"] = lambda ev_ids: " ".join(ev_ids)
    env.filters["fmt"] = lambda x: "—" if x is None or (isinstance(x, float) and math.isnan(x)) else x
    tpl = env.get_template("report.html.j2")
    html_out = tpl.render(
        f=f, top25=top25, evidence=ev_views,
        cell=f["cell_context"], cells_figure=cells_figure(f["cell_context"]), spatial_figure=spatial_figure(),
        volcano_rna=volcano("rna", set(top_genes)), volcano_protein=volcano("protein", set(top_genes)),
        network_json=json.dumps(network(top_genes, store)),
        plotly_js=get_plotlyjs(), cytoscape_js=(ASSETS / "cytoscape.min.js").read_text(encoding="utf-8"),
        verification=json.loads((findings_path.parent / "verification_final.json").read_text(encoding="utf-8"))
        if (findings_path.parent / "verification_final.json").exists() else None,
        graph_stats=json.loads((PROCESSED / "graph_stats.json").read_text()),
    )
    out = out or REPORTS / "ccrcc_report.html"
    out.write_text(html_out, encoding="utf-8")
    return out
