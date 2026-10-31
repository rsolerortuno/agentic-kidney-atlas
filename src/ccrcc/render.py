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
