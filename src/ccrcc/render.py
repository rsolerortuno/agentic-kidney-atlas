"""M7: deterministic HTML report from runs/<date>/findings.json (Jinja2 + inline Plotly + inline Cytoscape.js).

The LLM only wrote the narrative paragraphs inside findings.json; layout, tables, figures, citations and
provenance are produced here. The output is a single self-contained file (no server, no network needed).
"""
from __future__ import annotations

import html
import json
import math
from pathlib import Path

import pandas as pd
from jinja2 import Environment, FileSystemLoader, select_autoescape
from plotly.offline import get_plotlyjs

from . import report_figures as F
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
    elif t == "discovery_profile":
        link = e.get("url")
        detail = (f"{e['gene']} ({e['track']}, {e['status']}): novelty {e['novelty_score']}; {e['ccrcc_papers']} ccRCC / "
                  f"{e['total_papers']} total papers; clinical drugs {e.get('clinical_drugs') or 'none'}; tractability "
                  f"{e.get('modalities') or 'none'}; surface/secreted {e['surface_or_secreted']}; vital-organ max "
                  f"{e.get('vital_max_ntpm')} nTPM ({e.get('vital_max_tissue')}); "
                  f"{e.get('clinical_association') or 'no grade association'}")
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


def badges(f: dict) -> dict:
    """Per-gene one-line facts shown on the candidate cards (all from deterministic tables)."""
    cell = f["cell_context"]
    stage = ["PHASE_1", "PHASE_2", "PHASE_3", "PHASE_4", "APPROVAL"]
    out: dict[str, dict] = {}
    for r in cell["cell_rows"]:
        out.setdefault(r["gene"], {}).update(cell_label=r["top_label"], cell_share=r["tumor_cell_share"], cell_pct=r["top_pct"])
    for r in cell["pocket_rows"]:
        out.setdefault(r["gene"], {}).update(mgd=r["mgd_score"], plddt=r["mean_plddt"])
    for r in f["convergence"]["ranking"]:
        d = out.setdefault(r["gene"], {})
        d["rank"] = r["rank"]
        known = [x for x in r["drugs"] if x["max_stage"] in stage]
        d["drug"] = max(known, key=lambda x: stage.index(x["max_stage"]))["name"].title() if known else None
    return out


def render(findings_path: Path, out: Path | None = None) -> Path:
    f = json.loads(findings_path.read_text(encoding="utf-8"))
    store = GraphStore()
    ids = collect_ids(f)
    records = {e["evidence_id"]: e for e in store.get_evidence(ids)}
    for eid, e in records.items():  # chunk evidence: attach sentence text for display
        if e["type"] == "provenance":
            row = store.conn.execute("SELECT n.props FROM claims c JOIN nodes n ON n.id = c.subject "
                                     "WHERE c.evidence_id = ? LIMIT 1", (eid,)).fetchone()
            if row:
                e["text"] = json.loads(row[0]).get("text", "")
    ev_views = [evidence_view(records[i]) for i in ids if i in records]
    number = {v["id"]: k + 1 for k, v in enumerate(ev_views)}
    ranking = f["convergence"]["ranking"]
    top25 = ranking[:25]
    top_genes = [r["gene"] for r in top25]
    cell = f["cell_context"]
    st = json.loads((PROCESSED / "structure_summary.json").read_text())
    labels = cell["cell_groups"]
    report_genes = [r["gene"] for r in cell["cell_rows"]]
    disc_genes = json.loads((PROCESSED / "discovery_shortlist.json").read_text())
    dotplot_genes = report_genes[:29] + [g for g in disc_genes if g not in report_genes][:12]
    sw = f["convergence"]["config"]["stream_weights"]
    lit1, lit2 = F.literature_charts()
    figs = {
        "stream_stack": F.stream_stack(ranking, sw), "component_heatmap": F.component_heatmap(ranking),
        "volcano_rna": F.volcano("rna", set(top_genes)), "volcano_protein": F.volcano("protein", set(top_genes)),
        "rna_vs_protein": F.rna_vs_protein(set(top_genes)), "arm_loss": F.arm_loss(),
        "dotplot": F.celltype_dotplot(dotplot_genes, labels), "share_stack": F.share_stack(report_genes),
        "neighbors": F.neighbor_heatmap(report_genes + [g for g in disc_genes if g not in report_genes]),
        "map_core": F.spatial_map("core"), "map_interface": F.spatial_map("interface"),
        "mgd_bubble": F.mgd_bubble(report_genes, disc_genes),
        "mgd_components": F.mgd_components(cell["structure"]["config"]["mgd_weights"]),
        "ligand_validation": F.ligand_validation(st), "control_ranking": F.control_ranking([]),
        "discovery_scatter": F.discovery_scatter(),
        "discovery_bars": F.discovery_bars(cell["discovery"]["weights"]), "safety": F.safety_scatter(),
        "grade": F.grade_correlation(top_genes[:12] + [r["gene"] for r in cell["discovery_cards"] if r["track"] == "tumor"]),
        "lit_relations": lit1, "lit_years": lit2, "drug_stage": F.drug_stage(ranking),
    }
    pockets = pd.read_parquet(PROCESSED / "structure_pockets.parquet")
    best = pockets[pockets["rank"] == 1].sort_values("mgd_score", ascending=False).reset_index(drop=True)
    controls_view = []
    for c in st["results"]:
        if c.get("control") and c["modeled"] and c["pockets"] and (best["gene"] == c["gene"]).any():
            controls_view.append({"gene": c["gene"], "score": c["pockets"][0]["mgd_score"], "plddt": c["mean_plddt"],
                                  "rank": int(best.index[best["gene"] == c["gene"]][0]) + 1, "of": len(best)})
    kpis = {"universe": f["convergence"]["universe_size"],
            "claims": json.loads((PROCESSED / "graph_stats.json").read_text())["claims_total"],
            "evidence": len(ev_views), "cells": cell["sc"]["cells_used"], "spots": cell["spatial"]["spots"],
            "proteins": st["modeled"], "novel": cell["discovery"]["novel"], "pool": cell["discovery"]["pool"],
            "fic": cell["discovery"]["by_status"].get("first-in-class opportunity", 0)}

    env = Environment(loader=FileSystemLoader(ASSETS), autoescape=select_autoescape(["html"]))

    def cite(ev_ids):
        parts = [f'<a class="cite" href="#{anchor(i)}" title="{html.escape(i)}">{number[i]}</a>'
                 for i in ev_ids if i in number]
        return "<sup>[" + ",".join(parts) + "]</sup>" if parts else ""

    env.filters["cite"] = cite
    env.filters["evattr"] = lambda ev_ids: " ".join(ev_ids)
    env.filters["fmt"] = lambda x: "—" if x is None or (isinstance(x, float) and math.isnan(x)) else x
    tpl = env.get_template("report.html.j2")
    vpath = findings_path.parent / "verification_final.json"
    html_out = tpl.render(
        f=f, top25=top25, evidence=ev_views, cell=cell, figs=figs, kpis=kpis, badges=badges(f), controls=controls_view,
        network_json=json.dumps(network(top_genes, store)),
        plotly_js=get_plotlyjs(), cytoscape_js=(ASSETS / "cytoscape.min.js").read_text(encoding="utf-8"),
        verification=json.loads(vpath.read_text(encoding="utf-8")) if vpath.exists() else None,
        graph_stats=json.loads((PROCESSED / "graph_stats.json").read_text()),
    )
    out = out or REPORTS / "ccrcc_report.html"
    out.write_text(html_out, encoding="utf-8")
    return out
