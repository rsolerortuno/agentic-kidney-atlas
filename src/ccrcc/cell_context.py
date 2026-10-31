"""Deterministic single-cell / spatial / structure tables for findings.json and the report.

Everything here is read from data/processed and cites evidence ids that exist in the graph; no LLM is involved.
"""
from __future__ import annotations

import json

from .graph import GraphStore, evidence_id
from .paths import PROCESSED
from .singlecell import CELL_GROUPS


def _read(name: str):
    return json.loads((PROCESSED / name).read_text(encoding="utf-8"))


def build(store: GraphStore | None = None) -> dict:
    store = store or GraphStore()
    sc, sp_, st = _read("sc_summary.json"), _read("spatial_summary.json"), _read("structure_summary.json")

    def keep(ids):
        return store.check_evidence_ids(list(dict.fromkeys(ids)))["valid"]

    cell_rows, spatial_rows, pocket_rows = [], [], []
    for t in sc["targets"]:
        if not t["measured"]:
            continue
        g, top = t["gene"], t["top_groups"][0]
        ids = [evidence_id("sc_expr", g, x["cell_group"]) for x in t["top_groups"]]
        ids.append(evidence_id("sc_tvp", g))
        cell_rows.append({"gene": g, "top_group": top["cell_group"], "top_label": top["label"],
                          "top_pct": top["pct_expressing"], "top_share": top["share_of_tumor_transcripts"],
                          "tumor_cell_share": t["tumor_cell_share"], "tau": t["tau"], "tvp": t["tumor_vs_pt"],
                          "evidence_ids": keep(ids)})
    for t in sp_["targets"]:
        if not t["measured"]:
            continue
        g = t["gene"]
        ids = [evidence_id("spatial_expr", g, x) for x in (t["top_in_situ_expresser"], t["top_corrected_source"]) if x]
        ids += [evidence_id("spatial_nb", g, e["cell_group"]) for e in t.get("enriched_neighbors", [])]
        spatial_rows.append({"gene": g, "top_expresser": t["top_in_situ_expresser"], "corrected": t["top_corrected_source"],
                             "n_high": t["n_high_spots"], "neighbors": t.get("enriched_neighbors", [])[:5],
                             "evidence_ids": keep(ids)})
    for r in st["results"]:
        if not r["modeled"] or not r["pockets"]:
            continue
        v = r["ligand_site_validation"]
        for c in r["pockets"][:1]:
            ids = [evidence_id("mgd_site", r["uniprot"], c["cavity"]), evidence_id("alphafold", r["uniprot"], r["version"])]
            if v.get("checked"):
                ids.append(evidence_id("ligand_sites", r["uniprot"]))
            pocket_rows.append({"gene": r["gene"], "uniprot": r["uniprot"], "mean_plddt": r["mean_plddt"],
                                "rank": c["rank"], "volume_a3": c["volume_a3"], "max_depth_a": c["max_depth_a"],
                                "n_lysines": len(c["accessible_lysines"]), "loop_run": c["exposed_loop_run"],
                                "mgd_score": c["mgd_score"], "recall": v.get("recall") if v.get("checked") else None,
                                "recall_top": v.get("recall_in_top_pockets") if v.get("checked") else None,
                                "evidence_ids": keep(ids)})
    pocket_rows.sort(key=lambda r: -r["mgd_score"])
    return {"sc": sc, "spatial": {k: v for k, v in sp_.items() if k not in ("marker_genes", "targets")},
            "structure": {"config": st["config"], "tool": st["tool"]},
            "cell_rows": cell_rows, "spatial_rows": spatial_rows, "pocket_rows": pocket_rows,
            "cell_groups": {k: v[1] for k, v in CELL_GROUPS.items()}}
