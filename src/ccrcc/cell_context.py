"""Deterministic single-cell / spatial / structure tables for findings.json and the report.

Everything here is read from data/processed and cites evidence ids that exist in the graph; no LLM is involved.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd

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
            pocket_rows.append({"gene": r["gene"], "control": bool(r.get("control")), "uniprot": r["uniprot"], "mean_plddt": r["mean_plddt"],
                                "rank": c["rank"], "volume_a3": c["volume_a3"], "max_depth_a": c["max_depth_a"],
                                "n_lysines": len(c["accessible_lysines"]), "loop_run": c["exposed_loop_run"],
                                "mgd_score": c["mgd_score"], "recall": v.get("recall") if v.get("checked") else None,
                                "recall_top": v.get("recall_in_top_pockets") if v.get("checked") else None,
                                "evidence_ids": keep(ids)})
    pocket_rows.sort(key=lambda r: -r["mgd_score"])
    disc = discovery_block(store, st, sc)
    return {"sc": sc, "spatial": {k: v for k, v in sp_.items() if k not in ("marker_genes", "targets")},
            "structure": {"config": st["config"], "tool": st["tool"]},
            "cell_rows": cell_rows, "spatial_rows": spatial_rows, "pocket_rows": pocket_rows, **disc,
            "cell_groups": {k: v[1] for k, v in CELL_GROUPS.items()}}


def _round(x, n=2):
    return None if x is None or pd.isna(x) else round(float(x), n)


def discovery_block(store: GraphStore, st: dict, sc: dict) -> dict:
    """Rows for the discovery tab: filterable table (top 80), spotlight cards (top first-in-class per track)."""
    d = pd.read_parquet(PROCESSED / "discovery_candidates.parquet")
    summ = _read("discovery_summary.json")
    pockets = {r["gene"]: r["pockets"][0] for r in st["results"] if r["modeled"] and r["pockets"] and not r.get("control")}
    sc_t = {t["gene"]: t for t in sc["targets"] if t["measured"]}

    def valid(ids):
        return store.check_evidence_ids(list(dict.fromkeys(ids)))["valid"]

    def row(r):
        ids = [evidence_id("discovery", r.gene), evidence_id("omics_de", "rna", r.gene),
               evidence_id("omics_de", "protein", r.gene), evidence_id("sc_expr", r.gene, r.source_group),
               evidence_id("tcga_de", r.gene), evidence_id("purity", r.gene)]
        if getattr(r, "adverse_prognosis", False) is True:
            ids.append(evidence_id("tcga_survival", r.gene))
        pk = pockets.get(r.gene)
        if pk:
            ids.append(evidence_id("mgd_site", pk_acc(st, r.gene), pk["cavity"]))
        why = (f"{r.source_group} cells hold {r.source_share:.0%} of its tumour transcripts ({r.source_pct:.0f} % of cells express it); "
               f"tumor RNA log2FC {r.rna_log2fc:.2f} and protein log2FC {r.protein_log2fc:.2f}; "
               f"{int(r.ccrcc_papers)} ccRCC papers ({int(r.total_papers)} overall); "
               + ("no clinical-stage drug" if not r.clinical_drugs else f"clinical drug(s): {r.clinical_drugs}")
               + f"; tractability {r.modalities or 'none predicted'}"
               + ("; surface/secreted" if r.surface_or_secreted else "")
               + (f"; highest vital-organ expression {r.vital_max_ntpm:.0f} nTPM ({r.vital_max_tissue})"
                  if pd.notna(r.vital_max_ntpm) else "")
               + (f"; {r.clin_variable} rho {r.clin_rho:.2f}" if pd.notna(r.clin_rho) and r.clin_rho > 0 else "")
               + (f"; TCGA-KIRC log2FC {r.tcga_log2fc:.2f}" + ("" if r.replicated else " (not replicated)")
                  if pd.notna(getattr(r, "tcga_log2fc", np.nan)) else "")
               + (f"; log2FC at purity 1 {r.rna_lfc_at_purity1:.2f}" + (" (purity-robust)" if r.purity_robust else
                  " (not purity-robust)") if pd.notna(getattr(r, "rna_lfc_at_purity1", np.nan)) else "")
               + ("; higher expression, worse TCGA survival" if getattr(r, "adverse_prognosis", False) is True else "")
               + ("; selective DepMap dependency" if getattr(r, "selective_dependency", False) is True else ""))
        return {"gene": r.gene, "track": r.track, "status": r.status, "novelty": _round(r.novelty_score),
                "source_group": r.source_group, "source_share": _round(r.source_share), "source_pct": _round(r.source_pct, 0),
                "rna": _round(r.rna_log2fc), "protein": _round(r.protein_log2fc), "papers": int(r.ccrcc_papers),
                "total_papers": int(r.total_papers), "drugs": r.clinical_drugs, "modalities": r.modalities,
                "surface": bool(r.surface_or_secreted), "vital": _round(r.vital_max_ntpm, 0), "vital_tissue": r.vital_max_tissue,
                "grade": _round(r.clin_rho), "name": r.protein_name, "uniprot": r.uniprot,
                "mgd": _round(pk["mgd_score"]) if pk else None, "why": why,
                "components": {k: _round(getattr(r, f"c_{k}")) for k in ("omics", "cell_specificity", "understudied",
                                                                         "undrugged", "druggable", "safety", "aggressiveness",
                                                                         "replication", "tumor_intrinsic", "prognosis",
                                                                         "dependency") if hasattr(r, f"c_{k}")},
                "rank_in_track": int(getattr(r, "rank_in_track", 0) or 0),
                "tcga": _round(getattr(r, "tcga_log2fc", None)), "pure": _round(getattr(r, "rna_lfc_at_purity1", None)),
                "replicated": bool(getattr(r, "replicated", False) is True),
                "purity_robust": bool(getattr(r, "purity_robust", False) is True),
                "adverse_prognosis": bool(getattr(r, "adverse_prognosis", False) is True),
                "infiltration": bool(getattr(r, "infiltration_driven", False) is True),
                "evidence_ids": valid(ids)}

    # table: top 30 per track (the pool is ~75 % immune, so a global top 80 hid the tumour and stroma tracks)
    table = [row(r) for r in d.sort_values("novelty_score", ascending=False).groupby("track", sort=False).head(30)
             .sort_values(["track", "novelty_score"], ascending=[True, False]).itertuples()]
    cards = []
    for track, n in (("tumor", 6), ("stroma", 3), ("immune", 3)):
        sub = d[(d.track == track) & d.novel & d.druggable_core].sort_values("novelty_score", ascending=False).head(n)
        cards += [row(r) for r in sub.itertuples()]
    return {"discovery": {k: v for k, v in summ.items()}, "discovery_rows": table, "discovery_cards": cards}


def pk_acc(st: dict, gene: str) -> str:
    return next(r["uniprot"] for r in st["results"] if r["gene"] == gene)
