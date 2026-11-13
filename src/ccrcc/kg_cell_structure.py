"""Load the single-cell, spatial and structural layers into the knowledge graph (M11).

Every claim goes through an Evidence record, exactly like M2-M5:

- CellType nodes (Cell Ontology ids) with `expressed_in` / `enriched_near` / `neighbor_of` claims,
  evidence types sc_celltype_expression, sc_tumor_vs_pt, spatial_expresser, spatial_neighborhood.
- Pocket nodes (one per detected cavity) with `has_pocket` claims and structural_context evidence, plus
  alphafold_model evidence for the model itself and a ligand-site validation record where PDBe-KB has one.
"""
from __future__ import annotations

import json
from datetime import date

import pandas as pd

from .graph import GraphStore, evidence_id
from .paths import PROCESSED
from .singlecell import CELL_GROUPS

TODAY = date.today().isoformat()
DATASETS = {
    "scrna": ("DATASET:LI2022-SCRNA", "Li et al. 2022 ccRCC scRNA-seq (10x 5', CZ CELLxGENE Discover)"),
    "visium": ("DATASET:LI2022-VISIUM", "Li et al. 2022 ccRCC 10x Visium sections (CZ CELLxGENE Discover)"),
    "alphafold": ("DATASET:ALPHAFOLD-DB", "AlphaFold Protein Structure Database (EMBL-EBI)"),
}
SPATIAL_LABELS = {"T_cell": ("CL:0000084", "T cell (spot level)"), "DC": ("CL:0000451", "dendritic cell (spot level)")}
MIN_PCT = 10.0      # a cell group "expresses" a gene when >= 10 % of its tumor cells have >= 1 UMI
MIN_SHARE = 0.05    # and claims are only made for groups holding >= 5 % of the gene's tumor transcripts, or
                    # above MIN_PCT (the two capture "specific but rare" and "abundant source" separately)


def _celltype_node(store: GraphStore, group: str) -> str:
    cl, label = CELL_GROUPS.get(group) or SPATIAL_LABELS.get(group, (f"GROUP:{group}", group))
    nid = f"CELLTYPE:{group}"
    store.add_node(nid, "CellType", label, cell_ontology=cl, group=group,
                   url=f"https://www.ebi.ac.uk/ols4/ontologies/cl/classes?obo_id={cl}" if cl.startswith("CL:") else None)
    return nid


def load(store: GraphStore | None = None) -> dict:
    store = store or GraphStore()
    counts = {"celltype_expression": 0, "tumor_vs_pt": 0, "spatial_expresser": 0, "spatial_neighborhood": 0,
              "alphafold_model": 0, "structural_context": 0, "ligand_site_validation": 0}
    gene_ids = {r["name"]: r["id"] for r in store.conn.execute("SELECT id, name FROM nodes WHERE type = 'Gene'")}

    sc = json.loads((PROCESSED / "sc_summary.json").read_text())
    sp_sum = json.loads((PROCESSED / "spatial_summary.json").read_text())
    struct = json.loads((PROCESSED / "structure_summary.json").read_text())
    ct = pd.read_parquet(PROCESSED / "sc_celltype_expression.parquet")
    ct = ct[(ct["tissue"] == "tumor")]
    tvp = pd.read_parquet(PROCESSED / "sc_tumor_vs_pt.parquet").set_index("hgnc_symbol")
    sx = pd.read_parquet(PROCESSED / "spatial_expressers.parquet")
    nb = pd.read_parquet(PROCESSED / "spatial_neighbors.parquet")
    targets = [t["gene"] for t in sc["targets"]]

    with store.batch():
        for key, (did, name) in DATASETS.items():
            store.add_node(did, "Dataset", name, source=name)

        # ------------------------------------------------- scRNA: which cell types express the gene
        sub = ct[ct["hgnc_symbol"].isin(targets)]
        for r in sub.itertuples():
            gid = gene_ids.get(r.hgnc_symbol)
            if gid is None or (r.pct_expressing < MIN_PCT and r.share_of_tissue_transcripts < MIN_SHARE):
                continue
            cid = _celltype_node(store, r.cell_group)
            ev = store.add_evidence(evidence_id("sc_expr", r.hgnc_symbol, r.cell_group), "sc_celltype_expression",
                                    DATASETS["scrna"][0], TODAY, gene=r.hgnc_symbol, cell_group=r.cell_group,
                                    cell_ontology=r.cell_ontology, tissue="tumor",
                                    pct_expressing=round(float(r.pct_expressing), 1),
                                    mean_lognorm=round(float(r.mean_lognorm), 3),
                                    share_of_tumor_transcripts=round(float(r.share_of_tissue_transcripts), 3),
                                    pseudobulk_cpm=None if pd.isna(r.pseudobulk_cpm) else round(float(r.pseudobulk_cpm), 1),
                                    cpm_share=None if pd.isna(getattr(r, "cpm_share", None)) else round(float(r.cpm_share), 3),
                                    n_cells=int(r.n_cells), n_donors=int(r.n_donors))
            store.add_claim(gid, "expressed_in", cid, ev)
            counts["celltype_expression"] += 1

        # ------------------------------------------------- scRNA: tumor cells vs normal proximal tubule
        for g in targets:
            gid = gene_ids.get(g)
            if gid is None or g not in tvp.index:
                continue
            r = tvp.loc[g]
            r = r.iloc[0] if isinstance(r, pd.DataFrame) else r
            cid = _celltype_node(store, "RCC")
            ev = store.add_evidence(evidence_id("sc_tvp", g), "sc_tumor_vs_pt", DATASETS["scrna"][0], TODAY,
                                    gene=g, log2fc=round(float(r.log2fc), 3), padj=float(f"{r.padj:.3g}"),
                                    n_donors=int(r.n_donors), donors=r.donors,
                                    mean_log2cpm_tumor=round(float(r.mean_log2cpm_tumor), 3),
                                    mean_log2cpm_pt=round(float(r.mean_log2cpm_pt), 3),
                                    comparison="ccRCC tumor cells (tumour) vs proximal tubule (normal kidney), "
                                               "paired donor pseudobulk, BH-adjusted")
            # like the bulk layers: a direction is claimed only when the paired test is significant (fixed 2026-11-13;
            # before, the sign alone gave "upregulated_in" / "downregulated_in" to 42 of 51 non-significant genes)
            sig = pd.notna(r.padj) and r.padj < 0.05
            pred = ("upregulated_in" if r.log2fc > 0 else "downregulated_in") if sig else "not_changed_in"
            store.add_claim(gid, pred, cid, ev)
            counts["tumor_vs_pt"] += 1

        # ------------------------------------------------- scRNA: one gene-level summary (M14)
        # Same reason as the spatial summary below: a sentence comparing cell groups needs one id that holds them all.
        for g, gs in sub.groupby("hgnc_symbol"):
            gid = gene_ids.get(g)
            if gid is None:
                continue
            gs = gs.sort_values("share_of_tissue_transcripts", ascending=False)
            gs = gs[(gs.pct_expressing >= MIN_PCT) | (gs.share_of_tissue_transcripts >= MIN_SHARE)]
            if gs.empty:
                continue
            tv = tvp.loc[g] if g in tvp.index else None
            if isinstance(tv, pd.DataFrame):
                tv = tv.iloc[0]
            ev = store.add_evidence(
                evidence_id("sc_summary", g), "sc_summary", DATASETS["scrna"][0], TODAY, gene=g, tissue="tumor",
                groups=[{"cell_group": r.cell_group, "pct_expressing": round(float(r.pct_expressing), 1),
                         "mean_lognorm": round(float(r.mean_lognorm), 3),
                         "share_of_tumor_transcripts": round(float(r.share_of_tissue_transcripts), 3),
                         "cpm_share": None if pd.isna(r.cpm_share) else round(float(r.cpm_share), 3),
                         "n_cells": int(r.n_cells)} for r in gs.itertuples()],
                tumor_vs_pt=None if tv is None else {"log2fc": round(float(tv.log2fc), 3), "padj": float(f"{tv.padj:.3g}"),
                                                     "n_donors": int(tv.n_donors)},
                method="summary of sc_celltype_expression (tumour compartment) and sc_tumor_vs_pt for this gene")
            store.add_claim(gid, "cell_context_in", DATASETS["scrna"][0], ev)
            counts["sc_summary"] = counts.get("sc_summary", 0) + 1

        # ------------------------------------------------- spatial: in situ expressers
        for r in sx.itertuples():
            gid = gene_ids.get(r.gene)
            if gid is None or r.in_situ_contribution < MIN_SHARE:
                continue
            cid = _celltype_node(store, r.cell_group)
            ev = store.add_evidence(evidence_id("spatial_expr", r.gene, r.cell_group), "spatial_expresser",
                                    DATASETS["visium"][0], TODAY, gene=r.gene, cell_group=r.cell_group,
                                    in_situ_contribution=round(float(r.in_situ_contribution), 3),
                                    in_situ_beta_cpm=round(float(r.in_situ_beta_cpm), 2),
                                    composition_corrected_share=round(float(r.composition_corrected_share), 3),
                                    spearman_with_share=round(float(r.spearman_with_share), 3),
                                    method="NNLS of spot CPM on NNLS cell-group composition (markers exclude targets)")
            store.add_claim(gid, "expressed_in_situ_in", cid, ev)
            counts["spatial_expresser"] += 1

        # ------------------------------------------------- spatial: abundance-weighted neighbourhood
        for r in nb.itertuples():
            gid = gene_ids.get(r.gene)
            if gid is None or not (r.padj < 0.05 and r.log2_enrichment > 0):
                continue
            cid = _celltype_node(store, r.cell_group)
            ev = store.add_evidence(evidence_id("spatial_nb", r.gene, r.cell_group), "spatial_neighborhood",
                                    DATASETS["visium"][0], TODAY, gene=r.gene, cell_group=r.cell_group,
                                    observed_neighbor_share=round(float(r.observed_neighbor_share), 4),
                                    expected_neighbor_share=round(float(r.expected_neighbor_share), 4),
                                    log2_enrichment=round(float(r.log2_enrichment), 3), padj=float(f"{r.padj:.3g}"),
                                    n_high_spots=int(r.n_high_spots), n_sections=int(r.n_sections),
                                    method="ring-1 neighbours of target-high spots vs section-stratified "
                                           "permutations (1000), observed/expected, BH-adjusted")
            store.add_claim(gid, "enriched_near", cid, ev)
            counts["spatial_neighborhood"] += 1

        # ------------------------------------------------- spatial: one gene-level summary (M14)
        # Agents describing a niche quote several groups at once ("RCC ... and weakly macrophage") but tended to cite
        # only the main group's record, which the numeric check then rejected. One record per gene holds every
        # in-situ expresser and every significant enriched neighbour, so one id supports the whole sentence.
        for gene, g_nb in nb.groupby("gene"):
            gid = gene_ids.get(gene)
            if gid is None:
                continue
            sig = g_nb[(g_nb.padj < 0.05) & (g_nb.log2_enrichment > 0)].sort_values("log2_enrichment", ascending=False)
            g_sx = sx[(sx.gene == gene) & (sx.in_situ_contribution >= MIN_SHARE)].sort_values("in_situ_contribution",
                                                                                               ascending=False)
            if sig.empty and g_sx.empty:
                continue
            ev = store.add_evidence(
                evidence_id("spatial_summary", gene), "spatial_summary", DATASETS["visium"][0], TODAY, gene=gene,
                n_high_spots=int(g_nb.n_high_spots.iloc[0]), n_sections=int(g_nb.n_sections.iloc[0]),
                enriched_neighbors=[{"cell_group": r.cell_group, "log2_enrichment": round(float(r.log2_enrichment), 3),
                                     "observed_neighbor_share": round(float(r.observed_neighbor_share), 4),
                                     "expected_neighbor_share": round(float(r.expected_neighbor_share), 4),
                                     "padj": float(f"{r.padj:.3g}")} for r in sig.itertuples()],
                in_situ_expressers=[{"cell_group": r.cell_group,
                                     "in_situ_contribution": round(float(r.in_situ_contribution), 3),
                                     "composition_corrected_share": round(float(r.composition_corrected_share), 3),
                                     "spearman_with_share": round(float(r.spearman_with_share), 3)}
                                    for r in g_sx.itertuples()],
                method="summary of spatial_expresser and spatial_neighborhood records for this gene")
            store.add_claim(gid, "spatial_context_in", DATASETS["visium"][0], ev)
            counts["spatial_summary"] = counts.get("spatial_summary", 0) + 1

        # ------------------------------------------------- structure: AlphaFold monomer and computational cavities
        for res in struct["results"]:
            gid = gene_ids.get(res["gene"])
            if gid is None or not res["modeled"]:
                continue
            prot = f"UNIPROT:{res['uniprot']}"
            store.add_node(prot, "Protein", res["uniprot"], accession=res["uniprot"])
            ev_model = store.add_evidence(evidence_id("alphafold", res["uniprot"], res["version"]), "alphafold_model",
                                          DATASETS["alphafold"][0], res["created"], gene=res["gene"],
                                          accession=res["uniprot"], model=res["model_entity"], version=res["version"],
                                          mean_plddt=res["mean_plddt"], length=res["length"],
                                          fraction_confident_residues=res["fraction_confident_residues"],
                                          residues_kept=res["residues_kept_above_plddt"], url=res["url"])
            store.add_claim(gid, "has_structure_model", prot, ev_model)
            counts["alphafold_model"] += 1
            v = res["ligand_site_validation"]
            if v.get("checked"):
                ev_v = store.add_evidence(evidence_id("ligand_sites", res["uniprot"]), "ligand_site_validation",
                                          "PDBe-KB ligand sites (experimental holo structures)", TODAY,
                                          gene=res["gene"], accession=res["uniprot"], n_ligands=v["n_ligands"],
                                          n_known_residues=v["n_known_residues"], recall=v["recall"],
                                          recall_in_top_pockets=v["recall_in_top_pockets"],
                                          n_cavities=v["n_cavities"], best_cavity_rank=v["best_cavity_rank"],
                                          url=f"https://www.ebi.ac.uk/pdbe/pdbe-kb/proteins/{res['uniprot']}")
                store.add_claim(gid, "pocket_detection_validated_against", prot, ev_v)
                counts["ligand_site_validation"] += 1
            for c in res["pockets"]:
                pid = f"POCKET:{res['uniprot']}:{c['cavity']}"
                store.add_node(pid, "Pocket", f"{res['gene']} pocket {c['rank']} ({c['cavity']})",
                               uniprot=res["uniprot"], cavity=c["cavity"], rank=c["rank"],
                               volume_a3=c["volume_a3"], residues=c["residues"])
                ev = store.add_evidence(evidence_id("structural_context", res["uniprot"], c["cavity"]), "structural_context",
                                        DATASETS["alphafold"][0], TODAY, gene=res["gene"], accession=res["uniprot"],
                                        cavity=c["cavity"], rank=c["rank"], volume_a3=c["volume_a3"],
                                        area_a2=c["area_a2"], max_depth_a=c["max_depth_a"],
                                        avg_hydropathy=c["avg_hydropathy"], mean_plddt=c["mean_plddt"],
                                        n_residues=c["n_residues"], residues=c["residues"][:40],
                                        accessible_lysines=c["accessible_lysines"],
                                        n_accessible_lysines=len(c["accessible_lysines"]),
                                        exposed_loop_run=c["exposed_loop_run"], components={**{k:v for k,v in c["components"].items() if k != "degron"}, "exposed_loop":c["components"]["degron"]},
                                        structural_context_score=c["mgd_score"],
                                        method="pyKVFinder cavity detection on the pLDDT-trimmed AlphaFold model; "
                                               "Descriptive monomer heuristic = weighted confidence, geometry, nearby "
                                               "lysines and exposed loops (config.yaml); no degron or degradation prediction")
                store.add_claim(gid, "has_pocket", pid, ev)
                counts["structural_context"] += 1

    out = {"claims_added": counts, "cell_types": len(CELL_GROUPS),
           "sc_cells": sc["cells_used"], "spatial_spots": sp_sum["spots"], "modeled_proteins": struct["modeled"]}
    (PROCESSED / "kg_cell_structure_stats.json").write_text(json.dumps(out, indent=1))
    return out


if __name__ == "__main__":
    print(json.dumps(load(), indent=1))
