"""cell-mcp: single-cell, spatial, structural and target-discovery context (read-only, never raw matrices)."""
from __future__ import annotations

import json
from functools import lru_cache

import pandas as pd
from mcp.server.fastmcp import FastMCP

from ..graph import evidence_id
from ..paths import PROCESSED
from .common import applied_record, canonical_gene, clean, existing

mcp = FastMCP("cell-mcp")


@lru_cache(maxsize=1)
def _ct() -> pd.DataFrame:
    return pd.read_parquet(PROCESSED / "sc_celltype_expression.parquet")


@lru_cache(maxsize=1)
def _tvp() -> pd.DataFrame:
    return pd.read_parquet(PROCESSED / "sc_tumor_vs_pt.parquet")


@lru_cache(maxsize=1)
def _spatial() -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    return (pd.read_parquet(PROCESSED / "spatial_expressers.parquet"),
            pd.read_parquet(PROCESSED / "spatial_neighbors.parquet"),
            json.loads((PROCESSED / "spatial_summary.json").read_text()))


@lru_cache(maxsize=1)
def _structure() -> dict:
    return {r["gene"]: r for r in json.loads((PROCESSED / "structure_summary.json").read_text())["results"]}


@mcp.tool()
def celltype_expression(gene: str, tissue: str = "tumor", n: int = 6) -> list[dict]:
    """Which cell types express a gene in ccRCC (Li et al. 2022 scRNA-seq, 12 donors).

    tissue: 'tumor', 'interface' (tumour-normal) or 'normal'. Returns, per cell group ranked by mean
    log-normalized expression: % of cells with >= 1 UMI, donor-averaged pseudobulk CPM, and the group's share of
    all transcripts of that gene in the tissue — the last one answers whether a bulk signal is tumor-intrinsic
    or contributed by infiltrating cells. cpm_share is the abundance-independent version (per-cell level of the
    group relative to the other groups): use it to say which cell type expresses the gene most per cell, since
    transcript shares also depend on how many cells of each type were captured. evidence_id_summary (tumour tissue)
    holds every group of the gene plus the tumour-vs-proximal-tubule result: cite it when a sentence mentions more than
    one cell group. evidence_id is present for groups
    loaded into the graph."""
    g = canonical_gene(gene)
    d = _ct()
    sub = d[(d["hgnc_symbol"] == g) & (d["tissue"] == tissue)]
    if sub.empty:
        return [{"gene": g, "message": f"not measured in the scRNA-seq {tissue} compartment"}]
    sub = sub.sort_values("mean_lognorm", ascending=False).head(min(n, 20))
    return [clean({"gene": g, "tissue": tissue, "cell_group": r.cell_group, "cell_ontology": r.cell_ontology,
                   "pct_expressing": round(float(r.pct_expressing), 1),
                   "mean_lognorm": round(float(r.mean_lognorm), 3),
                   "pseudobulk_cpm": None if pd.isna(r.pseudobulk_cpm) else round(float(r.pseudobulk_cpm), 1),
                   "share_of_tissue_transcripts": round(float(r.share_of_tissue_transcripts), 3),
                   "cpm_share": None if pd.isna(getattr(r, "cpm_share", None)) else round(float(r.cpm_share), 3),
                   "n_cells": int(r.n_cells), "n_donors": int(r.n_donors),
                   "evidence_id_summary": existing(evidence_id("sc_summary", g)) if tissue == "tumor" else None,
                   "evidence_id": existing(evidence_id("sc_expr", g, r.cell_group)) if tissue == "tumor" else None})
            for r in sub.itertuples()]


@mcp.tool()
def tumor_cell_vs_proximal_tubule(gene: str) -> dict:
    """Pseudobulk differential expression of a gene between ccRCC tumor cells and normal proximal tubule cells
    (the cell of origin), paired by donor. This removes the stromal and immune contribution that bulk
    tumor-vs-normal comparisons mix in."""
    g = canonical_gene(gene)
    d = _tvp()
    sub = d[d["hgnc_symbol"] == g]
    if sub.empty:
        return {"gene": g, "message": "not testable (absent from the scRNA-seq pseudobulks)"}
    r = sub.iloc[0]
    return clean({"gene": g, "log2fc": round(float(r.log2fc), 3), "padj": float(f"{r.padj:.3g}"),
                  "pvalue": float(f"{r.pvalue:.3g}"), "n_donors": int(r.n_donors),
                  "mean_log2cpm_tumor_cells": round(float(r.mean_log2cpm_tumor), 3),
                  "mean_log2cpm_proximal_tubule": round(float(r.mean_log2cpm_pt), 3),
                  "evidence_id": existing(evidence_id("sc_tvp", g))})


@mcp.tool()
def spatial_context(gene: str) -> dict:
    """In situ context of a gene (Li et al. 2022 10x Visium, tumour core + tumour-normal interface).

    'in_situ_expressers': cell groups whose spot abundance explains the gene's spot expression (NNLS), with the
    composition-corrected share of its transcripts. 'enriched_neighbors': cell groups over-represented in the
    ring-1 neighbours of the gene's high-expression spots relative to section-matched random spots
    (observed/expected, so abundant cell types are not favoured), BH-adjusted over all targets x groups.
    'evidence_id_summary' covers every value returned here: cite it when a sentence mentions more than one cell group."""
    g = canonical_gene(gene)
    sx, nb, summ = _spatial()
    e, n = sx[sx["gene"] == g], nb[nb["gene"] == g]
    if e.empty and n.empty:
        return {"gene": g, "message": "not covered by the Visium sections"}
    return clean({
        "gene": g, "spots": summ["spots"], "sections": summ["sections"], "donors": summ["donors"],
        "evidence_id_summary": existing(evidence_id("spatial_summary", g)),
        "in_situ_expressers": [{"cell_group": r.cell_group,
                                "in_situ_contribution": round(float(r.in_situ_contribution), 3),
                                "composition_corrected_share": round(float(r.composition_corrected_share), 3),
                                "spearman_with_share": round(float(r.spearman_with_share), 3),
                                "evidence_id": existing(evidence_id("spatial_expr", g, r.cell_group))}
                               for r in e.sort_values("in_situ_contribution", ascending=False).head(5).itertuples()],
        "enriched_neighbors": [{"cell_group": r.cell_group, "log2_enrichment": round(float(r.log2_enrichment), 3),
                                "observed_neighbor_share": round(float(r.observed_neighbor_share), 4),
                                "expected_neighbor_share": round(float(r.expected_neighbor_share), 4),
                                "padj": float(f"{r.padj:.3g}"), "n_high_spots": int(r.n_high_spots),
                                "evidence_id": existing(evidence_id("spatial_nb", g, r.cell_group))}
                               for r in n[(n["padj"] < 0.05) & (n["log2_enrichment"] > 0)]
                               .sort_values("log2_enrichment", ascending=False).head(6).itertuples()],
    })


@mcp.tool()
def structure_pockets(gene: str) -> dict:
    """AlphaFold model quality and detected cavities with a descriptive structural-context heuristic in [0, 1].

    The score weights model confidence (pLDDT of the lining residues), pocket geometry (volume and depth),
    nearby solvent-accessible lysines and exposed non-helical loops. These describe a predicted monomer;
    they do not identify degrons, establish ligase compatibility or predict degradation.
    'ligand_site_validation' reports coverage of experimental ligand-contact residues by detected cavities;
    it validates cavity detection only. A low score does not exclude molecular-glue feasibility."""
    g = canonical_gene(gene)
    res = _structure().get(g)
    if not res:
        return {"gene": g, "message": "no structural analysis for this gene"}
    if not res["modeled"]:
        return {"gene": g, "modeled": False, "reason": res["reason"]}
    return clean({
        "gene": g, "uniprot": res["uniprot"], "alphafold": res["url"], "model_version": res["version"],
        "mean_plddt": res["mean_plddt"], "fraction_confident_residues": res["fraction_confident_residues"],
        "residues_kept_above_plddt50": res["residues_kept_above_plddt"],
        "model_evidence_id": existing(evidence_id("alphafold", res["uniprot"], res["version"])),
        "pockets": [{"rank": c["rank"], "cavity": c["cavity"], "volume_a3": c["volume_a3"],
                     "max_depth_a": c["max_depth_a"], "avg_hydropathy": c["avg_hydropathy"],
                     "mean_plddt": c["mean_plddt"], "n_residues": c["n_residues"],
                     "key_residues": c["residues"][:12], "accessible_lysines": c["accessible_lysines"],
                     "exposed_loop_run": c["exposed_loop_run"], "components": {**{k:v for k,v in c["components"].items() if k != "degron"}, "exposed_loop":c["components"]["degron"]},
                     "structural_context_score": c["mgd_score"],
                     "evidence_id": existing(evidence_id("structural_context", res["uniprot"], c["cavity"]))}
                    for c in res["pockets"]],
        "ligand_site_validation": {**res["ligand_site_validation"],
                                   "evidence_id": existing(evidence_id("ligand_sites", res["uniprot"]))},
    })


@mcp.tool()
def structural_context_ranking(n: int = 10) -> list[dict]:
    """Candidate genes ranked by their best structural-context score, with the cell type that expresses them in the tumor.
    Describes monomer cavities and sampled cell attribution. It does not prioritize validated molecular glues,
    identify degrons or predict degradation. n <= 30."""
    rows = []
    ct = _ct()
    tum = ct[ct["tissue"] == "tumor"]
    for g, res in _structure().items():
        if not res["modeled"] or not res["pockets"]:
            continue
        s = tum[tum["hgnc_symbol"] == g]
        top = s.sort_values("share_of_tissue_transcripts", ascending=False).head(1)
        best = res["pockets"][0]
        rows.append(clean({
            "gene": g, "uniprot": res["uniprot"], "best_structural_context_score": best["mgd_score"],
            "best_pocket_volume_a3": best["volume_a3"], "accessible_lysines": len(best["accessible_lysines"]),
            "mean_plddt": res["mean_plddt"],
            "main_expressing_cell_group": None if top.empty else top.iloc[0]["cell_group"],
            "share_of_tumor_transcripts": None if top.empty else round(float(top.iloc[0]["share_of_tissue_transcripts"]), 3),
            "evidence_id": existing(evidence_id("structural_context", res["uniprot"], best["cavity"])),
        }))
    rows.sort(key=lambda r: -(r["best_structural_context_score"] or 0))
    return rows[:min(n, 30)]



@lru_cache(maxsize=1)
def _discovery() -> pd.DataFrame:
    return pd.read_parquet(PROCESSED / "discovery_candidates.parquet")


def _disc_row(r) -> dict:
    return clean({
        "gene": r.gene, "track": r.track, "status": r.status, "novelty_score": round(float(r.novelty_score), 3),
        "cell_group": r.source_group, "share_of_tumor_transcripts": round(float(r.source_share), 3),
        "rna_log2fc": round(float(r.rna_log2fc), 2), "protein_log2fc": round(float(r.protein_log2fc), 2),
        "ccrcc_papers": int(r.ccrcc_papers), "total_papers": int(r.total_papers),
        "clinical_drugs": r.clinical_drugs or None, "tractable_modalities": r.modalities or None,
        "surface_or_secreted": bool(r.surface_or_secreted),
        "vital_organ_max_ntpm": None if pd.isna(r.vital_max_ntpm) else round(float(r.vital_max_ntpm), 1),
        "vital_organ_tissue": None if pd.isna(r.vital_max_tissue) else r.vital_max_tissue,
        "grade_association": None if pd.isna(r.clin_rho) else f"{r.clin_variable} rho {r.clin_rho:.2f}",
        "rank_in_track": getattr(r, "rank_in_track", None),
        "tcga_log2fc": _f(getattr(r, "tcga_log2fc", None)), "replicated_in_tcga": _b(getattr(r, "replicated", None)),
        "rna_log2fc_at_purity1": _f(getattr(r, "rna_lfc_at_purity1", None)),
        "purity_robust": _b(getattr(r, "purity_robust", None)),
        "infiltration_driven": _b(getattr(r, "infiltration_driven", None)),
        "tcga_cox_log_hr": _f(getattr(r, "tcga_cox_log_hr", None)), "adverse_prognosis": _b(getattr(r, "adverse_prognosis", None)),
        "depmap_mean_effect_ccrcc": _f(getattr(r, "depmap_mean_effect_ccrcc", None)),
        "selective_dependency": _b(getattr(r, "selective_dependency", None)),
        "evidence_id": existing(evidence_id("discovery", r.gene)),
        "evidence_id_tcga": existing(evidence_id("tcga_de", r.gene)),
        "evidence_id_purity": existing(evidence_id("purity", r.gene)),
        "evidence_id_survival": existing(evidence_id("tcga_survival", r.gene))})


def _f(x, nd=3):
    return None if x is None or pd.isna(x) else round(float(x), nd)


def _b(x):
    return None if x is None or (isinstance(x, float) and pd.isna(x)) else bool(x)


@mcp.tool()
def novel_targets(track: str | None = None, status: str = "first-in-class opportunity", n: int = 10) -> list[dict]:
    """Under-studied, undrugged target candidates from the discovery screen, ranked by novelty score.

    The screen covers genes up in ccRCC tumors in both CPTAC RNA and protein that scRNA-seq attributes mostly to one cell
    type. track: 'tumor' (expressed by tumor cells), 'stroma' (endothelium / fibroblasts / normal epithelium) or 'immune';
    omit for all. status: 'first-in-class opportunity' (no clinical drug, <= 15 ccRCC papers, small-molecule or antibody
    tractable, low vital-organ expression, replicated in TCGA-KIRC and, for tumor / stroma genes, purity-robust),
    'novel, not replicated in TCGA', 'novel, not purity-robust (may be stromal / immune signal)',
    'novel, vital-organ safety unknown (no HPA data)', 'novel, off-tumor risk (vital-organ expression)',
    'novel, harder to drug',
    'has clinical drug (repurposing)', 'under-studied vs. established', or 'any' for no filter. n <= 30. These are
    data-ranked hypotheses: concordant tumor-up expression does not show that the gene is required for tumor growth.
    Ranked within the chosen track (the pool is mostly immune genes); each row carries the M14 layers: TCGA replication,
    log2FC at purity 1 (infiltration check), TCGA survival and DepMap dependency, with their evidence ids."""
    d = _discovery()
    if track:
        d = d[d["track"] == track]
    if status != "any":
        d = d[d["status"] == status]
    return [_disc_row(r) for r in d.head(min(n, 30)).itertuples()]


@mcp.tool()
def discovery_profile(gene: str) -> dict:
    """Discovery-screen profile of one gene: literature volume (ccRCC and overall), clinical drugs, tractability, location,
    vital-organ expression (Human Protein Atlas), grade association, novelty score and status."""
    g = canonical_gene(gene)
    d = _discovery()
    sub = d[d["gene"] == g]
    if sub.empty:
        return {"gene": g, "message": "not in the discovery pool (needs concordant tumor-up RNA and protein and a dominant cell type)"}
    return _disc_row(next(sub.itertuples()))


@mcp.tool()
def treatment_response(gene: str) -> dict:
    """Does the gene predict response to treatment in kidney cancer trials? IMmotion150 (263 patients: atezolizumab +
    bevacizumab, atezolizumab, sunitinib) per arm: log hazard ratio of progression per SD of expression (negative =
    longer benefit with higher expression), adjusted p-value and sample size, plus the difference between the
    combination and sunitinib, the 16-patient Choueiri 2016 nivolumab trial (direction only) and CheckMate 025
    (250 patients: nivolumab vs everolimus, per arm and the difference between them). No gene differs between the arms
    of either randomised trial after correction, so say 'associated with' and never 'predicts'. Cite evidence_id."""
    return applied_record(gene, "treatment_response", "gene not tested in the treatment-response trials (too little expression variation)")


@mcp.tool()
def subtype_profile(gene: str) -> dict:
    """Does the gene differ between tumour subgroups? Mutant minus wild-type log2 difference and adjusted p-value for
    VHL, PBRM1, BAP1 and SETD2 mutation (CPTAC RNA, CPTAC protein, TCGA-KIRC RNA) and the ccB minus ccA difference
    (ClearCode34 split) with a flag for subgroup target candidates. Immune genes can differ because of infiltration."""
    return applied_record(gene, "subtype_profile", "gene has no subtype comparison (not measured or too few tumours)")


@mcp.tool()
def adc_profile(gene: str) -> dict:
    """Suitability as an antibody-drug conjugate or radioligand target: surface topology (UniProt), tumour increase,
    staining in vital normal tissues (Human Protein Atlas immunohistochemistry, levels 0 none to 3 high), kidney staining,
    share of tumour cells expressing it in the weakest quarter of donors, the 0-1 suitability score and its rank."""
    return applied_record(gene, "adc_profile", "gene is not among the tumour-up candidates screened for ADC suitability")


@mcp.tool()
def immune_profile(gene: str) -> dict:
    """Immune-target profile: share of the gene's tumour transcripts made by immune cells, Spearman correlation across
    CPTAC tumours with CD8 T cells, regulatory T cells, M2 macrophages and the xCell immune score, log hazard ratio on
    atezolizumab + bevacizumab, and the immune-target score and rank. Checkpoint genes are controls."""
    return applied_record(gene, "immune_profile", "gene is not in the immune-profile focus set")


@mcp.tool()
def sc_replication_check(gene: str) -> dict:
    """Is the cell-type attribution the same in a second single-cell dataset (Zhang 2021, 7 ccRCC donors)? Dominant cell
    group and tumour-cell share in Li 2022 and in Zhang 2021 (seven shared groups), and the tumour vs proximal tubule
    log2 fold change in both."""
    return applied_record(gene, "sc_replication", "gene not measured in the second single-cell dataset")


if __name__ == "__main__":
    mcp.run()
