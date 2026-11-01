"""Target discovery: promising, under-studied and undrugged candidates (deterministic, no LLM).

Pool: genes up in tumor in BOTH CPTAC RNA and protein, that scRNA-seq attributes mostly to one cell group (tumor
cell = "tumor" track, immune cells = "immune", endothelium / fibroblast / normal epithelium = "stroma"). Each is annotated with
  - literature volume: Europe PMC records for the gene (and its previous / alias symbols) AND ccRCC terms, and overall,
  - drug status: Open Targets drug / clinical candidates and tractability modalities,
  - location: UniProt subcellular location (surface / secreted proteins are antibody-accessible),
  - safety: Human Protein Atlas expression in vital organs (on-target, off-tumor risk),
  - aggressiveness: correlation of tumor RNA / protein with histologic grade and pathologic stage (CPTAC clinical file).
Novelty score = weighted sum of components in [0, 1] (weights in config.yaml):
  omics (strength of concordant up-regulation), cell_specificity (share of tumour transcripts in the dominant
  cell group x support from tumor-vs-proximal-tubule), understudied (few ccRCC papers), undrugged (no clinical-stage
  drug), druggable (tractable modality + surface/secreted), safety (low expression in vital organs) and
  aggressiveness (positive, significant association with grade or stage).
A candidate is "novel" when it has no clinical-stage drug and <= novel_max_ccrcc_papers ccRCC papers; "first-in-class
opportunity" additionally needs small-molecule or antibody tractability and <= safe_max_ntpm nTPM in every vital organ.
Spatial and structure results (MGD pocket score) are attached to the shortlist as annotations, not scored.

v2 (M14) adds the independent M13 layers: TCGA-KIRC replication, purity extrapolation (does the up-regulation survive in
a pure tumour?), TCGA overall survival and DepMap dependency. They enter the novelty score, and "first-in-class" now
also requires replication and, for tumour / stroma genes, a purity-robust change. Immune-track genes keep their own
list (myeloid / lymphoid targets are legitimate), but those whose bulk change disappears at purity 1 are labelled
infiltration-driven. Rankings are reported within each track, because the pool is dominated by immune genes (the
scRNA dissociation over-represents immune cells, so they hold most tumour-compartment transcripts).
"""
from __future__ import annotations

import json
import math
import re
from datetime import date

import numpy as np
import pandas as pd

from .config import load
from .external import CCRCC_TERMS, europepmc_count, ot_target, uniprot
from .graph import GraphStore, evidence_id
from .paths import PROCESSED, RAW
from .singlecell import CELL_GROUPS

TODAY = date.today().isoformat()
SURFACE = re.compile(r"membrane|secreted|extracellular|cell surface", re.I)
IMMUNE = {"CD8_T", "CD4_T", "Treg", "T_other", "NK", "B", "Plasma", "Macrophage", "Monocyte", "cDC", "pDC", "Mast"}
CLINICAL = {"APPROVAL", "PHASE_4", "PHASE_3", "PHASE_2", "PHASE_1", "PHASE_2_3", "PHASE_1_2"}


def _gene_terms(symbol: str, hg: pd.DataFrame) -> str:
    """Europe PMC TITLE_ABS term for a symbol plus its previous/alias symbols that are specific enough."""
    names = {symbol}
    if symbol in hg.index:
        for col in ("prev_symbol", "alias_symbol"):
            v = hg.loc[symbol, col]
            for a in (v.split("|") if isinstance(v, str) else []):
                a = a.strip()
                if len(a) >= 4 and (re.search(r"\d", a) or len(a) >= 6) and a.isupper():
                    names.add(a)
    return "(" + " OR ".join(f'TITLE_ABS:"{n}"' for n in sorted(names)) + ")"


def pool() -> pd.DataFrame:
    cfg = load()["discovery"]
    u = pd.read_parquet(PROCESSED / "gene_universe.parquet").set_index("gene")
    u = u[(u.rna_padj < 0.05) & (u.rna_log2fc >= cfg["min_rna_log2fc"]) & (u.protein_padj < 0.05)
          & (u.protein_log2fc >= cfg["min_protein_log2fc"])]
    ct = pd.read_parquet(PROCESSED / "sc_celltype_expression.parquet")
    ct = ct[ct["tissue"] == "tumor"]
    top = (ct.sort_values("share_of_tissue_transcripts", ascending=False)
           .drop_duplicates("hgnc_symbol").set_index("hgnc_symbol"))
    tvp = pd.read_parquet(PROCESSED / "sc_tumor_vs_pt.parquet").drop_duplicates("hgnc_symbol").set_index("hgnc_symbol")
    rcc = ct[ct["cell_group"] == "RCC"].drop_duplicates("hgnc_symbol").set_index("hgnc_symbol")
    d = u.join(top[["cell_group", "pct_expressing", "share_of_tissue_transcripts"]]
               .rename(columns={"cell_group": "source_group", "pct_expressing": "source_pct",
                                "share_of_tissue_transcripts": "source_share"}), how="inner")
    d = d.join(rcc[["pct_expressing"]].rename(columns={"pct_expressing": "rcc_pct"}), how="left")
    d = d.join(tvp[["log2fc", "padj"]].rename(columns={"log2fc": "sc_tvp_log2fc", "padj": "sc_tvp_padj"}), how="left")
    d = d[(d.source_share >= cfg["min_source_share"]) & (d.source_pct >= cfg["min_source_pct"])]
    d["track"] = np.select([d["source_group"] == "RCC", d["source_group"].isin(IMMUNE)], ["tumor", "immune"], "stroma")
    return d.reset_index().rename(columns={"index": "gene"})


def annotate(d: pd.DataFrame) -> pd.DataFrame:
    hg = pd.read_csv(RAW / "hgnc_complete_set.txt", sep="\t", dtype=str,
                     usecols=["symbol", "ensembl_gene_id", "prev_symbol", "alias_symbol"]).dropna(subset=["symbol"])
    hg = hg.drop_duplicates("symbol").set_index("symbol")
    ot_assoc = pd.read_parquet(PROCESSED / "external_ot_associations.parquet").set_index("symbol")
    rows = []
    for g in d["gene"]:
        terms = _gene_terms(g, hg)
        rec = {"gene": g,
               "ccrcc_papers": europepmc_count(f"{terms} AND {CCRCC_TERMS}"),
               "total_papers": europepmc_count(terms)}
        ens = hg.loc[g, "ensembl_gene_id"] if g in hg.index else None
        t = ot_target(ens) if isinstance(ens, str) else None
        drugs = (t or {}).get("drugs", [])
        rec.update({"ot_available": t is not None,
                    "clinical_drugs": [x["name"] for x in drugs if x["max_stage"] in CLINICAL],
                    "any_drugs": [x["name"] for x in drugs],
                    "max_drug_stage": max((x["max_stage"] for x in drugs), key=lambda s: (s in CLINICAL, s), default=None),
                    "modalities": (t or {}).get("tractable_modalities", []),
                    "ot_ccrcc_score": float(ot_assoc["score"].get(g, 0.0)),
                    "ot_genetic": float(ot_assoc["dt_genetic_association"].get(g, 0.0))})
        u = uniprot(g)
        loc = (u or {}).get("subcellular_location", [])
        rec.update({"uniprot": (u or {}).get("accession"), "protein_name": (u or {}).get("name"),
                    "locations": loc, "surface_or_secreted": any(SURFACE.search(x) for x in loc)})
        rec["ensembl"] = ens if isinstance(ens, str) else None
        rows.append(rec)
    out = d.merge(pd.DataFrame(rows), on="gene")
    ntis = pd.read_parquet(PROCESSED / "normal_tissue.parquet")
    out = out.merge(ntis, on="ensembl", how="left")
    clin = pd.read_parquet(PROCESSED / "clinical_associations.parquet")
    clin = clin[clin.padj < 0.1].sort_values("rho", ascending=False).drop_duplicates("hgnc_symbol")
    clin = clin.set_index("hgnc_symbol")
    out["clin_rho"] = out["gene"].map(clin["rho"])
    out["clin_padj"] = out["gene"].map(clin["padj"])
    out["clin_variable"] = out["gene"].map(clin["variable"] + " (" + clin["layer"] + ")")
    return out


def annotate_m13(d: pd.DataFrame) -> pd.DataFrame:
    """Attach the M13 layers (computed by layers.py) to the discovery pool. Missing layers leave NaN."""
    d = d.copy().set_index("gene")

    def take(name, cols, key="gene"):
        p = PROCESSED / name
        if not p.exists():
            for c in cols.values():
                d[c] = np.nan
            return
        t = pd.read_parquet(p).dropna(subset=[key]).drop_duplicates(key).set_index(key)
        for src, dst in cols.items():
            d[dst] = d.index.map(t[src]) if src in t.columns else np.nan

    take("layer_tcga.parquet", {"log2fc": "tcga_log2fc", "padj": "tcga_padj", "cox_log_hr_per_sd": "tcga_cox_log_hr",
                                 "cox_padj": "tcga_cox_padj"})
    take("layer_purity.parquet", {"rna_lfc_at_purity1": "rna_lfc_at_purity1", "rna_lfc_at_purity1_se": "rna_lfc_at_purity1_se",
                                  "rna_purity_rho": "rna_purity_rho"}, key="hgnc_symbol")
    take("layer_dependency.parquet", {"mean_effect_ccrcc": "depmap_mean_effect_ccrcc", "selectivity": "depmap_selectivity",
                                      "frac_ccrcc_dependent": "depmap_frac_dependent", "padj": "depmap_padj",
                                      "common_essential": "depmap_common_essential", "n_ccrcc_lines": "depmap_n_lines"})
    return d.reset_index().rename(columns={"index": "gene"})


M13_COLS = ("tcga_log2fc", "tcga_padj", "tcga_cox_log_hr", "tcga_cox_padj", "rna_lfc_at_purity1", "rna_lfc_at_purity1_se",
            "rna_purity_rho", "depmap_mean_effect_ccrcc", "depmap_selectivity", "depmap_frac_dependent", "depmap_padj",
            "depmap_common_essential", "depmap_n_lines")


def m13_flags(d: pd.DataFrame) -> pd.DataFrame:
    """Boolean gates derived from the M13 values (missing values -> False)."""
    d = d.copy()
    for c in M13_COLS:
        if c not in d.columns:
            d[c] = np.nan
    d["replicated"] = (d["tcga_padj"] < 0.05) & (d["tcga_log2fc"] > 0)
    d["purity_robust"] = (d["rna_lfc_at_purity1"] - 1.96 * d["rna_lfc_at_purity1_se"]) > 0
    d["adverse_prognosis"] = (d["tcga_cox_padj"] < 0.05) & (d["tcga_cox_log_hr"] > 0)
    dc = load()["convergence"]["dependency"]
    ess = d["depmap_common_essential"].astype("boolean").fillna(False).astype(bool)
    d["selective_dependency"] = ((d["depmap_padj"] < dc["padj"]) & (d["depmap_selectivity"] > 0)
                                 & (d["depmap_frac_dependent"] > 0) & ~ess)
    return d


def score(d: pd.DataFrame) -> pd.DataFrame:
    cfg = load()["discovery"]
    w = cfg["weights"]
    d = d.copy()
    if "tcga_log2fc" not in d.columns:
        d = annotate_m13(d)
    if "track" not in d.columns:
        d["track"] = "tumor"
    d = m13_flags(d)
    omics = 0.5 * (d.protein_log2fc / 2).clip(0, 1) + 0.5 * (d.rna_log2fc / 4).clip(0, 1)
    sc = (d.sc_tvp_log2fc.fillna(0) / 4).clip(0, 1) * np.where(d.sc_tvp_padj.fillna(1) < 0.1, 1.0, 0.5)
    d["c_omics"] = omics
    d["c_cell_specificity"] = 0.6 * d.source_share.clip(upper=1) + 0.4 * sc
    d["c_understudied"] = 1 - (np.log10(1 + d.ccrcc_papers) / math.log10(1 + cfg["saturation_papers"])).clip(upper=1)
    d["c_undrugged"] = (d.clinical_drugs.map(len) == 0).astype(float)
    # Open Targets modality codes: SM small molecule, AB antibody, PR PROTAC/degrader, OC other clinical. PR is called for
    # ~90 % of genes, so it only adds a little; SM/AB tractability and surface/secreted location carry the component.
    core = d.modalities.map(lambda m: bool({"SM", "AB"} & set(m)))
    degrader = d.modalities.map(lambda m: "PR" in set(m))
    d["c_druggable"] = 0.5 * core + 0.2 * degrader + 0.3 * d.surface_or_secreted
    d["druggable_core"] = core
    d["c_safety"] = (1 - (d.vital_max_ntpm / 50).clip(0, 1)).fillna(0.5)
    d["c_aggressiveness"] = (d.clin_rho.fillna(0) / 0.4).clip(0, 1)
    # ---- M14: independent layers
    caps = load()["convergence"]["omics_caps"]
    d["c_replication"] = ((d.tcga_log2fc / caps["rna"]).clip(0, 1) * d.replicated.fillna(False)).fillna(0.0)
    gated = d.track.isin(cfg.get("require_purity_robust_tracks", ["tumor", "stroma"]))
    intrinsic = ((d.rna_lfc_at_purity1 / caps["rna"]).clip(0, 1) * d.purity_robust.fillna(False)).fillna(0.0)
    # immune-track genes are scored on the other components only (the purity component would penalise every one)
    d["c_tumor_intrinsic"] = np.where(gated, intrinsic, 0.0)
    d["c_prognosis"] = ((d.tcga_cox_log_hr / 0.5).clip(0, 1) * d.adverse_prognosis.fillna(False)).fillna(0.0)
    d["c_dependency"] = ((d.depmap_selectivity / 0.5).clip(0, 1) * d.selective_dependency.fillna(False)).fillna(0.0)
    wsum = {k: v for k, v in w.items() if not (k == "tumor_intrinsic")}
    base = sum(wt * d[f"c_{k}"] for k, wt in wsum.items())
    # immune track: same score without the purity component (normalised by its own weight total, so tracks compare)
    d["novelty_score"] = np.where(gated, (base + w.get("tumor_intrinsic", 0) * d["c_tumor_intrinsic"]) / sum(w.values()),
                                  base / sum(wsum.values()))
    d["novel"] = (d.clinical_drugs.map(len) == 0) & (d.ccrcc_papers <= cfg["novel_max_ccrcc_papers"])
    safe = d.vital_max_ntpm.fillna(0) <= cfg["safe_max_ntpm"]
    d["safe"] = safe
    repl_ok = d.replicated.fillna(False) if cfg.get("require_replication", True) else True
    purity_ok = np.where(gated, d.purity_robust.fillna(False), True)
    d["infiltration_driven"] = (d.track == "immune") & ~d.purity_robust.fillna(False) & (d.rna_purity_rho < 0)
    fic = d.novel & d.druggable_core & safe & repl_ok & purity_ok
    d["status"] = np.select(
        [d.clinical_drugs.map(len) > 0, fic, d.novel & d.druggable_core & safe & ~repl_ok,
         d.novel & d.druggable_core & safe & ~purity_ok, d.novel & d.druggable_core, d.novel],
        ["has clinical drug (repurposing)", "first-in-class opportunity", "novel, not replicated in TCGA",
         "novel, not purity-robust (may be stromal / immune signal)", "novel, off-tumor risk (vital-organ expression)",
         "novel, harder to drug"],
        default="under-studied vs. established")
    d["rank_in_track"] = d.groupby("track")["novelty_score"].rank(ascending=False, method="first").astype(int)
    return d.sort_values("novelty_score", ascending=False).reset_index(drop=True)


def run() -> dict:
    cfg = load()["discovery"]
    d = score(annotate(pool()))
    keep = d.drop(columns=["clinical_drugs", "any_drugs", "modalities", "locations"]).assign(
        clinical_drugs=d.clinical_drugs.map("; ".join), any_drugs=d.any_drugs.map("; ".join),
        modalities=d.modalities.map("; ".join), locations=d.locations.map("; ".join))
    keep.to_parquet(PROCESSED / "discovery_candidates.parquet", index=False)
    short = []
    for track, n in cfg["shortlist_per_track"].items():
        short += d[d.novel & d.druggable_core & (d.track == track)].head(n)["gene"].tolist()
    (PROCESSED / "discovery_shortlist.json").write_text(json.dumps(short))
    fic = d[d.status == "first-in-class opportunity"]
    summary = {"pool": int(len(d)), "novel": int(d.novel.sum()), "shortlist": short,
               "by_status": d.status.value_counts().to_dict(), "by_track": d.track.value_counts().to_dict(),
               "first_in_class_by_track": fic.track.value_counts().to_dict(),
               "top_by_track": {t: g[g.status == "first-in-class opportunity"].head(5)["gene"].tolist()
                                for t, g in d.groupby("track")},
               "m13": {"replicated": int(d.replicated.fillna(False).sum()),
                       "purity_robust": int(d.purity_robust.fillna(False).sum()),
                       "adverse_prognosis": int(d.adverse_prognosis.fillna(False).sum()),
                       "selective_dependency": int(d.selective_dependency.fillna(False).sum()),
                       "infiltration_driven": int(d.infiltration_driven.sum())},
               "thresholds": {k: cfg[k] for k in ("min_rna_log2fc", "min_protein_log2fc", "min_source_share",
                                                   "min_source_pct", "novel_max_ccrcc_papers", "saturation_papers")},
               "weights": cfg["weights"]}
    (PROCESSED / "discovery_summary.json").write_text(json.dumps(summary, indent=1))
    return summary


def load_graph(store: GraphStore | None = None) -> dict:
    """Evidence + claims for every candidate in the pool."""
    store = store or GraphStore()
    d = pd.read_parquet(PROCESSED / "discovery_candidates.parquet")
    gene_ids = {r["name"]: r["id"] for r in store.conn.execute("SELECT id, name FROM nodes WHERE type = 'Gene'")}
    n = 0
    with store.batch():
        store.add_node("DATASET:EUROPEPMC", "Dataset", "Europe PMC literature counts", source="Europe PMC REST API")
        for r in d.itertuples():
            gid = gene_ids.get(r.gene)
            if gid is None:
                continue
            ev = store.add_evidence(
                evidence_id("discovery", r.gene), "discovery_profile", "Europe PMC + Open Targets + UniProt + CPTAC + scRNA-seq",
                TODAY, gene=r.gene, track=r.track, status=r.status, novelty_score=round(float(r.novelty_score), 4),
                ccrcc_papers=int(r.ccrcc_papers), total_papers=int(r.total_papers),
                clinical_drugs=r.clinical_drugs or None, max_drug_stage=r.max_drug_stage, modalities=r.modalities,
                surface_or_secreted=bool(r.surface_or_secreted), locations=r.locations,
                source_group=r.source_group, source_share=round(float(r.source_share), 3),
                source_pct=round(float(r.source_pct), 1), rna_log2fc=round(float(r.rna_log2fc), 3),
                protein_log2fc=round(float(r.protein_log2fc), 3),
                components={k: round(float(getattr(r, f"c_{k}")), 3) for k in
                            ("omics", "cell_specificity", "understudied", "undrugged", "druggable", "safety",
                             "aggressiveness", "replication", "tumor_intrinsic", "prognosis", "dependency")},
                rank_in_track=int(r.rank_in_track),
                tcga_log2fc=None if pd.isna(r.tcga_log2fc) else round(float(r.tcga_log2fc), 3),
                tcga_padj=None if pd.isna(r.tcga_padj) else float(f"{r.tcga_padj:.3g}"),
                tcga_cox_log_hr=None if pd.isna(r.tcga_cox_log_hr) else round(float(r.tcga_cox_log_hr), 3),
                tcga_cox_padj=None if pd.isna(r.tcga_cox_padj) else float(f"{r.tcga_cox_padj:.3g}"),
                rna_lfc_at_purity1=None if pd.isna(r.rna_lfc_at_purity1) else round(float(r.rna_lfc_at_purity1), 3),
                rna_lfc_at_purity1_se=None if pd.isna(r.rna_lfc_at_purity1_se) else round(float(r.rna_lfc_at_purity1_se), 3),
                depmap_mean_effect_ccrcc=None if pd.isna(r.depmap_mean_effect_ccrcc) else round(float(r.depmap_mean_effect_ccrcc), 3),
                replicated=bool(r.replicated) if pd.notna(r.replicated) else None,
                purity_robust=bool(r.purity_robust) if pd.notna(r.purity_robust) else None,
                infiltration_driven=bool(r.infiltration_driven),
                vital_max_ntpm=None if pd.isna(r.vital_max_ntpm) else round(float(r.vital_max_ntpm), 1),
                vital_max_tissue=None if pd.isna(r.vital_max_tissue) else r.vital_max_tissue,
                clinical_association=None if pd.isna(r.clin_rho) else f"{r.clin_variable} rho {r.clin_rho:.2f}",
                url=f"https://europepmc.org/search?query={r.gene}%20AND%20%22renal%20cell%20carcinoma%22")
            store.add_claim(gid, "novel_candidate_for" if r.novel else "discovery_scored_for", "MONDO:0005005", ev)
            n += 1
    return {"discovery_evidence": n}


if __name__ == "__main__":
    s = run()
    print(json.dumps({k: v for k, v in s.items() if k != "weights"}, indent=1))
