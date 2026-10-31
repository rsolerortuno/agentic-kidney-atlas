"""Spatial context (Li et al. 2022 10x Visium, tumour core + tumour-normal interface sections).

Deterministic, no LLM:
1. Deconvolution: each spot's CPM profile over marker genes is fit by non-negative least squares to the
   scRNA-seq cell-group profiles of the same study (singlecell.py). The weights, normalized to 1, estimate the
   share of the spot's transcripts coming from each cell group. Candidate targets are never used as markers,
   so target expression is an independent read-out of the deconvolution.
2. Who expresses the target in situ: NNLS of the target's spot CPM on the spot compositions (beta_k = in situ
   expression per unit of group k; contribution_k = beta_k * mean share_k).
3. Who the target-expressing spots are next to, weighted by abundance: target-high spots are the top decile of
   a section (>= 1 UMI). For every cell group, the mean composition of their ring-1 neighbours (the 6 adjacent
   spots, ~100 um centre to centre) is compared with the expectation for random spots of the same sections
   (observed / expected, so common cell types are not favoured). Section-stratified permutations give
   two-sided empirical p-values; Benjamini-Hochberg across targets x groups.
4. Composition-corrected source of the bulk signal: share of a gene's transcripts expected from each group
   = in situ group share (tumour-core sections) x scRNA profile, which removes the dissociation bias of
   scRNA-seq (immune cells over-captured, tumour cells under-captured).
"""
from __future__ import annotations

import json

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.optimize import nnls
from scipy.spatial import cKDTree
from statsmodels.stats.multitest import multipletests

from .config import load
from .paths import PROCESSED, RAW

FILES = {"core": RAW / "singlecell" / "li2022_visium_core.h5ad",
         "interface": RAW / "singlecell" / "li2022_visium_interface.h5ad"}
# Visium spots hold several cells: T/NK and DC subsets are merged into groups the spot data can separate
SPATIAL_GROUPS = {
    "RCC": ["RCC"], "PT": ["PT"], "Epi_other": ["Epi_other"], "Endothelial": ["Endothelial"],
    "Fibroblast": ["Fibroblast"], "Macrophage": ["Macrophage"], "Monocyte": ["Monocyte"], "DC": ["cDC", "pDC"],
    "Mast": ["Mast"], "T_cell": ["CD8_T", "CD4_T", "Treg", "T_other"], "NK": ["NK"], "B": ["B"], "Plasma": ["Plasma"],
}


def reference() -> pd.DataFrame:
    """Merged cell-group CPM (genes x spatial groups, Ensembl index)."""
    ref = pd.read_parquet(PROCESSED / "sc_reference.parquet")
    sym = ref.pop("hgnc_symbol")
    merged = pd.DataFrame({g: ref[parts].sum(axis=1) for g, parts in SPATIAL_GROUPS.items()})
    cpm = merged / merged.sum(axis=0) * 1e6
    cpm.insert(0, "hgnc_symbol", sym)
    return cpm


def select_markers(ref: pd.DataFrame, exclude: set[str], genes_in_spots: set[str]) -> dict[str, list[str]]:
    cfg = load()["cell_context"]["spatial"]
    x = ref.drop(columns="hgnc_symbol")
    ok = ~ref["hgnc_symbol"].isin(exclude) & ref.index.isin(genes_in_spots) & ref["hgnc_symbol"].notna()
    ok &= ~ref["hgnc_symbol"].fillna("").str.match(r"^(MT-|RPL|RPS)")
    x = x[ok]
    out = {}
    for g in x.columns:
        others = x.drop(columns=g).max(axis=1)
        score = np.log2((x[g] + 1) / (others + 1))
        cand = score[(x[g] >= cfg["min_marker_cpm"]) & (score > 1)].sort_values(ascending=False)
        out[g] = cand.head(cfg["markers_per_group"]).index.tolist()
    return out


def load_spots() -> tuple[pd.DataFrame, sp.csr_matrix, pd.Index]:
    obs, mats, var = [], [], None
    for region, path in FILES.items():
        a = ad.read_h5ad(path)
        var = a.var_names if var is None else var
        assert (a.var_names == var).all()
        o = pd.DataFrame({"spot": a.obs_names, "sample": a.obs["sample"].astype(str).values,
                          "donor": a.obs["donor_id"].astype(str).values, "region": region,
                          "x": a.obsm["spatial"][:, 0], "y": a.obsm["spatial"][:, 1]})
        obs.append(o)
        mats.append(sp.csr_matrix(a.X, dtype=np.float64))
    o = pd.concat(obs, ignore_index=True)
    o["section"] = o["region"] + ":" + o["sample"]
    return o, sp.vstack(mats).tocsr(), var


def deconvolve(counts: sp.csr_matrix, var: pd.Index, ref: pd.DataFrame, markers: dict[str, list[str]]) -> np.ndarray:
    genes = list(dict.fromkeys(g for m in markers.values() for g in m))
    col = pd.Index(var).get_indexer(genes)
    lib = np.asarray(counts.sum(axis=1)).ravel()
    y = (sp.diags(1e6 / np.maximum(lib, 1)) @ counts[:, col]).toarray()
    R = ref.loc[genes].drop(columns="hgnc_symbol").to_numpy()
    scale = R.mean(axis=1, keepdims=True)          # equalize marker genes so a few abundant ones do not dominate
    Rs, ys = R / scale, y / scale.T
    W = np.vstack([nnls(Rs, ys[i])[0] for i in range(len(ys))])
    s = W.sum(axis=1, keepdims=True)
    return np.divide(W, s, out=np.zeros_like(W), where=s > 0)


def neighbor_matrix(spots: pd.DataFrame) -> sp.csr_matrix:
    """Row-normalized ring-1 adjacency within each section (centre distance <= 1.25 x median NN distance)."""
    rows, cols = [], []
    for _, idx in spots.groupby("section").indices.items():
        xy = spots[["x", "y"]].to_numpy()[idx]
        tree = cKDTree(xy)
        d, _ = tree.query(xy, k=2)
        r = 1.25 * np.median(d[:, 1])
        for i, nb in enumerate(tree.query_ball_point(xy, r)):
            for j in nb:
                if j != i:
                    rows.append(idx[i])
                    cols.append(idx[j])
    A = sp.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(spots), len(spots)))
    deg = np.asarray(A.sum(axis=1)).ravel()
    return sp.diags(1 / np.maximum(deg, 1)) @ A


def target_expression(counts: sp.csr_matrix, var: pd.Index, ens: str) -> tuple[np.ndarray, np.ndarray]:
    j = pd.Index(var).get_loc(ens)
    lib = np.asarray(counts.sum(axis=1)).ravel()
    umi = counts[:, j].toarray().ravel()
    return umi, umi / np.maximum(lib, 1) * 1e6


def neighborhood_enrichment(high: np.ndarray, sections: np.ndarray, NP: np.ndarray, has_nb: np.ndarray,
                            rng: np.random.Generator, n_perm: int) -> dict:
    """Observed vs section-stratified expected neighbour composition of target-high spots."""
    obs_sum, exp_sum, null = np.zeros(NP.shape[1]), np.zeros(NP.shape[1]), np.zeros((n_perm, NP.shape[1]))
    n_high = 0
    for s in np.unique(sections[high]):
        pool = np.flatnonzero((sections == s) & has_nb)
        h = np.flatnonzero(high & (sections == s) & has_nb)
        if len(h) == 0:
            continue
        n_high += len(h)
        obs_sum += NP[h].sum(axis=0)
        exp_sum += len(h) * NP[pool].mean(axis=0)
        draws = np.vstack([rng.choice(pool, size=len(h), replace=False) for _ in range(n_perm)])
        null += NP[draws].sum(axis=1)
    obs, exp, null = obs_sum / n_high, exp_sum / n_high, null / n_high
    eps = 1e-4
    lr = np.log2((obs + eps) / (exp + eps))
    lr_null = np.log2((null + eps) / (exp + eps))
    p = (1 + (np.abs(lr_null) >= np.abs(lr)).sum(axis=0)) / (1 + n_perm)
    return {"observed": obs, "expected": exp, "log2_enrichment": lr, "p": p, "n_high": n_high}


def run() -> dict:
    cfg_all = load()["cell_context"]
    cfg = cfg_all["spatial"]
    targets = cfg_all["targets"]
    rng = np.random.default_rng(cfg["seed"])
    ref = reference()
    spots, counts, var = load_spots()
    markers = select_markers(ref, set(targets), set(var))
    groups = list(SPATIAL_GROUPS)
    P = deconvolve(counts, var, ref, markers)
    dec = spots.copy()
    for k, g in enumerate(groups):
        dec[g] = P[:, k].astype("float32")
    dec.to_parquet(PROCESSED / "spatial_deconvolution.parquet", index=False)

    A = neighbor_matrix(spots)
    NP = A @ P
    has_nb = np.asarray((A > 0).sum(axis=1)).ravel() > 0
    sections = spots["section"].to_numpy()
    sym2ens = ref["hgnc_symbol"].dropna().reset_index().drop_duplicates("hgnc_symbol").set_index("hgnc_symbol")
    sym2ens = sym2ens.iloc[:, 0].to_dict()
    core = (spots["region"] == "core").to_numpy()
    mean_share_core = P[core].mean(axis=0)
    ref_groups = ref.drop(columns="hgnc_symbol")

    nb_rows, ex_rows, summary_t = [], [], []
    for gene in targets:
        ens = sym2ens.get(gene)
        if ens is None or ens not in var:
            summary_t.append({"gene": gene, "measured": False})
            continue
        umi, cpm = target_expression(counts, var, ens)
        high = np.zeros(len(spots), bool)
        used_sections = []
        for s, idx in spots.groupby("section").indices.items():
            pos = idx[umi[idx] > 0]
            if len(pos) < cfg["min_high_spots"]:
                continue
            thr = np.quantile(cpm[idx], cfg["high_quantile"])
            h = idx[(cpm[idx] >= thr) & (umi[idx] > 0)]
            if len(h) >= cfg["min_high_spots"]:
                high[h] = True
                used_sections.append(s)
        beta = nnls(P, cpm)[0]
        contrib = beta * P.mean(axis=0)
        contrib = contrib / contrib.sum() if contrib.sum() > 0 else contrib
        r = ref_groups.loc[ens].to_numpy()
        corrected = mean_share_core * r / (mean_share_core * r).sum() if (mean_share_core * r).sum() > 0 else r * 0
        for k, g in enumerate(groups):
            ex_rows.append({"gene": gene, "cell_group": g, "in_situ_beta_cpm": float(beta[k]),
                            "in_situ_contribution": float(contrib[k]), "composition_corrected_share": float(corrected[k]),
                            "spearman_with_share": float(pd.Series(cpm).corr(pd.Series(P[:, k]), method="spearman"))})
        entry = {"gene": gene, "measured": True, "spots_with_umi": int((umi > 0).sum()),
                 "sections_used": len(used_sections), "n_high_spots": int(high.sum()),
                 "top_in_situ_expresser": groups[int(np.argmax(contrib))] if contrib.sum() > 0 else None,
                 "top_corrected_source": groups[int(np.argmax(corrected))] if corrected.sum() > 0 else None}
        if used_sections:
            e = neighborhood_enrichment(high, sections, NP, has_nb, rng, cfg["permutations"])
            for k, g in enumerate(groups):
                nb_rows.append({"gene": gene, "cell_group": g, "observed_neighbor_share": float(e["observed"][k]),
                                "expected_neighbor_share": float(e["expected"][k]),
                                "log2_enrichment": float(e["log2_enrichment"][k]), "pvalue": float(e["p"][k]),
                                "n_high_spots": int(e["n_high"]), "n_sections": len(used_sections)})
        summary_t.append(entry)

    nb = pd.DataFrame(nb_rows)
    nb["padj"] = multipletests(nb["pvalue"], method="fdr_bh")[1] if len(nb) else []
    nb.to_parquet(PROCESSED / "spatial_neighbors.parquet", index=False)
    ex = pd.DataFrame(ex_rows)
    ex.to_parquet(PROCESSED / "spatial_expressers.parquet", index=False)
    for t in summary_t:
        if not t["measured"]:
            continue
        s = nb[nb["gene"] == t["gene"]]
        t["enriched_neighbors"] = [{"cell_group": r.cell_group, "log2_enrichment": round(r.log2_enrichment, 3),
                                    "padj": float(f"{r.padj:.3g}")}
                                   for r in s[(s["padj"] < 0.05) & (s["log2_enrichment"] > 0)]
                                   .sort_values("log2_enrichment", ascending=False).itertuples()]

    def corr(gene, group):
        ens = sym2ens.get(gene)
        if ens is None or ens not in var:
            return None
        _, cpm = target_expression(counts, var, ens)
        return round(float(pd.Series(cpm).corr(pd.Series(P[:, groups.index(group)]), method="spearman")), 3)

    validation = {
        "mean_share_core": dict(zip(groups, np.round(mean_share_core, 3).tolist())),
        "mean_share_interface": dict(zip(groups, np.round(P[~core].mean(axis=0), 3).tolist())),
        # targets were excluded from the markers: their correlation with the matching group is an independent check
        "spearman_CA9_vs_RCC": corr("CA9", "RCC"), "spearman_KDR_vs_Endothelial": corr("KDR", "Endothelial"),
        "spearman_CD74_vs_Macrophage": corr("CD74", "Macrophage"), "spearman_PTPRC_vs_T_cell": corr("PTPRC", "T_cell"),
    }
    summary = {"spots": int(len(spots)), "sections": int(spots["section"].nunique()),
               "donors": int(spots["donor"].nunique()), "groups": groups,
               "markers": {g: len(m) for g, m in markers.items()}, "marker_genes": markers,
               "validation": validation, "targets": summary_t}
    (PROCESSED / "spatial_summary.json").write_text(json.dumps(summary, indent=1))
    return summary


if __name__ == "__main__":
    s = run()
    print(json.dumps({k: v for k, v in s.items() if k not in ("marker_genes", "targets")}, indent=1))
    for t in s["targets"]:
        print(t["gene"], t.get("top_in_situ_expresser"), t.get("top_corrected_source"), t.get("n_high_spots"),
              [(e["cell_group"], e["log2_enrichment"]) for e in t.get("enriched_neighbors", [])][:5])
