"""Replicate the single-cell attribution in a second, independent ccRCC dataset.

Deterministic, no LLM. The cell-type calls (which cell group holds a gene's tumour transcripts; tumour cell vs proximal
tubule) come from one study, Li et al. 2022. Here they are recomputed in Zhang et al. 2021 (PNAS; CZ CELLxGENE):
7 ccRCC donors (20,509 cells) and 5 normal-kidney donors (6,044 cells), 4 donors with both.

Cell groups differ between the studies, so both are mapped to seven shared groups: RCC, Endothelial, Stroma (fibroblast,
pericyte, vascular smooth muscle), Myeloid (macrophage, monocyte, dendritic), T/NK, B/Plasma and Mast. For every gene:
- share of its tumour-tissue transcripts per shared group, in each study, and whether the dominant group agrees;
- tumour cell vs proximal tubule log2 fold change of donor pseudobulk CPM (Zhang: paired over the 4 shared donors).
The summary reports agreement for the discovery pool by track and for the configured report targets.
"""
from __future__ import annotations

import json

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats

from .config import load
from .omics import load_hgnc
from .paths import PROCESSED, RAW

SC = RAW / "singlecell"
SHARED = ["RCC", "Endothelial", "Stroma", "Myeloid", "T_NK", "B_Plasma", "Mast"]
LI_TO_SHARED = {"RCC": "RCC", "Endothelial": "Endothelial", "Fibroblast": "Stroma", "Macrophage": "Myeloid",
                "Monocyte": "Myeloid", "cDC": "Myeloid", "pDC": "Myeloid", "CD8_T": "T_NK", "CD4_T": "T_NK",
                "Treg": "T_NK", "T_other": "T_NK", "NK": "T_NK", "B": "B_Plasma", "Plasma": "B_Plasma", "Mast": "Mast"}
ZHANG_TO_SHARED = {"Tumor": "RCC", "Macro": "Myeloid", "Macro_MKI67": "Myeloid", "Mono": "Myeloid",
                   "Endo_PLVAP": "Endothelial", "Endo_ACKR1": "Endothelial", "vSMC": "Stroma", "Peri": "Stroma",
                   "Tcell": "T_NK", "Tcell_CD8": "T_NK", "NKcell": "T_NK", "Bcell": "B_Plasma", "Plasma": "B_Plasma",
                   "Mast": "Mast"}


def _symbols(var_names: pd.Index) -> list:
    hg = load_hgnc().dropna(subset=["ensembl_gene_id"]).drop_duplicates("ensembl_gene_id")
    sym = dict(zip(hg["ensembl_gene_id"], hg["symbol"]))
    return [sym.get(v) for v in var_names]


def _load(name: str) -> tuple[sp.csr_matrix, pd.DataFrame, list]:
    a = ad.read_h5ad(SC / name)
    X = sp.csr_matrix(a.raw.X if a.raw is not None else a.X)
    obs = pd.DataFrame({"donor": a.obs["donor_id"].astype(str).values,
                        "group": a.obs["author_cell_type"].astype(str).values})
    return X, obs, _symbols(a.raw.var_names if a.raw is not None else a.var_names)


def zhang_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    X, obs, sym = _load("zhang2021_ccrcc.h5ad")
    obs["shared"] = obs["group"].map(ZHANG_TO_SHARED)
    keep = obs["shared"].notna().to_numpy()
    G = pd.get_dummies(obs.loc[keep, "shared"]).reindex(columns=SHARED, fill_value=False).to_numpy(float)
    sums = (sp.csr_matrix(G.T) @ X[keep]).toarray()          # shared group x gene, tumour tissue
    share = pd.DataFrame(sums / np.maximum(sums.sum(axis=0, keepdims=True), 1), index=SHARED, columns=sym).T
    share = share[share.index.notna()]
    share = share.groupby(level=0).sum()
    # tumour cells vs proximal tubule, paired by donor
    Xn, on, symn = _load("zhang2021_normal.h5ad")
    assert symn == sym
    def pb(Xm, o, mask):
        out = {}
        for d in o.loc[mask, "donor"].unique():
            m = (mask & (o["donor"] == d)).to_numpy()
            if m.sum() >= 20:
                v = np.asarray(Xm[m].sum(axis=0)).ravel()
                out[d] = np.log2(v / v.sum() * 1e6 + 1)
        return out
    tum = pb(X, obs, obs["group"].eq("Tumor"))
    pt = pb(Xn, on, on["group"].str.startswith("PT"))
    donors = sorted(set(tum) & set(pt))
    D = np.array([tum[d] - pt[d] for d in donors])
    t, p = stats.ttest_1samp(D, 0, axis=0) if len(donors) >= 3 else (np.full(D.shape[1], np.nan),) * 2
    tvp = pd.DataFrame({"gene": sym, "zhang_tvp_log2fc": D.mean(axis=0), "zhang_tvp_pvalue": p,
                        "zhang_tvp_donors": len(donors)}).dropna(subset=["gene"]).drop_duplicates("gene").set_index("gene")
    return share, tvp


def li_share() -> pd.DataFrame:
    ct = pd.read_parquet(PROCESSED / "sc_celltype_expression.parquet")
    ct = ct[ct["tissue"] == "tumor"].copy()
    ct["shared"] = ct["cell_group"].map(LI_TO_SHARED)
    s = ct.groupby(["hgnc_symbol", "shared"])["share_of_tissue_transcripts"].sum().unstack().reindex(columns=SHARED)
    return s.fillna(0.0)


def run() -> dict:
    zs, ztvp = zhang_tables()
    ls = li_share()
    genes = ls.index.intersection(zs.index)
    expressed = zs.loc[genes].sum(axis=1) > 0
    genes = genes[expressed.to_numpy()]
    L, Z = ls.loc[genes], zs.loc[genes]
    tvp = pd.read_parquet(PROCESSED / "sc_tumor_vs_pt.parquet").drop_duplicates("hgnc_symbol").set_index("hgnc_symbol")
    out = pd.DataFrame({"gene": genes, "li_top": L.idxmax(axis=1).values, "zhang_top": Z.idxmax(axis=1).values,
                        "li_rcc_share": L["RCC"].values, "zhang_rcc_share": Z["RCC"].values,
                        "li_top_share": L.max(axis=1).values, "zhang_top_share": Z.max(axis=1).values})
    out["agree_top"] = out["li_top"] == out["zhang_top"]
    rho = [stats.spearmanr(L.loc[g], Z.loc[g])[0] if L.loc[g].std() > 0 and Z.loc[g].std() > 0 else np.nan for g in genes]
    out["share_profile_rho"] = rho
    out = out.set_index("gene").join(ztvp, how="left")
    out["li_tvp_log2fc"] = tvp["log2fc"].reindex(out.index)
    out.index.name = "gene"
    out = out.reset_index()
    out.to_parquet(PROCESSED / "sc_replication.parquet", index=False)

    disc = pd.read_parquet(PROCESSED / "discovery_candidates.parquet", columns=["gene", "track"])
    track_map = {"tumor": "RCC", "stroma": {"Endothelial", "Stroma"}, "immune": {"Myeloid", "T_NK", "B_Plasma", "Mast"}}
    m = disc.merge(out, on="gene", how="left")
    by_track = {}
    for tr, grp in m.groupby("track"):
        want = track_map[tr]
        ok = grp["zhang_top"].map(lambda z: z == want if isinstance(want, str) else z in want)
        by_track[tr] = {"genes": int(len(grp)), "measured": int(grp["zhang_top"].notna().sum()),
                        "same_track_in_zhang": int(ok.sum())}
    targets = load()["cell_context"]["targets"]
    tt = out[out["gene"].isin(targets)]
    both = out.dropna(subset=["zhang_tvp_log2fc", "li_tvp_log2fc"])
    summary = {"zhang": {"tumour_cells": 20509, "tumour_donors": 7, "normal_donors": 5,
                         "paired_donors": int(ztvp["zhang_tvp_donors"].max())},
               "genes_compared": int(len(out)),
               "dominant_group_agreement_all": round(float(out["agree_top"].mean()), 3),
               "rcc_share_spearman": round(float(stats.spearmanr(out["li_rcc_share"], out["zhang_rcc_share"])[0]), 3),
               "tvp_log2fc_spearman": round(float(stats.spearmanr(both["li_tvp_log2fc"], both["zhang_tvp_log2fc"])[0]), 3),
               "discovery_by_track": by_track,
               "targets": tt.set_index("gene")[["li_top", "zhang_top", "li_rcc_share", "zhang_rcc_share",
                                                "zhang_tvp_log2fc", "li_tvp_log2fc"]].round(3).to_dict("index")}
    (PROCESSED / "sc_replication_summary.json").write_text(json.dumps(summary, indent=1, default=float))
    return summary


if __name__ == "__main__":
    print(json.dumps(run(), indent=1, default=float))
