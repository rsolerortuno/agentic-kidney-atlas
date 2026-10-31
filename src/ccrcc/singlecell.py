"""Single-cell context (Li et al. 2022 scRNA-seq): which cell types express each gene in ccRCC.

Deterministic, no LLM. Raw UMI counts are streamed from the h5ad in row chunks (the matrix is never loaded
whole) and summed into donor x tissue x cell-group pseudobulks. Outputs, for all genes:

- sc_celltype_expression.parquet: per gene, tissue and cell group: % expressing cells, mean log-normalized
  expression, donor-averaged pseudobulk CPM, and the group's share of all transcripts of that gene in the tissue
  (answers "is the bulk tumor signal tumor-intrinsic or contributed by infiltrating cells?").
- sc_tumor_vs_pt.parquet: per donor, pseudobulk of ccRCC tumor cells (tumour region) vs normal proximal tubule
  cells (normal kidney; the cell of origin of ccRCC); paired t-test on log2(CPM + 1) + Benjamini-Hochberg.
- sc_reference.parquet: summed UMI counts per cell group over all kidney tissues (spatial deconvolution reference).
"""
from __future__ import annotations

import json

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats
from statsmodels.stats.multitest import multipletests

from .config import analysis_targets, load
from .omics import load_hgnc
from .paths import PROCESSED, RAW

SC_FILE = RAW / "singlecell" / "li2022_scrna.h5ad"
CHUNK = 20000

# cell group -> (Cell Ontology id, readable label)
CELL_GROUPS = {
    "RCC": ("CL:0001063", "ccRCC tumor cell"),
    "PT": ("CL:0002306", "proximal tubule epithelial cell"),
    "Epi_other": ("CL:0000066", "other renal epithelial cell"),
    "Endothelial": ("CL:0000115", "endothelial cell"),
    "Fibroblast": ("CL:0000057", "fibroblast / pericyte"),
    "Macrophage": ("CL:0000235", "macrophage (TAM)"),
    "Monocyte": ("CL:0000576", "monocyte"),
    "cDC": ("CL:0000990", "conventional dendritic cell"),
    "pDC": ("CL:0000784", "plasmacytoid dendritic cell"),
    "Mast": ("CL:0000097", "mast cell"),
    "CD8_T": ("CL:0000625", "CD8 T cell"),
    "CD4_T": ("CL:0000624", "CD4 T cell"),
    "Treg": ("CL:0000815", "regulatory T cell"),
    "T_other": ("CL:0000084", "other T cell (MAIT, gamma-delta)"),
    "NK": ("CL:0000623", "natural killer cell"),
    "B": ("CL:0000236", "B cell"),
    "Plasma": ("CL:0000786", "plasma cell"),
}
_BROAD = {"Epi_PT": "PT", "Epi_non-PT": "Epi_other", "EC": "Endothelial", "Fibro": "Fibroblast", "Mast": "Mast",
          "pDC": "pDC", "NK": "NK", "B-cell": "B", "Plasma": "Plasma"}


def cell_group(broad: str, annotation: str, cell_type: str) -> str | None:
    """Map the authors' broad type / fine annotation to one analysis group (None = not used)."""
    if broad == "RCC":
        # the authors' "MHCII" RCC sub-cluster was curated by CELLxGENE as MHC-II+ monocytes: excluded so that
        # MHC-II genes (e.g. CD74) are not attributed to tumor cells by possible myeloid doublets
        return None if annotation == "MHCII" else "RCC"
    if broad == "T-cell":
        if annotation.startswith("CD8"):
            return "CD8_T"
        if annotation == "CD4+Treg":
            return "Treg"
        return "CD4_T" if annotation.startswith("CD4") else "T_other"
    if broad == "Myeloid":
        if "macrophage" in cell_type:
            return "Macrophage"
        if "monocyte" in cell_type:
            return "Monocyte"
        return "cDC" if "dendritic" in cell_type else None
    return _BROAD.get(broad)


def cell_metadata() -> pd.DataFrame:
    cfg = load()["cell_context"]["scrna"]
    obs = ad.read_h5ad(SC_FILE, backed="r").obs
    o = pd.DataFrame({"row": np.arange(len(obs)), "donor": obs["donor_id"].astype(str).values,
                      "tissue": obs["summaryDescription"].astype(str).map(cfg["tissues"]).values,
                      "disease": obs["disease"].astype(str).values})
    o["group"] = [cell_group(b, a, c) for b, a, c in zip(obs["broad_type"].astype(str), obs["annotation"].astype(str),
                                                         obs["cell_type"].astype(str))]
    keep = o["tissue"].notna() & o["group"].notna() & ~o["disease"].isin(cfg["exclude_disease"])
    return o[keep].reset_index(drop=True)


def gene_table(h: h5py.File) -> pd.DataFrame:
    ens = h["raw/var/_index"].asstr()[:]
    hg = load_hgnc().dropna(subset=["ensembl_gene_id"]).drop_duplicates("ensembl_gene_id")
    sym = dict(zip(hg["ensembl_gene_id"], hg["symbol"]))
    return pd.DataFrame({"ensembl": ens, "hgnc_symbol": [sym.get(e) for e in ens]})


def pseudobulk(meta: pd.DataFrame) -> dict:
    """Stream raw counts; return per (donor, tissue, group) sums of counts, expressing cells and log-norm values."""
    keys = meta[["donor", "tissue", "group"]].drop_duplicates().sort_values(["donor", "tissue", "group"])
    keys = keys.reset_index(drop=True)
    kidx = {tuple(r): i for i, r in enumerate(keys.itertuples(index=False))}
    meta = meta.assign(key=[kidx[t] for t in zip(meta["donor"], meta["tissue"], meta["group"])])
    with h5py.File(SC_FILE, "r") as h:
        genes = gene_table(h)
        g = h["raw/X"]
        indptr = g["indptr"][:]
        n_genes = len(genes)
        counts = np.zeros((len(keys), n_genes))
        expressing = np.zeros((len(keys), n_genes))
        lognorm = np.zeros((len(keys), n_genes))
        rows = meta["row"].to_numpy()
        for start in range(0, len(rows), CHUNK):
            sel = rows[start:start + CHUNK]
            lo, hi = indptr[sel[0]], indptr[sel[-1] + 1]
            block = sp.csr_matrix((g["data"][lo:hi], g["indices"][lo:hi], indptr[sel[0]:sel[-1] + 2] - lo),
                                  shape=(sel[-1] - sel[0] + 1, n_genes))[sel - sel[0]]
            onehot = sp.csr_matrix((np.ones(len(sel)), (meta["key"].to_numpy()[start:start + CHUNK],
                                                         np.arange(len(sel)))), shape=(len(keys), len(sel)))
            counts += (onehot @ block).toarray()
            expressing += (onehot @ (block > 0).astype(np.float64)).toarray()
            lib = np.asarray(block.sum(axis=1)).ravel()
            norm = sp.diags(1e4 / np.maximum(lib, 1)) @ block
            norm.data = np.log1p(norm.data)
            lognorm += (onehot @ norm).toarray()
    keys["n_cells"] = meta.groupby("key").size().reindex(range(len(keys))).to_numpy()
    return {"keys": keys, "genes": genes, "counts": counts, "expressing": expressing, "lognorm": lognorm}


def celltype_table(pb: dict) -> pd.DataFrame:
    cfg = load()["cell_context"]["scrna"]
    keys, genes = pb["keys"], pb["genes"]
    cpm = pb["counts"] / np.maximum(pb["counts"].sum(axis=1, keepdims=True), 1) * 1e6
    out = []
    for tissue in sorted(keys["tissue"].unique()):
        kt = keys[keys["tissue"] == tissue]
        tissue_total = pb["counts"][kt.index].sum(axis=0)
        for group in CELL_GROUPS:
            kg = kt[kt["group"] == group]
            if kg.empty:
                continue
            n = kg["n_cells"].sum()
            usable = kg[kg["n_cells"] >= cfg["min_cells_pseudobulk"]]
            out.append(pd.DataFrame({
                "ensembl": genes["ensembl"], "hgnc_symbol": genes["hgnc_symbol"], "tissue": tissue, "cell_group": group,
                "cell_ontology": CELL_GROUPS[group][0], "n_cells": int(n), "n_donors": int(len(usable)),
                "pct_expressing": 100 * pb["expressing"][kg.index].sum(axis=0) / n,
                "mean_lognorm": pb["lognorm"][kg.index].sum(axis=0) / n,
                "pseudobulk_cpm": cpm[usable.index].mean(axis=0) if len(usable) else np.nan,
                "share_of_tissue_transcripts": pb["counts"][kg.index].sum(axis=0) / np.maximum(tissue_total, 1),
            }))
    df = pd.concat(out, ignore_index=True)
    df = df[df["hgnc_symbol"].notna()]
    for c in ("pct_expressing", "mean_lognorm", "pseudobulk_cpm", "share_of_tissue_transcripts"):
        df[c] = df[c].astype("float32")
    return df


def tau(values: np.ndarray) -> float:
    """Tissue/cell-type specificity index (Yanai 2005): 0 = ubiquitous, 1 = one group only."""
    v = np.asarray(values, dtype=float)
    if v.max() <= 0:
        return float("nan")
    return float((1 - v / v.max()).sum() / (len(v) - 1))


def tumor_vs_pt(pb: dict) -> pd.DataFrame:
    """Paired pseudobulk: ccRCC tumor cells (tumour region) vs normal proximal tubule (normal kidney)."""
    cfg = load()["cell_context"]["scrna"]
    keys, genes = pb["keys"], pb["genes"]
    cpm = pb["counts"] / np.maximum(pb["counts"].sum(axis=1, keepdims=True), 1) * 1e6
    ok = keys["n_cells"] >= cfg["min_cells_pseudobulk"]
    t = pd.Series(keys.index, index=keys["donor"])[(ok & (keys["tissue"] == "tumor") & (keys["group"] == "RCC")).values]
    n = pd.Series(keys.index, index=keys["donor"])[(ok & (keys["tissue"] == "normal") & (keys["group"] == "PT")).values]
    donors = sorted(set(t.index) & set(n.index))
    lt = np.log2(cpm[t.loc[donors].to_numpy()] + 1)
    ln = np.log2(cpm[n.loc[donors].to_numpy()] + 1)
    diff = lt - ln
    with np.errstate(invalid="ignore", divide="ignore"):
        res = stats.ttest_rel(lt, ln, axis=0)
    p = np.asarray(res.pvalue, dtype=float)
    testable = np.isfinite(p) & ((lt > 0) | (ln > 0)).any(axis=0)
    padj = np.full(len(p), np.nan)
    padj[testable] = multipletests(p[testable], method="fdr_bh")[1]
    df = pd.DataFrame({"ensembl": genes["ensembl"], "hgnc_symbol": genes["hgnc_symbol"],
                       "log2fc": diff.mean(axis=0), "pvalue": p, "padj": padj, "n_donors": len(donors),
                       "mean_log2cpm_tumor": lt.mean(axis=0), "mean_log2cpm_pt": ln.mean(axis=0),
                       "donors": ",".join(donors)})
    return df[df["hgnc_symbol"].notna() & testable]


def reference_profiles(pb: dict) -> pd.DataFrame:
    """Summed UMI counts per cell group over all kidney tissues (groups can be merged exactly downstream)."""
    keys, genes = pb["keys"], pb["genes"]
    rows = {}
    for group in CELL_GROUPS:
        idx = keys.index[keys["group"] == group]
        if len(idx):
            rows[group] = pb["counts"][idx].sum(axis=0)
    ref = pd.DataFrame(rows, index=genes["ensembl"])
    ref.insert(0, "hgnc_symbol", genes["hgnc_symbol"].to_numpy())
    return ref


def target_summary(ct: pd.DataFrame, tvp: pd.DataFrame, targets: list[str]) -> list[dict]:
    cfg = load()["cell_context"]["scrna"]
    out = []
    tum = ct[ct["tissue"] == "tumor"]
    for g in targets:
        s = tum[tum["hgnc_symbol"] == g].set_index("cell_group")
        if s.empty:
            out.append({"gene": g, "measured": False})
            continue
        order = s.sort_values("mean_lognorm", ascending=False)
        d = tvp[tvp["hgnc_symbol"] == g]
        out.append({
            "gene": g, "measured": True,
            "top_groups": [{"cell_group": k, "label": CELL_GROUPS[k][1], "mean_lognorm": round(float(r.mean_lognorm), 3),
                            "pct_expressing": round(float(r.pct_expressing), 1),
                            "share_of_tumor_transcripts": round(float(r.share_of_tissue_transcripts), 3)}
                           for k, r in order.head(3).iterrows()],
            "expressing_groups": [k for k, r in order.iterrows() if r.pct_expressing >= cfg["expressing_min_pct"]],
            "tumor_cell_share": round(float(s.loc["RCC", "share_of_tissue_transcripts"]), 3) if "RCC" in s.index else None,
            "tau": round(tau(s["mean_lognorm"].to_numpy()), 3),
            "tumor_vs_pt": None if d.empty else {"log2fc": round(float(d.iloc[0].log2fc), 3),
                                                 "padj": float(f"{d.iloc[0].padj:.3g}"), "n_donors": int(d.iloc[0].n_donors)},
        })
    return out


def write_summary() -> dict:
    """Per-target summary (light: reads the parquet tables). Re-run after the discovery shortlist changes."""
    ct = pd.read_parquet(PROCESSED / "sc_celltype_expression.parquet")
    tvp = pd.read_parquet(PROCESSED / "sc_tumor_vs_pt.parquet")
    meta = json.loads((PROCESSED / "sc_meta.json").read_text())
    summary = {**meta, "tumor_vs_pt_donors": int(tvp["n_donors"].iloc[0]) if len(tvp) else 0,
               "targets": target_summary(ct, tvp, analysis_targets())}
    (PROCESSED / "sc_summary.json").write_text(json.dumps(summary, indent=1))
    return summary


def run() -> dict:
    meta = cell_metadata()
    pb = pseudobulk(meta)
    ct = celltype_table(pb)
    ct.to_parquet(PROCESSED / "sc_celltype_expression.parquet", index=False)
    tvp = tumor_vs_pt(pb)
    tvp.to_parquet(PROCESSED / "sc_tumor_vs_pt.parquet", index=False)
    reference_profiles(pb).to_parquet(PROCESSED / "sc_reference.parquet")
    pb["keys"].to_parquet(PROCESSED / "sc_pseudobulk_keys.parquet", index=False)
    (PROCESSED / "sc_meta.json").write_text(json.dumps({
        "cells_used": int(len(meta)), "cells_by_group": meta["group"].value_counts().to_dict(),
        "cells_by_tissue": meta["tissue"].value_counts().to_dict(), "donors": int(meta["donor"].nunique())}))
    return write_summary()


if __name__ == "__main__":
    print(json.dumps(run(), indent=1)[:4000])
