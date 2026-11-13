"""Immune microenvironment as its own line of evidence.

Deterministic, no LLM. Tumour-vs-normal scores treat immune genes as noise from infiltration; here infiltration is the
signal. Inputs: CIBERSORT (absolute mode, 22 leukocyte types) and xCell (64 cell types plus immune / stroma scores)
computed by WashU for the 110 CPTAC discovery tumours (cptac Zenodo record; adjacent normal samples are left out), the scRNA-seq cell attribution, the Visium
neighbourhoods and the IMmotion150 response layer.

1. Tumour immune classes: xCell ImmuneScore tertiles (cold / intermediate / hot), with the CD8 T-cell share and the
   M2 / M1 macrophage ratio reported per class.
2. Gene-level infiltration coupling: Spearman correlation of CPTAC tumour RNA with CIBERSORT CD8 T cells, regulatory
   T cells, M2 macrophages and the xCell ImmuneScore (BH within cell type).
3. Immune-target profile, weights fixed before ranking, components in [0, 1]:
     immune_source  1.0  share of the gene's tumour transcripts in immune cells (scRNA-seq)
     coupling       1.0  max |rho| with CD8 T / Treg / M2 / ImmuneScore (q < 0.05), saturating at 0.6
     io_response    1.0  favourable PFS on atezolizumab + bevacizumab (q < 0.25): min(1, -log HR / 0.5)
     tumour_up      1.0  CPTAC RNA log2FC / 3 when significant (the target is enriched in tumours)
     tractable      0.5  predicted small-molecule or antibody tractability (Open Targets)
   Checkpoint genes with approved or late-stage drugs (PDCD1, CD274, CTLA4, LAG3, HAVCR2, TIGIT, IDO1) are positive
   controls, not inputs: no drug information enters the profile.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd

from .layers import EXTRA, _rowwise_spearman, bh
from .omics import LAYERS, read_matrix
from .paths import PROCESSED

WEIGHTS = {"immune_source": 1.0, "coupling": 1.0, "io_response": 1.0, "tumour_up": 1.0, "tractable": 0.5}
CONTROLS = ["PDCD1", "CD274", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "IDO1"]
CELLS = {"T cells CD8": "cd8", "T cells regulatory (Tregs)": "treg", "Macrophages M2": "m2", "ImmuneScore": "immune_score"}
IMMUNE_GROUPS = {"CD8_T", "CD4_T", "Treg", "T_other", "NK", "B", "Plasma", "Macrophage", "Monocyte", "cDC", "pDC", "Mast"}


def deconvolution() -> pd.DataFrame:
    # The WashU files hold the 110 tumours (-T) and 75 adjacent normal kidney samples (-A). Only tumours are used:
    # with the normals left in, the "cold" third of the immune classes was made of normal kidney (fixed 2026-11-13).
    cb = pd.read_csv(EXTRA / "cptac_ccrcc_cibersort_abs.txt.gz", sep="\t", index_col=0)
    cb = cb[cb.index.str.endswith("-T")]
    cb.index = cb.index.str.replace(r"-T$", "", regex=True)
    cb = cb.drop(columns=[c for c in ("P-value", "Correlation", "RMSE", "Absolute score (sig.score)") if c in cb.columns])
    xc = pd.read_csv(EXTRA / "cptac_ccrcc_xcell.txt.gz", sep="\t", index_col=0).T
    xc = xc[xc.index.str.endswith("-T")]
    xc.index = xc.index.str.replace(r"-T$", "", regex=True)
    cb.columns = ["cibersort:" + c for c in cb.columns]
    xc.columns = ["xcell:" + c for c in xc.columns]
    d = cb.join(xc, how="outer")
    d.index.name = "case"
    return d


def tumour_classes(d: pd.DataFrame) -> pd.DataFrame:
    s = d["xcell:ImmuneScore"]
    q = s.quantile([1 / 3, 2 / 3])
    cls = np.select([s <= q.iloc[0], s <= q.iloc[1]], ["cold", "intermediate"], "hot")
    out = pd.DataFrame({"case": d.index, "immune_class": cls, "immune_score": s.values,
                        "cd8": d["cibersort:T cells CD8"].values, "m2": d["cibersort:Macrophages M2"].values,
                        "m1": d["cibersort:Macrophages M1"].values, "treg": d["cibersort:T cells regulatory (Tregs)"].values})
    return out


def coupling(d: pd.DataFrame) -> pd.DataFrame:
    rna = read_matrix(LAYERS["rna"][0])
    parts = []
    for col, key in CELLS.items():
        src = "xcell:" + col if col == "ImmuneScore" else "cibersort:" + col
        v = d[src].reindex(rna.columns)
        P = pd.DataFrame(np.tile(v.to_numpy(float), (len(rna), 1)), index=rna.index, columns=rna.columns)
        r = _rowwise_spearman(rna, P, min_n=40).set_index("gene")
        r["padj"] = bh(r["pvalue"])
        parts.append(r[["rho", "padj"]].add_prefix(f"{key}_"))
    out = pd.concat(parts, axis=1)
    out.index.name = "dataset_symbol"
    return out.reset_index()


def profile(cp: pd.DataFrame) -> pd.DataFrame:
    from .scores import tractability_table
    u = pd.read_parquet(PROCESSED / "gene_universe.parquet").set_index("gene")
    ct = pd.read_parquet(PROCESSED / "sc_celltype_expression.parquet")
    ct = ct[ct["tissue"] == "tumor"]
    imm = ct[ct["cell_group"].isin(IMMUNE_GROUPS)].groupby("hgnc_symbol")["share_of_tissue_transcripts"].sum()
    cp = cp[cp["dataset_symbol"].isin(u.index)].set_index("dataset_symbol")
    t = pd.DataFrame(index=u.index)
    t["immune_share"] = imm.reindex(t.index)
    t["c_immune_source"] = t["immune_share"].clip(0, 1).fillna(0)
    rhos = []
    for key in CELLS.values():
        r = cp[f"{key}_rho"].where(cp[f"{key}_padj"] < 0.05).abs()
        t[f"{key}_rho"] = cp[f"{key}_rho"].reindex(t.index)
        rhos.append(r.reindex(t.index))
    t["c_coupling"] = (pd.concat(rhos, axis=1).max(axis=1) / 0.6).clip(0, 1).fillna(0)
    tr = pd.read_parquet(PROCESSED / "layer_treatment.parquet") if (PROCESSED / "layer_treatment.parquet").exists() else None
    if tr is not None:
        ab = tr[tr["arm"] == "atezo_bev"].set_index("gene")
        t["io_log_hr"], t["io_padj"] = ab["cox_log_hr"].reindex(t.index), ab["cox_padj"].reindex(t.index)
        t["c_io_response"] = ((-t["io_log_hr"] / 0.5).clip(0, 1) * (t["io_padj"] < 0.25)).fillna(0)
    else:
        t["c_io_response"] = 0.0
    t["rna_log2fc"] = u["rna_log2fc"]
    t["c_tumour_up"] = ((u["rna_log2fc"] / 3).clip(0, 1) * (u["rna_padj"] < 0.05)).fillna(0)
    trt = tractability_table()
    t["c_tractable"] = t.index.map(lambda g: float(trt["predicted_sm_ab"].get(g, False)) if len(trt) else 0.0)
    t["immune_target_score"] = sum(w * t[f"c_{k}"] for k, w in WEIGHTS.items())
    t = t.sort_values("immune_target_score", ascending=False)
    t["rank"] = np.arange(1, len(t) + 1)
    t["control"] = t.index.isin(CONTROLS)
    t.index.name = "gene"
    return t.reset_index()


def niche(genes: list[str]) -> pd.DataFrame:
    """Immune neighbours of the gene's high spots (Visium), where the gene was analysed spatially."""
    p = PROCESSED / "spatial_neighbors.parquet"
    if not p.exists():
        return pd.DataFrame()
    n = pd.read_parquet(p)
    n = n[n["hgnc_symbol"].isin(genes)] if "hgnc_symbol" in n.columns else n[n["gene"].isin(genes)]
    return n


def run() -> dict:
    d = deconvolution()
    d.to_parquet(PROCESSED / "immune_deconvolution.parquet")
    cls = tumour_classes(d)
    cls.to_parquet(PROCESSED / "immune_classes.parquet", index=False)
    cp = coupling(d)
    cp.to_parquet(PROCESSED / "immune_coupling.parquet", index=False)
    prof = profile(cp)
    prof.to_parquet(PROCESSED / "immune_profile.parquet", index=False)
    ctrl = prof[prof["control"]].set_index("gene")
    summary = {
        "tumours": int(len(d)), "classes": cls["immune_class"].value_counts().to_dict(),
        "class_means": cls.groupby("immune_class")[["cd8", "m2", "m1", "treg", "immune_score"]].mean().round(4).to_dict("index"),
        "coupled_genes": {k: int((cp[f"{k}_padj"] < 0.05).sum()) for k in CELLS.values()},
        "weights": WEIGHTS,
        "controls": {g: {"rank": int(r["rank"]), "score": round(float(r["immune_target_score"]), 3)} for g, r in ctrl.iterrows()},
        "controls_in_top100": int((ctrl["rank"] <= 100).sum()), "universe": int(len(prof)),
        "top": prof.head(20)["gene"].tolist(),
    }
    (PROCESSED / "immune_summary.json").write_text(json.dumps(summary, indent=1, default=float))
    return summary


if __name__ == "__main__":
    print(json.dumps(run(), indent=1, default=float))
