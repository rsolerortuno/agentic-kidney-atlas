"""Deterministic per-gene omics evidence table (input to literature selection and convergence)."""
from __future__ import annotations

import numpy as np
import pandas as pd

from .config import load
from .paths import PROCESSED


def omics_gene_table() -> pd.DataFrame:
    cfg = load()["omics"]
    de = pd.read_parquet(PROCESSED / "de_results.parquet")
    de = de[de["hgnc_symbol"].notna()]
    # one row per (layer, HGNC symbol): keep the most significant if several dataset symbols map to it
    de = de.sort_values("padj").drop_duplicates(["layer", "hgnc_symbol"])
    wide = {}
    for layer in ("rna", "protein", "phospho"):
        d = de[de["layer"] == layer].set_index("hgnc_symbol")
        thr = cfg["min_abs_log2fc"][layer]
        wide[f"{layer}_log2fc"] = d["log2fc"]
        wide[f"{layer}_padj"] = d["padj"]
        wide[f"{layer}_sig"] = (d["padj"] < cfg["padj"]) & (d["log2fc"].abs() >= thr)
    t = pd.DataFrame(wide)
    for layer in ("rna", "protein", "phospho"):
        t[f"{layer}_sig"] = t[f"{layer}_sig"].fillna(False).astype(bool)
    t["concordant"] = (t["rna_sig"] & t["protein_sig"]
                       & (np.sign(t["rna_log2fc"]) == np.sign(t["protein_log2fc"])))
    cnv = pd.read_parquet(PROCESSED / "cnv_frequencies.parquet")
    cnv = cnv[cnv["hgnc_symbol"].notna()].drop_duplicates("hgnc_symbol").set_index("hgnc_symbol")
    t = t.join(cnv[["loss_freq", "gain_freq"]], how="outer")
    for col in ("rna_sig", "protein_sig", "phospho_sig", "concordant"):
        t[col] = t[col].fillna(False).astype(bool)
    t["cnv_freq"] = t[["loss_freq", "gain_freq"]].max(axis=1)
    t.index.name = "gene"
    return t.reset_index()


def top_omics_candidates(n: int) -> list[str]:
    """Genes significant in RNA and protein with concordant direction, ranked by |protein log2FC| then |RNA log2FC|."""
    t = omics_gene_table()
    c = t[t["concordant"]].copy()
    c["_p"] = c["protein_log2fc"].abs()
    c["_r"] = c["rna_log2fc"].abs()
    return c.sort_values(["_p", "_r"], ascending=False)["gene"].head(n).tolist()
