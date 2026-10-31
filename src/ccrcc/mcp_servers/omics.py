"""omics-mcp: read-only access to CPTAC ccRCC omics results (never raw matrices)."""
from __future__ import annotations

from functools import lru_cache

import pandas as pd
from mcp.server.fastmcp import FastMCP

from ..config import load
from ..graph import evidence_id
from ..paths import PROCESSED
from .common import canonical_gene, clean, existing

mcp = FastMCP("omics-mcp")
LAYERS = ("rna", "protein", "phospho")


@lru_cache(maxsize=1)
def _de() -> pd.DataFrame:
    de = pd.read_parquet(PROCESSED / "de_results.parquet").dropna(subset=["hgnc_symbol"])
    return de.sort_values("padj").drop_duplicates(["layer", "hgnc_symbol"])


def _row(r) -> dict:
    thr = load()["omics"]
    sig = bool(r["padj"] < thr["padj"] and abs(r["log2fc"]) >= thr["min_abs_log2fc"][r["layer"]])
    return clean({"gene": r["hgnc_symbol"], "dataset_symbol": r["gene"], "layer": r["layer"],
                  "log2fc": round(r["log2fc"], 3), "padj": float(f"{r['padj']:.3g}"), "method": r["method"],
                  "n_pairs": r["n_pairs"], "n_tumor": r["n_tumor"], "n_normal": r["n_normal"], "unit": r["unit"],
                  "significant": sig,
                  "evidence_id": existing(evidence_id("omics_de", r["layer"], r["hgnc_symbol"]))})


@mcp.tool()
def get_de_results(gene: str, layer: str | None = None) -> list[dict]:
    """Tumor-vs-normal differential expression for one gene in CPTAC ccRCC.

    layer: 'rna', 'protein' or 'phospho' (default: all layers). Returns log2 fold change (tumor - normal),
    BH-adjusted p-value, test used, sample sizes, whether it passes the configured thresholds, and the
    evidence_id to cite."""
    g = canonical_gene(gene)
    de = _de()
    sub = de[de["hgnc_symbol"] == g]
    if layer:
        sub = sub[sub["layer"] == layer]
    if sub.empty:
        return [{"gene": g, "message": "not measured in the requested layer(s)"}]
    return [_row(r) for _, r in sub.iterrows()]


@mcp.tool()
def top_genes(layer: str = "protein", n: int = 20, direction: str = "up") -> list[dict]:
    """Top differentially expressed genes in one layer ('rna', 'protein', 'phospho'), ranked by |log2FC|
    among genes passing padj and fold-change thresholds. direction: 'up', 'down' or 'both'. n <= 100."""
    if layer not in LAYERS:
        raise ValueError(f"layer must be one of {LAYERS}")
    thr = load()["omics"]
    d = _de()
    d = d[(d["layer"] == layer) & (d["padj"] < thr["padj"]) & (d["log2fc"].abs() >= thr["min_abs_log2fc"][layer])]
    if direction == "up":
        d = d[d["log2fc"] > 0]
    elif direction == "down":
        d = d[d["log2fc"] < 0]
    d = d.reindex(d["log2fc"].abs().sort_values(ascending=False).index).head(min(n, 100))
    return [_row(r) for _, r in d.iterrows()]


@mcp.tool()
def rna_protein_correlation(gene: str) -> dict:
    """Spearman correlation between RNA and protein abundance across CPTAC ccRCC tumors for one gene."""
    g = canonical_gene(gene)
    c = pd.read_parquet(PROCESSED / "rna_protein_correlation.parquet")
    sub = c[c["hgnc_symbol"] == g]
    if sub.empty:
        return {"gene": g, "message": "not available (gene missing from RNA or proteome, or < 20 paired tumors)"}
    r = sub.iloc[0]
    return clean({"gene": g, "spearman_rho": round(r["spearman_rho"], 3), "padj": float(f"{r['padj']:.3g}"),
                  "n_tumors": r["n"], "evidence_id": existing(evidence_id("omics_rna_protein", g))})


@mcp.tool()
def cnv_frequency(gene: str) -> dict:
    """Fraction of CPTAC ccRCC tumors with copy-number loss / gain (|log2 ratio| > 0.2) for one gene,
    plus its chromosome arm. evidence_id is present when loss or gain frequency >= 0.3."""
    g = canonical_gene(gene)
    c = pd.read_parquet(PROCESSED / "cnv_frequencies.parquet")
    sub = c[c["hgnc_symbol"] == g]
    if sub.empty:
        return {"gene": g, "message": "not in CNV data"}
    r = sub.iloc[0]
    return clean({"gene": g, "loss_freq": round(r["loss_freq"], 3), "gain_freq": round(r["gain_freq"], 3),
                  "arm": r["arm"], "n_tumors": r["n_tumor"], "evidence_id": existing(evidence_id("omics_cnv", g))})


if __name__ == "__main__":
    mcp.run()
