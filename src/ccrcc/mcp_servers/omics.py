"""omics-mcp: read-only access to CPTAC ccRCC omics results (never raw matrices)."""
from __future__ import annotations

from functools import lru_cache

import pandas as pd
from mcp.server.fastmcp import FastMCP

from ..config import load
from ..graph import evidence_id
from ..paths import PROCESSED
from .common import applied_record, canonical_gene, clean, existing

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


# ---------------------------------------------------------------- M13 layers
def _layer(name: str, col: str = "gene") -> pd.DataFrame:
    p = PROCESSED / name
    if not p.exists():
        return pd.DataFrame()
    return pd.read_parquet(p)


def _one(df: pd.DataFrame, g: str, col: str = "gene"):
    if df.empty:
        return None
    sub = df[df[col] == g]
    return None if sub.empty else sub.iloc[0]


def _num(x, nd=3):
    if x is None or pd.isna(x):
        return None
    return float(f"{x:.3g}") if abs(x) < 1e-3 and x != 0 else round(float(x), nd)


@mcp.tool()
def mutation_recurrence(gene: str) -> dict:
    """Somatic nonsynonymous mutation recurrence in ccRCC (CPTAC discovery patients + TCGA-KIRC, 511 tumours),
    tested against a background rate scaled by the gene's coding length (so long genes such as TTN are not called
    drivers). Returns mutated patients, frequency, observed/expected and BH q-value, with the evidence_id."""
    g = canonical_gene(gene)
    r = _one(_layer("layer_mutations.parquet"), g)
    if r is None:
        return {"gene": g, "message": "no nonsynonymous mutation in either cohort (or no coding length)"}
    return clean({"gene": g, "mutated_patients": int(r["mutated"]), "patients": int(r["patients"]),
                  "freq": _num(r["freq"]), "cptac_mutated": int(r["cptac_mutated"]), "tcga_mutated": int(r["tcga_mutated"]),
                  "obs_exp": _num(r["obs_exp"]), "padj": _num(r["padj"]), "recurrent": bool(r["padj"] < 0.05),
                  "evidence_id": existing(evidence_id("somatic_recurrence", g))})


@mcp.tool()
def cis_copy_number(gene: str) -> dict:
    """Does the gene's copy-number change move its own RNA / protein in CPTAC tumours? Spearman rho and BH q for
    copy number vs RNA and vs protein. A frequent loss with a positive cis effect (e.g. 3p genes) is a real dosage
    effect; a frequent loss without one is a passenger of the arm-level event."""
    g = canonical_gene(gene)
    r = _one(_layer("layer_cis_dosage.parquet"), g)
    if r is None:
        return {"gene": g, "message": "not in copy-number data"}
    return clean({"gene": g, "rna_rho": _num(r["rna_rho"]), "rna_padj": _num(r["rna_padj"]),
                  "protein_rho": _num(r["protein_rho"]), "protein_padj": _num(r["protein_padj"]),
                  "evidence_id": existing(evidence_id("cis_dosage", g))})


@mcp.tool()
def depmap_dependency(gene: str) -> dict:
    """DepMap Public 26Q1 CRISPR (Chronos) dependency: mean gene effect in clear cell RCC lines vs all other lines
    (more negative = more essential), fraction of ccRCC lines with effect < -0.5, selectivity (other - ccRCC),
    Welch BH q-value and whether the gene is a common essential (dependent in >= 90 % of all lines: toxicity risk)."""
    g = canonical_gene(gene)
    r = _one(_layer("layer_dependency.parquet"), g)
    if r is None:
        return {"gene": g, "message": "not screened in DepMap 26Q1"}
    return clean({"gene": g, "n_ccrcc_lines": int(r["n_ccrcc_lines"]), "mean_effect_ccrcc": _num(r["mean_effect_ccrcc"]),
                  "mean_effect_other": _num(r["mean_effect_other"]), "selectivity": _num(r["selectivity"]),
                  "frac_ccrcc_dependent": _num(r["frac_ccrcc_dependent"]), "padj": _num(r["padj"]),
                  "common_essential": bool(r["common_essential"]),
                  "evidence_id": existing(evidence_id("depmap_dependency", g))})


@mcp.tool()
def tcga_replication(gene: str) -> dict:
    """External replication in TCGA-KIRC: paired tumour vs adjacent normal RNA (72 pairs; log2FC, BH q) and
    univariate Cox regression of overall survival on standardized tumour expression (log HR per SD, BH q).
    Positive log HR = higher expression, worse survival."""
    g = canonical_gene(gene)
    r = _one(_layer("layer_tcga.parquet"), g)
    if r is None:
        return {"gene": g, "message": "not in TCGA-KIRC RNA-seq"}
    return clean({"gene": g, "log2fc": _num(r["log2fc"]), "padj": _num(r["padj"]), "n_pairs": _num(r["n_pairs"], 0),
                  "cox_log_hr_per_sd": _num(r["cox_log_hr_per_sd"]), "cox_padj": _num(r["cox_padj"]),
                  "n_survival": _num(r["n_survival"], 0),
                  "evidence_id_de": existing(evidence_id("tcga_de", g)),
                  "evidence_id_survival": existing(evidence_id("tcga_survival", g))})


@mcp.tool()
def purity_effect(gene: str) -> dict:
    """Tumour purity (ESTIMATE) check: Spearman of tumour expression with purity (negative = comes with stroma /
    immune cells) and the paired log2FC extrapolated to a pure tumour (purity = 1) with its standard error, for RNA
    and protein. A large bulk log2FC that shrinks towards 0 at purity 1 is infiltration, not tumour-cell expression."""
    g = canonical_gene(gene)
    df = _layer("layer_purity.parquet")
    r = _one(df, g, col="hgnc_symbol")
    if r is None:
        return {"gene": g, "message": "not in CPTAC RNA / protein"}
    keys = ["rna_purity_rho", "rna_purity_padj", "rna_lfc_at_purity1", "rna_lfc_at_purity1_se",
            "protein_purity_rho", "protein_purity_padj", "protein_lfc_at_purity1", "protein_lfc_at_purity1_se"]
    return clean({"gene": g, **{k: _num(r[k]) for k in keys}, "evidence_id": existing(evidence_id("purity", g))})


@mcp.tool()
def normal_tissue_protein(gene: str) -> dict:
    """Protein in vital normal organs (GTEx proteome, Jiang et al. 2020; 16 vital tissues, no kidney): how many vital
    tissues detect it, the highest median log2 abundance relative to the run reference and where, and the highest
    tissue-specificity (TS) score. TS >= 4 = enriched in that organ, an on-target, off-tumour safety concern. Absence
    from this proteome is not proof of absence (detection limit). Cite evidence_id."""
    g = canonical_gene(gene)
    r = _one(_layer("normal_proteome.parquet"), g)
    if r is None:
        return {"gene": g, "message": "not detected in the GTEx normal-tissue proteome"}
    return clean({"gene": g, "vital_tissues_detected": _num(r["vital_protein_detected"], 0),
                  "vital_protein_max_rel": _num(r["vital_protein_max_rel"]), "vital_protein_max_tissue": r["vital_protein_max_tissue"],
                  "vital_protein_ts_max": _num(r["vital_protein_ts_max"]), "vital_protein_ts_tissue": r["vital_protein_ts_tissue"],
                  "vital_protein_enriched": bool(r["vital_protein_enriched"]),
                  "evidence_id": existing(evidence_id("normal_proteome", g))})


@mcp.tool()
def drug_response(gene: str) -> dict:
    """Do compounds against this gene's product kill ccRCC cell lines more than other lines? PRISM Repurposing 20Q2
    primary screen (2.5 uM, viability log2 fold change; 10-14 ccRCC lines vs about 550 others). Returns the most selective
    compound annotated to the gene (preferring compounds with <= 3 annotated targets), its mechanism, selectivity
    (mean LFC other - mean LFC ccRCC; positive = ccRCC lines more sensitive), BH q and how many compounds target the gene.
    Cell lines in 2D culture miss HIF-2 biology (belzutifan-type drugs are not selective here). Cite evidence_id."""
    g = canonical_gene(gene)
    r = _one(_layer("layer_prism.parquet"), g)
    if r is None:
        return {"gene": g, "message": "no PRISM compound is annotated to this gene"}
    return clean({"gene": g, "compound": r["compound"], "moa": r["moa"], "phase": r["phase"],
                  "compound_targets": _num(r["compound_targets"], 0), "n_compounds": _num(r["n_compounds"], 0),
                  "selectivity": _num(r["selectivity"]), "mean_lfc_ccrcc": _num(r["mean_lfc_ccrcc"]),
                  "mean_lfc_other": _num(r["mean_lfc_other"]), "padj": _num(r["padj"]),
                  "n_ccrcc_lines": _num(r["n_ccrcc_lines"], 0), "evidence_id": existing(evidence_id("prism_response", g))})


@mcp.tool()
def phosphosites(gene: str, n: int = 10) -> list[dict]:
    """Site-level tumour vs normal phosphorylation for one gene (CPTAC BCM pipeline, paired t-test, BH), with the
    site change adjusted for the protein's own change. evidence_id exists for sites with BH < 0.05 and |log2FC| >= 1."""
    g = canonical_gene(gene)
    df = _layer("layer_phosphosite.parquet")
    sub = df[df["gene"] == g].sort_values("padj").head(min(n, 30)) if not df.empty else df
    if sub.empty:
        return [{"gene": g, "message": "no quantified phosphosite"}]
    return [clean({"gene": g, "site": r["site"], "log2fc": _num(r["log2fc"]),
                   "protein_adjusted_log2fc": _num(r.get("protein_adjusted_log2fc")), "padj": _num(r["padj"]),
                   "n_pairs": int(r["n_pairs"]), "evidence_id": existing(evidence_id("phosphosite_de", r["site_id"]))})
            for _, r in sub.iterrows()]


@mcp.tool()
def kinase_activity(gene: str | None = None, n: int = 15) -> list[dict]:
    """Inferred kinase activity in tumours vs normal (KSEA on protein-adjusted phosphosite changes, OmniPath
    kinase-substrate relations). With gene: that kinase only; without: the top n by z (activated first)."""
    df = _layer("layer_kinase_activity.parquet")
    if df.empty:
        return [{"message": "kinase activity not computed"}]
    if gene:
        df = df[df["gene"] == canonical_gene(gene)]
        if df.empty:
            return [{"gene": canonical_gene(gene), "message": "fewer than 5 quantified substrate sites"}]
    return [clean({"gene": r["gene"], "ksea_z": _num(r["z"]), "padj": _num(r["padj"]),
                   "n_substrate_sites": int(r["n_substrate_sites"]),
                   "evidence_id": existing(evidence_id("kinase_activity", r["gene"]))})
            for _, r in df.head(min(n, 50)).iterrows()]


@mcp.tool()
def confirmatory_replication(gene: str) -> dict:
    """Does the protein (and phosphosite) change replicate in the confirmatory CPTAC cohort? That cohort has different
    patients and a different technology (label-free DIA mass spectrometry, Proteomic Data Commons PDC000411), so
    absolute fold changes differ from the discovery TMT data: judge direction and significance. Returns the
    confirmatory and discovery protein log2 fold change, adjusted p-value, pairs and whether it replicates, plus the
    number of phosphosites of the gene that were measured in both cohorts and replicate."""
    return applied_record(gene, "confirmatory_replication", "gene not tested in the confirmatory cohort: absent or too sparse there (fewer than 10 matched tumour-normal pairs)")


if __name__ == "__main__":
    mcp.run()
