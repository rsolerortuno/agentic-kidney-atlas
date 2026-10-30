"""M2: deterministic omics analysis (tumor vs normal DE, CNV frequencies, RNA-protein correlation).

All statistics are computed here; no LLM is involved.
"""
from __future__ import annotations

import json
import re

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

from .paths import PROCESSED, RAW

LAYERS = {
    "rna": ("HS_CPTAC_CCRCC_RNAseq_fpkm_log2_Tumor.cct", "HS_CPTAC_CCRCC_RNAseq_fpkm_log2_Normal.cct"),
    "protein": ("HS_CPTAC_CCRCC_proteome_Tumor.cct", "HS_CPTAC_CCRCC_proteome_Normal.cct"),
    "phospho": ("HS_CPTAC_CCRCC_phosphoproteome_gene_Tumor.cct", "HS_CPTAC_CCRCC_phosphoproteome_gene_Normal.cct"),
}
UNITS = {"rna": "log2 FPKM", "protein": "TMT log2 ratio", "phospho": "TMT log2 ratio (gene-level)"}
MIN_PAIRS = 10
MIN_GROUP = 5
CNV_THRESHOLD = 0.2  # |log2 ratio| above which a gene is called gained / lost


def read_matrix(name: str) -> pd.DataFrame:
    df = pd.read_csv(RAW / name, sep="\t", index_col=0)
    df.index = df.index.astype(str)
    return df.apply(pd.to_numeric, errors="coerce")


# ---------------------------------------------------------------- ID mapping
def load_hgnc() -> pd.DataFrame:
    cols = ["hgnc_id", "symbol", "name", "locus_group", "location", "alias_symbol",
            "prev_symbol", "entrez_id", "ensembl_gene_id", "uniprot_ids"]
    hg = pd.read_csv(RAW / "hgnc_complete_set.txt", sep="\t", usecols=cols, dtype=str, low_memory=False)
    return hg[hg["hgnc_id"].notna()]


def build_symbol_map(symbols, hgnc: pd.DataFrame) -> pd.DataFrame:
    """Map dataset symbols to HGNC: approved symbol first, then unique previous/alias symbol."""
    approved = hgnc.set_index("symbol")
    secondary: dict[str, set[str]] = {}
    for col in ("prev_symbol", "alias_symbol"):
        for sym, aliases in zip(hgnc["symbol"], hgnc[col]):
            if isinstance(aliases, str):
                for a in aliases.split("|"):
                    secondary.setdefault(a.strip(), set()).add(sym)
    rows = []
    for s in symbols:
        if s in approved.index:
            target, how = s, "approved"
        elif len(secondary.get(s, ())) == 1:
            target, how = next(iter(secondary[s])), "previous_or_alias"
        else:
            target, how = None, "unmapped"
        rec = {"dataset_symbol": s, "symbol": target, "mapping": how}
        if target is not None:
            r = approved.loc[target]
            rec.update({k: r[k] for k in ("hgnc_id", "name", "location", "entrez_id", "ensembl_gene_id", "uniprot_ids")})
        rows.append(rec)
    return pd.DataFrame(rows)


def chromosome_arm(location) -> str | None:
    if not isinstance(location, str):
        return None
    m = re.match(r"^(\d+|X|Y)([pq])", location)
    return f"{m.group(1)}{m.group(2)}" if m else None


# ---------------------------------------------------------------- DE
def differential(tumor: pd.DataFrame, normal: pd.DataFrame) -> pd.DataFrame:
    """Tumor vs normal per gene. Paired t-test on matched cases when >= MIN_PAIRS pairs,
    otherwise Welch t-test when each group has >= MIN_GROUP values. Wilcoxon reported alongside."""
    genes = tumor.index.intersection(normal.index)
    paired_cases = tumor.columns.intersection(normal.columns)
    out = []
    for g in genes:
        t = tumor.loc[g]
        n = normal.loc[g]
        tp, np_ = t[paired_cases], n[paired_cases]
        ok = tp.notna() & np_.notna()
        rec = {"gene": g, "n_tumor": int(t.notna().sum()), "n_normal": int(n.notna().sum())}
        if ok.sum() >= MIN_PAIRS:
            d = (tp[ok] - np_[ok]).to_numpy()
            rec["method"] = "paired_t"
            rec["n_pairs"] = int(ok.sum())
            rec["log2fc"] = float(d.mean())
            if np.allclose(d, d[0]):
                rec["stat"], rec["pvalue"], rec["wilcoxon_p"] = np.nan, np.nan, np.nan
            else:
                res = stats.ttest_rel(tp[ok], np_[ok])
                rec["stat"], rec["pvalue"] = float(res.statistic), float(res.pvalue)
                rec["wilcoxon_p"] = float(stats.wilcoxon(d).pvalue)
        elif t.notna().sum() >= MIN_GROUP and n.notna().sum() >= MIN_GROUP:
            res = stats.ttest_ind(t.dropna(), n.dropna(), equal_var=False)
            rec["method"] = "welch_t"
            rec["n_pairs"] = int(ok.sum())
            rec["log2fc"] = float(t.mean() - n.mean())
            rec["stat"], rec["pvalue"] = float(res.statistic), float(res.pvalue)
            rec["wilcoxon_p"] = float(stats.mannwhitneyu(t.dropna(), n.dropna()).pvalue)
        else:
            continue
        rec["mean_tumor"] = float(t.mean())
        rec["mean_normal"] = float(n.mean())
        out.append(rec)
    df = pd.DataFrame(out)
    df = df[df["pvalue"].notna()].copy()
    df["padj"] = multipletests(df["pvalue"], method="fdr_bh")[1]
    df["direction"] = np.where(df["log2fc"] > 0, "up", "down")
    return df.sort_values("padj").reset_index(drop=True)


# ---------------------------------------------------------------- CNV
def cnv_frequencies(cnv: pd.DataFrame) -> pd.DataFrame:
    n = cnv.notna().sum(axis=1)
    return pd.DataFrame({
        "gene": cnv.index,
        "n_tumor": n.values,
        "loss_freq": ((cnv < -CNV_THRESHOLD).sum(axis=1) / n).values,
        "gain_freq": ((cnv > CNV_THRESHOLD).sum(axis=1) / n).values,
        "mean_log2ratio": cnv.mean(axis=1).values,
    })


# ---------------------------------------------------------------- RNA-protein
def rna_protein_correlation(rna_t: pd.DataFrame, prot_t: pd.DataFrame, min_n: int = 20) -> pd.DataFrame:
    genes = rna_t.index.intersection(prot_t.index)
    cases = rna_t.columns.intersection(prot_t.columns)
    rows = []
    for g in genes:
        r, p = rna_t.loc[g, cases], prot_t.loc[g, cases]
        ok = r.notna() & p.notna()
        if ok.sum() < min_n or r[ok].nunique() < 3 or p[ok].nunique() < 3:
            continue
        rho, pv = stats.spearmanr(r[ok], p[ok])
        rows.append({"gene": g, "n": int(ok.sum()), "spearman_rho": float(rho), "pvalue": float(pv)})
    df = pd.DataFrame(rows)
    df["padj"] = multipletests(df["pvalue"], method="fdr_bh")[1]
    return df


# ---------------------------------------------------------------- pipeline
def _annotate(df: pd.DataFrame, symmap: pd.DataFrame) -> pd.DataFrame:
    m = symmap.rename(columns={"dataset_symbol": "gene", "symbol": "hgnc_symbol"})
    return df.merge(m[["gene", "hgnc_symbol", "hgnc_id", "ensembl_gene_id", "uniprot_ids", "location", "mapping"]],
                    on="gene", how="left")


def run() -> dict:
    hgnc = load_hgnc()
    all_symbols: set[str] = set()
    summary: dict = {"layers": {}}

    de_tables = {}
    tumors = {}
    for layer, (tf, nf) in LAYERS.items():
        t, n = read_matrix(tf), read_matrix(nf)
        tumors[layer] = t
        de = differential(t, n)
        de.insert(0, "layer", layer)
        de["unit"] = UNITS[layer]
        de_tables[layer] = de
        all_symbols |= set(de["gene"])
        summary["layers"][layer] = {
            "genes_tested": int(len(de)),
            "significant_padj_0.05": int((de["padj"] < 0.05).sum()),
            "up": int(((de["padj"] < 0.05) & (de["log2fc"] > 0)).sum()),
            "down": int(((de["padj"] < 0.05) & (de["log2fc"] < 0)).sum()),
            "paired_cases": int(len(t.columns.intersection(n.columns))),
        }

    cnv = read_matrix("HS_CPTAC_CCRCC_CNV_gene_Tumor.cct")
    cnvf = cnv_frequencies(cnv)
    all_symbols |= set(cnvf["gene"])
    corr = rna_protein_correlation(tumors["rna"], tumors["protein"])

    symmap = build_symbol_map(sorted(all_symbols), hgnc)
    symmap.to_parquet(PROCESSED / "symbol_map.parquet", index=False)
    summary["id_mapping"] = symmap["mapping"].value_counts().to_dict()

    de_all = pd.concat([_annotate(d, symmap) for d in de_tables.values()], ignore_index=True)
    de_all.to_parquet(PROCESSED / "de_results.parquet", index=False)

    cnvf = _annotate(cnvf, symmap)
    cnvf["arm"] = cnvf["location"].map(chromosome_arm)
    cnvf.to_parquet(PROCESSED / "cnv_frequencies.parquet", index=False)
    arm = (cnvf.dropna(subset=["arm"]).groupby("arm")
           .agg(n_genes=("gene", "size"), median_loss_freq=("loss_freq", "median"),
                median_gain_freq=("gain_freq", "median")).reset_index())
    arm.to_parquet(PROCESSED / "cnv_arm_summary.parquet", index=False)
    summary["cnv_top_arm_losses"] = arm.sort_values("median_loss_freq", ascending=False).head(5).to_dict("records")
    summary["cnv_top_arm_gains"] = arm.sort_values("median_gain_freq", ascending=False).head(5).to_dict("records")

    corr = _annotate(corr, symmap)
    corr.to_parquet(PROCESSED / "rna_protein_correlation.parquet", index=False)
    summary["rna_protein_correlation"] = {"genes": int(len(corr)), "median_rho": float(corr["spearman_rho"].median())}

    (PROCESSED / "omics_summary.json").write_text(json.dumps(summary, indent=2, default=str))
    return summary


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, default=str))
