"""Replication of protein and phosphosite changes in the confirmatory CPTAC ccRCC cohort.

Deterministic, no LLM. The discovery cohort (110 tumours, TMT proteomics) is where every protein and phosphosite result
of this project comes from. The confirmatory cohort is a separate set of patients measured by label-free DIA mass
spectrometry at Johns Hopkins (Li et al. 2023, Cancer Cell). Its processed matrices come from the Proteomic Data Commons
(study PDC000411, supplementary archive for PMID 36563681), downloaded to data/raw/confirmatory/.

Only patients that are not in the discovery cohort are used, so the test is independent. Intensities are log2
transformed and tumour is compared with the matched normal of the same patient (paired t-test, at least 10 pairs, BH
within each table), exactly as in the discovery analysis. The intratumour-heterogeneity replicates (CPT... _ITH) and
duplicate columns are not used. For every gene or site measured in both cohorts the table keeps the discovery log2FC
next to the confirmatory one and says whether the call replicates: same direction and q < 0.05 in the confirmatory
cohort. Phosphosites are matched on gene and residue (e.g. PABPN1_S19); multi-site peptides are skipped.

The different technology (TMT in discovery, DIA here) means absolute fold changes differ (TMT compresses them), so the
comparison reports direction, significance and rank correlation, not equality of effect sizes.
"""
from __future__ import annotations

import json
import re

import numpy as np
import pandas as pd
from scipy import stats

from .layers import MIN_PAIRS, bh, paired_de
from .omics import LAYERS, load_hgnc, read_matrix
from .paths import PROCESSED, RAW

CONF = RAW / "confirmatory"
PROT = CONF / "CPTAC_ccRCC_global_proteome_DIA.tsv"
PHOS = CONF / "CPTAC_ccRCC_phosphoproteome_DIA_site.tsv"
CASE = re.compile(r"^(C3[LN]-\d+)_([TN])$")


def _split(df: pd.DataFrame, exclude_cases: set[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    tum, nor = {}, {}
    for c in df.columns:
        m = CASE.match(c)
        if not m or m.group(1) in exclude_cases:
            continue
        (tum if m.group(2) == "T" else nor)[m.group(1)] = c
    t = df[list(tum.values())].copy()
    t.columns = list(tum)
    n = df[list(nor.values())].copy()
    n.columns = list(nor)
    return np.log2(t.where(t > 0)), np.log2(n.where(n > 0))


def _symbols() -> set[str]:
    return set(load_hgnc()["symbol"])


def proteome() -> pd.DataFrame:
    disc_cases = set(read_matrix(LAYERS["protein"][0]).columns)
    raw = pd.read_csv(PROT, sep="\t", index_col=0)
    raw = raw.groupby(level=0).mean(numeric_only=True)
    t, n = _split(raw, disc_cases)
    de = paired_de(t, n).rename(columns={"feature": "gene", "log2fc": "conf_log2fc", "padj": "conf_padj",
                                         "n_pairs": "conf_n_pairs", "pvalue": "conf_pvalue", "t": "conf_t"})
    d = pd.read_parquet(PROCESSED / "de_results.parquet")
    d = d[d["layer"] == "protein"].dropna(subset=["hgnc_symbol"]).drop_duplicates("hgnc_symbol").set_index("hgnc_symbol")
    out = de.set_index("gene").join(d[["log2fc", "padj"]].rename(columns={"log2fc": "disc_log2fc", "padj": "disc_padj"}),
                                    how="inner").reset_index()
    out["disc_sig"] = out["disc_padj"] < 0.05
    out["replicates"] = out["disc_sig"] & (out["conf_padj"] < 0.05) & (np.sign(out["conf_log2fc"]) == np.sign(out["disc_log2fc"]))
    out["n_tumour_conf"] = t.shape[1]
    out.to_parquet(PROCESSED / "confirmatory_protein.parquet", index=False)
    return out


def phospho() -> pd.DataFrame:
    disc_cases = set(read_matrix(LAYERS["protein"][0]).columns)
    raw = pd.read_csv(PHOS, sep="\t", index_col=0)
    sym = _symbols()
    keep, ids = [], []
    for name in raw.index:
        parts = str(name).split("_")           # NP_001347480_PABPN1_S19  (accession has an underscore)
        if len(parts) < 4 or ";" in parts[-1]:
            continue
        gene, site = parts[-2], parts[-1]
        if gene in sym and re.match(r"^[STY]\d+$", site):
            keep.append(name)
            ids.append(f"{gene}_{site}")
    raw = raw.loc[keep]
    raw.index = ids
    raw = raw.groupby(level=0).mean(numeric_only=True)
    t, n = _split(raw, disc_cases)
    de = paired_de(t, n).rename(columns={"feature": "site_id", "log2fc": "conf_log2fc", "padj": "conf_padj",
                                         "n_pairs": "conf_n_pairs", "pvalue": "conf_pvalue", "t": "conf_t"})
    d = pd.read_parquet(PROCESSED / "layer_phosphosite.parquet").sort_values("padj").drop_duplicates("site_id")
    d = d.set_index("site_id")[["log2fc", "padj", "gene"]].rename(columns={"log2fc": "disc_log2fc", "padj": "disc_padj"})
    out = de.set_index("site_id").join(d, how="inner").reset_index()
    out["disc_sig"] = out["disc_padj"] < 0.05
    out["replicates"] = out["disc_sig"] & (out["conf_padj"] < 0.05) & (np.sign(out["conf_log2fc"]) == np.sign(out["disc_log2fc"]))
    out.to_parquet(PROCESSED / "confirmatory_phospho.parquet", index=False)
    return out


def _stats(d: pd.DataFrame, key: str) -> dict:
    sig = d[d["disc_sig"]]
    rho = float(stats.spearmanr(d["disc_log2fc"], d["conf_log2fc"])[0])
    big = sig[sig["disc_log2fc"].abs() >= 0.5]
    return {"measured_in_both": int(len(d)), "discovery_significant": int(len(sig)),
            "replicated": int(sig["replicates"].sum()),
            "replication_rate": round(float(sig["replicates"].mean()), 3),
            "replication_rate_large_effect": round(float(big["replicates"].mean()), 3) if len(big) else None,
            "n_large_effect": int(len(big)),
            "direction_agreement_significant": round(float((np.sign(sig["conf_log2fc"]) == np.sign(sig["disc_log2fc"])).mean()), 3),
            "spearman_log2fc": round(rho, 3), "pairs_median": int(d["conf_n_pairs"].median())}


def run() -> dict:
    if not PROT.exists() or not PHOS.exists():
        raise FileNotFoundError("confirmatory matrices missing: run `ccrcc download-extra`")
    p, s = proteome(), phospho()
    cand = ["CA9", "NDUFA4L2", "ANGPTL4", "CD70", "NNMT", "ENPP3", "SCARB1", "SLC16A3", "HK2", "PSAT1", "PCK1", "SLC6A19",
            "ALDH6A1", "EPAS1", "MET", "EGFR", "KDR", "ALDOB", "SLC22A8", "PIGR", "UMOD"]
    pc = p[p["gene"].isin(cand)].set_index("gene")
    summary = {"cohort": {"proteome_tumours": int(p["n_tumour_conf"].iloc[0]), "median_pairs_protein": int(p["conf_n_pairs"].median())},
               "protein": _stats(p, "gene"), "phosphosite": _stats(s, "site_id"),
               "candidates": {g: {"conf_log2fc": round(float(r["conf_log2fc"]), 3), "conf_padj": float(f"{r['conf_padj']:.3g}"),
                                  "disc_log2fc": round(float(r["disc_log2fc"]), 3), "replicates": bool(r["replicates"])}
                              for g, r in pc.iterrows()}}
    (PROCESSED / "confirmatory_summary.json").write_text(json.dumps(summary, indent=1, default=float))
    return summary


if __name__ == "__main__":
    print(json.dumps(run(), indent=1, default=float))
