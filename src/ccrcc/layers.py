"""M13: extra evidence layers (deterministic, no LLM).

Each function reads data/raw/extra/ (see extra_download.py) and writes one table to data/processed/:

- mutations      layer_mutations.parquet      somatic recurrence in CPTAC (discovery patients) + TCGA-KIRC, corrected
                                              for coding length (binomial test against a per-cohort background rate)
- cis_dosage     layer_cis_dosage.parquet     CPTAC copy number -> RNA / protein Spearman correlation in tumors
                                              (does the copy-number change actually move the gene product?)
- dependency     layer_dependency.parquet     DepMap 26Q1 Chronos gene effect: ccRCC lines vs all other lines
- replication    layer_tcga.parquet           TCGA-KIRC paired tumor vs adjacent normal RNA + Cox overall survival
- phosphosites   layer_phosphosite.parquet    CPTAC site-level paired tumor vs normal (BCM pipeline)
                 layer_kinase_activity.parquet  KSEA-style kinase activity from OmniPath kinase -> site relations
- purity         layer_purity.parquet         ESTIMATE purity: correlation with tumor expression and the paired
                                              log2FC extrapolated to a pure tumor (purity = 1)

All p-values are Benjamini-Hochberg adjusted within each table.
"""
from __future__ import annotations

import gzip
import json
import re

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

from .config import load as load_cfg
from .omics import LAYERS, load_hgnc, read_matrix, to_hgnc
from .paths import PROCESSED, RAW

EXTRA = RAW / "extra"
NONSYN = {"Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site",
          "In_Frame_Del", "In_Frame_Ins", "Nonstop_Mutation", "Translation_Start_Site"}
HYPERMUTATOR = 500     # tumours with more nonsynonymous mutations are excluded from the background rate
MIN_PAIRS = 10


def bh(p) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    out = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    if ok.any():
        out[ok] = multipletests(p[ok], method="fdr_bh")[1]
    return out


def paired_de(tumor: pd.DataFrame, normal: pd.DataFrame, min_pairs: int = MIN_PAIRS) -> pd.DataFrame:
    """Vectorized paired t-test (tumor - normal) on matched columns, rows = features."""
    feats = tumor.index.intersection(normal.index)
    cases = tumor.columns.intersection(normal.columns)
    d = tumor.loc[feats, cases].to_numpy(float) - normal.loc[feats, cases].to_numpy(float)
    n = np.isfinite(d).sum(axis=1)
    mean = np.nanmean(np.where(n[:, None] > 0, d, 0), axis=1) if d.size else np.array([])
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = np.nansum(d, axis=1) / np.maximum(n, 1)
        sd = np.sqrt(np.nansum((d - mean[:, None]) ** 2, axis=1) / np.maximum(n - 1, 1))
        t = mean / (sd / np.sqrt(np.maximum(n, 1)))
        p = 2 * stats.t.sf(np.abs(t), df=np.maximum(n - 1, 1))
    p[(n < min_pairs) | ~np.isfinite(t)] = np.nan
    out = pd.DataFrame({"feature": feats, "n_pairs": n, "log2fc": mean, "t": t, "pvalue": p})
    out = out[out["pvalue"].notna()].copy()
    out["padj"] = bh(out["pvalue"])
    return out.reset_index(drop=True)


# ---------------------------------------------------------------- somatic mutations
def cds_lengths() -> pd.Series:
    """Union of CDS intervals per gene (bp) from the GENCODE v22 GTF."""
    iv: dict[str, list[tuple[int, int]]] = {}
    with gzip.open(EXTRA / "gencode.v22.annotation.gtf.gz", "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.split("\t", 8)
            if p[2] != "CDS":
                continue
            m = re.search(r'gene_name "([^"]+)"', p[8])
            if m:
                iv.setdefault(m.group(1), []).append((int(p[3]), int(p[4])))
    out = {}
    for g, xs in iv.items():
        xs.sort()
        total, (s, e) = 0, xs[0]
        for a, b in xs[1:]:
            if a > e:
                total += e - s + 1
                s, e = a, b
            else:
                e = max(e, b)
        out[g] = total + e - s + 1
    return pd.Series(out, name="cds_bp")


def mutations() -> pd.DataFrame:
    hg = load_hgnc()
    approved = set(hg["symbol"])
    prev = {}
    for s, al in zip(hg["symbol"], hg["prev_symbol"]):
        if isinstance(al, str):
            for a in al.split("|"):
                prev.setdefault(a.strip(), s)
    fix = lambda g: g if g in approved else prev.get(g, g)  # noqa: E731
    cds = cds_lengths()
    cds.index = [fix(g) for g in cds.index]
    cds = cds.groupby(level=0).max()
    discovery = set(read_matrix(LAYERS["rna"][0]).columns) | set(read_matrix(LAYERS["protein"][0]).columns)
    rows, cohorts = [], {}
    for key in ("cptac", "tcga"):
        m = pd.read_parquet(EXTRA / f"mutations_{key}.parquet")
        samples = json.loads((EXTRA / f"mutations_{key}.samples.json").read_text())
        patients = {s.rsplit("-", 1)[0] if key == "cptac" else s[:12] for s in samples}
        if key == "cptac":
            patients &= discovery          # only the CPTAC discovery patients analysed everywhere else
        m = m[m["mutation_type"].isin(NONSYN) & m["patient"].isin(patients) & m["gene"].notna()].copy()
        m["gene"] = m["gene"].map(fix)
        burden = m.groupby("patient").size()
        hyper = set(burden[burden > HYPERMUTATOR].index)
        keep = patients - hyper
        m = m[m["patient"].isin(keep)]
        genes_with_cds = m[m["gene"].isin(cds.index)]
        rate = len(genes_with_cds) / len(keep) / cds.sum()   # nonsynonymous mutations per coding bp per tumour
        k = m.groupby("gene")["patient"].nunique()
        cohorts[key] = {"patients": len(keep), "excluded_hypermutators": len(hyper), "rate_per_mb": rate * 1e6,
                        "k": k}
    genes = sorted(set(cohorts["cptac"]["k"].index) | set(cohorts["tcga"]["k"].index))
    for g in genes:
        if g not in cds.index:
            continue
        L = cds[g]
        k_tot, n_tot, exp_tot, rec = 0, 0, 0.0, {"gene": g, "cds_bp": int(L)}
        for key, c in cohorts.items():
            k = int(c["k"].get(g, 0))
            p0 = 1 - np.exp(-c["rate_per_mb"] / 1e6 * L)      # P(tumour carries >= 1 nonsyn. mutation by chance)
            rec[f"{key}_mutated"] = k
            rec[f"{key}_freq"] = k / c["patients"]
            k_tot += k
            n_tot += c["patients"]
            exp_tot += p0 * c["patients"]
        rec["mutated"] = k_tot
        rec["patients"] = n_tot
        rec["expected"] = exp_tot
        rec["obs_exp"] = k_tot / exp_tot if exp_tot > 0 else np.nan
        rec["pvalue"] = float(stats.binom.sf(k_tot - 1, n_tot, exp_tot / n_tot)) if k_tot else 1.0
        rows.append(rec)
    df = pd.DataFrame(rows)
    df["padj"] = bh(df["pvalue"])
    df["freq"] = df["mutated"] / df["patients"]
    df = df.sort_values("pvalue").reset_index(drop=True)
    df.to_parquet(PROCESSED / "layer_mutations.parquet", index=False)
    meta = {k: {kk: (vv if not isinstance(vv, pd.Series) else None) for kk, vv in v.items() if kk != "k"}
            for k, v in cohorts.items()}
    (PROCESSED / "layer_mutations_meta.json").write_text(json.dumps(meta, indent=1))
    return df


# ---------------------------------------------------------------- copy number cis effect
def _rowwise_spearman(a: pd.DataFrame, b: pd.DataFrame, min_n: int = 30) -> pd.DataFrame:
    genes = a.index.intersection(b.index)
    cases = a.columns.intersection(b.columns)
    A = a.loc[genes, cases].rank(axis=1).to_numpy(float)
    B = b.loc[genes, cases].rank(axis=1).to_numpy(float)
    mask = np.isfinite(A) & np.isfinite(B)
    n = mask.sum(axis=1)
    A = np.where(mask, A, np.nan)
    B = np.where(mask, B, np.nan)
    # re-rank within the jointly observed values
    A = pd.DataFrame(A).rank(axis=1).to_numpy()
    B = pd.DataFrame(B).rank(axis=1).to_numpy()
    with np.errstate(invalid="ignore", divide="ignore"):
        am = A - np.nanmean(A, axis=1, keepdims=True)
        bm = B - np.nanmean(B, axis=1, keepdims=True)
        rho = np.nansum(am * bm, axis=1) / np.sqrt(np.nansum(am ** 2, axis=1) * np.nansum(bm ** 2, axis=1))
        t = rho * np.sqrt((n - 2) / np.maximum(1 - rho ** 2, 1e-12))
        p = 2 * stats.t.sf(np.abs(t), df=np.maximum(n - 2, 1))
    p[n < min_n] = np.nan
    rho[n < min_n] = np.nan
    return pd.DataFrame({"gene": genes, "rho": rho, "n": n, "pvalue": p})


def cis_dosage() -> pd.DataFrame:
    cnv = read_matrix("HS_CPTAC_CCRCC_CNV_gene_Tumor.cct")
    parts = []
    for layer in ("rna", "protein"):
        r = _rowwise_spearman(cnv, read_matrix(LAYERS[layer][0]))
        r["padj"] = bh(r["pvalue"])
        parts.append(r.set_index("gene").add_prefix(f"{layer}_"))
    df = pd.concat(parts, axis=1).reset_index().rename(columns={"index": "gene"})
    df = to_hgnc(df, "gene", rank_by="rna_n")
    df.to_parquet(PROCESSED / "layer_cis_dosage.parquet", index=False)
    return df


# ---------------------------------------------------------------- DepMap
def ccrcc_lines() -> pd.DataFrame:
    """Renal cell carcinoma lines with a DepMap id, excluding papillary / chromophobe / translocation subtypes
    (generic 'renal cell carcinoma' lines such as 786-O, A-498 and 769-P are VHL-null clear cell models)."""
    rows, cur = [], {}
    for line in (EXTRA / "cellosaurus_rcc.txt").read_text().splitlines():
        if line.startswith("ID   "):
            cur = {"name": line[5:].strip(), "di": [], "ach": None}
        elif line.startswith("DI   "):
            cur["di"].append(line.split(";")[-1].strip())
        elif line.startswith("DR   DepMap;"):
            cur["ach"] = line.split(";")[1].strip()
        elif line.startswith("//") and cur:
            rows.append(cur)
            cur = {}
    df = pd.DataFrame(rows)
    bad = ("Papillary", "chromophobe", "Chromophobe", "translocation", "TFE3", "leiomyomatosis", "MiT family")
    df["clear_cell"] = df["di"].map(lambda d: any("Clear cell" in x or x == "Renal cell carcinoma" for x in d)
                                    and not any(b in x for x in d for b in bad))
    return df[df["clear_cell"] & df["ach"].notna()].drop_duplicates("ach")


def selective_dependency(dep: pd.DataFrame, cfg: dict | None = None) -> pd.Series:
    """The one definition of "selective DepMap dependency of ccRCC" used by the scores, the graph and discovery.

    ccRCC lines depend on the gene more than other lines (BH q < padj, positive selectivity), the gene is not common
    essential, and at least `min_dependent_lines` ccRCC lines are dependent (Chronos < -0.5). The last condition was
    added after looking at the M14 results: a "selective" call carried by one or two lines out of 24 is too fragile
    (documented as a post-hoc criterion change in evals/history.md)."""
    dc = cfg or load_cfg()["convergence"]["dependency"]
    ess = dep["common_essential"].astype("boolean").fillna(False).astype(bool)
    n_dep = dep["n_ccrcc_dependent"] if "n_ccrcc_dependent" in dep else dep["frac_ccrcc_dependent"] * dep["n_ccrcc_lines"]
    return ((dep["n_ccrcc_lines"] >= dc["min_lines"]) & (dep["padj"] < dc["padj"]) & (dep["selectivity"] > 0)
            & (n_dep.round() >= dc.get("min_dependent_lines", 1)) & ~ess).fillna(False).astype(bool)


def dependency() -> pd.DataFrame:
    ge = pd.read_csv(EXTRA / "depmap_26Q1_gene_effect.csv", index_col=0, dtype={0: str})
    ge.columns = [c.split(" (")[0] for c in ge.columns]
    ge = ge.astype("float32")
    lines = ccrcc_lines()
    rcc = [a for a in lines["ach"] if a in ge.index]
    other = ge.index.difference(rcc)
    R, O = ge.loc[rcc], ge.loc[other]
    res = stats.ttest_ind(R, O, axis=0, equal_var=False, nan_policy="omit")
    df = pd.DataFrame({
        "gene": ge.columns, "n_ccrcc_lines": R.notna().sum().to_numpy(), "mean_effect_ccrcc": R.mean().to_numpy(),
        "mean_effect_other": O.mean().to_numpy(), "frac_ccrcc_dependent": (R < -0.5).mean().to_numpy(),
        "n_ccrcc_dependent": (R < -0.5).sum().to_numpy(),
        "frac_all_dependent": (ge < -0.5).mean().to_numpy(), "pvalue": np.asarray(res.pvalue, float)})
    df["selectivity"] = df["mean_effect_other"] - df["mean_effect_ccrcc"]
    df["padj"] = bh(df["pvalue"])
    df["common_essential"] = df["frac_all_dependent"] >= 0.9
    df = to_hgnc(df, "gene", rank_by="n_ccrcc_lines")
    df.to_parquet(PROCESSED / "layer_dependency.parquet", index=False)
    (PROCESSED / "layer_dependency_lines.json").write_text(json.dumps(
        {"release": "DepMap Public 26Q1 (Chronos)", "ccrcc_lines": lines[lines["ach"].isin(rcc)][["name", "ach"]]
         .to_dict("records"), "other_lines": int(len(other))}, indent=1))
    return df


# ---------------------------------------------------------------- TCGA-KIRC replication + survival
def replication() -> pd.DataFrame:
    x = pd.read_csv(EXTRA / "tcga_kirc_hiseqv2.gz", sep="\t", index_col=0)
    tum = x[[c for c in x.columns if c.endswith("-01")]]
    nor = x[[c for c in x.columns if c.endswith("-11")]]
    tum.columns = [c[:12] for c in tum.columns]
    nor.columns = [c[:12] for c in nor.columns]
    tum = tum.loc[:, ~tum.columns.duplicated()]
    nor = nor.loc[:, ~nor.columns.duplicated()]
    de = paired_de(tum, nor).rename(columns={"feature": "gene"})
    surv = pd.read_csv(EXTRA / "tcga_kirc_survival.txt", sep="\t")
    surv = surv[surv["sample"].str.endswith("-01")].dropna(subset=["OS", "OS.time"])
    surv["patient"] = surv["sample"].str[:12]
    surv = surv.drop_duplicates("patient").set_index("patient")
    common = tum.columns.intersection(surv.index)
    T = surv.loc[common, "OS.time"].to_numpy(float)
    E = surv.loc[common, "OS"].to_numpy(int)
    X = tum[common]
    X = X[(X.mean(axis=1) >= 1) & (X.std(axis=1) > 0.1)]
    Z = ((X.T - X.mean(axis=1)) / X.std(axis=1)).to_numpy(float)
    hr, pv = cox_univariate(Z, T, E)
    cox = pd.DataFrame({"gene": X.index, "cox_log_hr_per_sd": hr, "cox_pvalue": pv})
    cox["cox_padj"] = bh(cox["cox_pvalue"])
    df = de.merge(cox, on="gene", how="outer")
    df["n_survival"] = len(common)
    # Xena's TCGA matrix uses 2016 symbols; re-key on current HGNC like the CPTAC tables.
    df = to_hgnc(df, "gene", rank_by="n_pairs")
    df.to_parquet(PROCESSED / "layer_tcga.parquet", index=False)
    return df


def cox_univariate(Z: np.ndarray, T: np.ndarray, E: np.ndarray, iters: int = 25, return_se: bool = False):
    """Vectorized univariate Cox regression (Breslow ties) for many covariates at once (Newton-Raphson).
    Z: n x g standardized covariates. Returns (beta, Wald p-value) per covariate, plus the standard error if asked."""
    order = np.argsort(-T, kind="mergesort")            # descending time: risk sets are prefixes
    T, E, Z = T[order], E[order], Z[order]
    beta = np.zeros(Z.shape[1])
    for _ in range(iters):
        w = np.exp(np.clip(Z * beta, -50, 50))
        s0 = np.cumsum(w, axis=0)
        s1 = np.cumsum(w * Z, axis=0)
        s2 = np.cumsum(w * Z * Z, axis=0)
        # tie handling: risk set of a time includes everyone with T >= t -> take the last index of each tied block
        last = np.r_[np.nonzero(np.diff(T))[0], len(T) - 1]
        idx = np.searchsorted(last, np.arange(len(T)))
        s0, s1, s2 = s0[last][idx], s1[last][idx], s2[last][idx]
        ev = E[:, None] == 1
        grad = np.where(ev, Z - s1 / s0, 0).sum(axis=0)
        info = np.where(ev, s2 / s0 - (s1 / s0) ** 2, 0).sum(axis=0)
        step = grad / np.maximum(info, 1e-12)
        beta = beta + step
        if np.nanmax(np.abs(step)) < 1e-8:
            break
    se = 1 / np.sqrt(np.maximum(info, 1e-12))
    p = 2 * stats.norm.sf(np.abs(beta / se))
    return (beta, p, se) if return_se else (beta, p)


# ---------------------------------------------------------------- phosphosites and kinase activity
def phosphosites() -> pd.DataFrame:
    t = pd.read_csv(EXTRA / "cptac_ccrcc_phosphosite_tumor.txt.gz", sep="\t", index_col=0)
    n = pd.read_csv(EXTRA / "cptac_ccrcc_phosphosite_normal.txt.gz", sep="\t", index_col=0)
    de = paired_de(t, n)
    parts = de["feature"].str.split("|", expand=True)
    hg = load_hgnc().dropna(subset=["ensembl_gene_id"]).drop_duplicates("ensembl_gene_id")
    sym = dict(zip(hg["ensembl_gene_id"], hg["symbol"]))
    de["gene"] = parts[0].str.split(".").str[0].map(sym)
    de["site"] = parts[2]
    de["flank"] = parts[3]
    de = de.dropna(subset=["gene"])
    de["site_id"] = de["gene"] + "_" + de["site"]
    de.to_parquet(PROCESSED / "layer_phosphosite.parquet", index=False)
    return de


def protein_kinases() -> set[str]:
    """HGNC genes whose name or gene group says kinase, minus regulatory / anchoring / inhibitor proteins
    (OmniPath 'enzymes' also include ligands such as EGF, which act through receptors, not on the site)."""
    hg = pd.read_csv(RAW / "hgnc_complete_set.txt", sep="\t", usecols=["symbol", "name", "gene_group"], dtype=str)
    txt = (hg["name"].fillna("") + " | " + hg["gene_group"].fillna("")).str.lower()
    kin = txt.str.contains("kinase")
    bad = hg["name"].fillna("").str.lower().str.contains(
        "regulatory subunit|anchor|inhibitor|interacting|binding|substrate|activator|adaptor|phosphatase")
    return set(hg.loc[kin & ~bad, "symbol"])


def kinase_activity(sites: pd.DataFrame | None = None, min_substrates: int = 5) -> pd.DataFrame:
    """KSEA (Casado et al. 2013): z = (mean log2FC of a kinase's substrate sites - mean of all sites) * sqrt(m) / sd.

    Site log2FC is first adjusted for the change of its protein (site log2FC - protein log2FC, CPTAC proteome), so a
    kinase is not called active only because its substrates are more abundant in tumors."""
    sites = phosphosites() if sites is None else sites
    s = sites.sort_values("padj").drop_duplicates("site_id").set_index("site_id")
    prot = pd.read_parquet(PROCESSED / "de_results.parquet")
    prot = prot[prot["layer"] == "protein"].dropna(subset=["hgnc_symbol"]).drop_duplicates("hgnc_symbol")
    plfc = prot.set_index("hgnc_symbol")["log2fc"]
    s = (s["log2fc"] - s["gene"].map(plfc)).dropna()   # sites whose protein is not quantified are not used
    sites["protein_adjusted_log2fc"] = sites["log2fc"] - sites["gene"].map(plfc)
    sites.to_parquet(PROCESSED / "layer_phosphosite.parquet", index=False)
    es = pd.read_csv(EXTRA / "omnipath_enzsub.tsv", sep="\t")
    es = es[(es["modification"] == "phosphorylation") & es["enzyme_genesymbol"].isin(protein_kinases())]
    es["site_id"] = es["substrate_genesymbol"] + "_" + es["residue_type"] + es["residue_offset"].astype(str)
    mu, sd = s.mean(), s.std()
    rows = []
    for kin, grp in es.groupby("enzyme_genesymbol"):
        sub = s.reindex(grp["site_id"].unique()).dropna()
        if len(sub) < min_substrates:
            continue
        z = (sub.mean() - mu) * np.sqrt(len(sub)) / sd
        rows.append({"gene": kin, "n_substrate_sites": len(sub), "mean_substrate_log2fc": sub.mean(), "z": z,
                     "pvalue": 2 * stats.norm.sf(abs(z))})
    df = pd.DataFrame(rows)
    df["padj"] = bh(df["pvalue"])
    df = df.sort_values("z", ascending=False).reset_index(drop=True)
    df.to_parquet(PROCESSED / "layer_kinase_activity.parquet", index=False)
    return df


# ---------------------------------------------------------------- tumour purity
def purity() -> pd.DataFrame:
    pur = pd.read_csv(EXTRA / "cptac_pancan_purity_estimate.tsv.gz", sep="\t")
    pur = pur[pur["Sample_ID"].str.endswith("-T")]
    pur["case"] = pur["Sample_ID"].str[:-2]
    pur = pur.drop_duplicates("case").set_index("case")["TumorPurity"]
    out = []
    for layer in ("rna", "protein"):
        tf, nf = LAYERS[layer]
        t, n = read_matrix(tf), read_matrix(nf)
        cases = t.columns.intersection(pur.index)
        # (1) correlation of tumour expression with purity
        P = pd.DataFrame(np.tile(pur[cases].to_numpy(float), (len(t), 1)), index=t.index, columns=cases)
        c = _rowwise_spearman(t[cases], P, min_n=30).set_index("gene")
        # (2) paired log2FC regressed on purity: d_i = a + b * purity_i ; extrapolated log2FC at purity 1 = a + b
        pc = cases.intersection(n.columns)
        d = (t.loc[:, pc] - n.reindex(columns=pc)).to_numpy(float)
        x = pur[pc].to_numpy(float)
        ok = np.isfinite(d)
        nn = ok.sum(axis=1)
        X = np.where(ok, x[None, :], np.nan)
        with np.errstate(invalid="ignore", divide="ignore"):
            xm = np.nanmean(X, axis=1, keepdims=True)
            ym = np.nanmean(np.where(ok, d, np.nan), axis=1, keepdims=True)
            sxx = np.nansum((X - xm) ** 2, axis=1)
            sxy = np.nansum((X - xm) * (np.where(ok, d, np.nan) - ym), axis=1)
            b = sxy / sxx
            a = ym[:, 0] - b * xm[:, 0]
            resid = np.where(ok, d, np.nan) - (a[:, None] + b[:, None] * X)
            s2 = np.nansum(resid ** 2, axis=1) / np.maximum(nn - 2, 1)
            se_pred = np.sqrt(s2 * (1 / np.maximum(nn, 1) + (1 - xm[:, 0]) ** 2 / sxx))
        fit = pd.DataFrame({"lfc_at_purity1": a + b, "lfc_at_purity1_se": se_pred, "lfc_purity_slope": b,
                            "n_pairs_purity": nn}, index=t.index)
        fit.loc[fit["n_pairs_purity"] < MIN_PAIRS, ["lfc_at_purity1", "lfc_at_purity1_se", "lfc_purity_slope"]] = np.nan
        part = pd.DataFrame({"purity_rho": c["rho"], "purity_pvalue": c["pvalue"]}).join(fit, how="outer")
        part["purity_padj"] = bh(part["purity_pvalue"])
        part.columns = [f"{layer}_{col}" for col in part.columns]
        out.append(part)
    df = pd.concat(out, axis=1)
    df.index.name = "dataset_symbol"
    df = df.reset_index()
    # previous / alias symbols map too (approved only dropped renamed genes such as NDUFA4L2 -> COXFA4L2)
    df = to_hgnc(df, "dataset_symbol", rank_by="rna_n_pairs_purity", out="hgnc_symbol")
    df["n_tumors_with_purity"] = int(len(pur.index.intersection(read_matrix(LAYERS["rna"][0]).columns)))
    df.to_parquet(PROCESSED / "layer_purity.parquet", index=False)
    return df


def run() -> dict:
    mut = mutations()
    cis = cis_dosage()
    dep = dependency()
    tcga = replication()
    sites = phosphosites()
    kin = kinase_activity(sites)
    pur = purity()
    summary = {
        "mutations": {"genes": int(len(mut)), "significant": int((mut["padj"] < 0.05).sum()),
                      "top": mut.head(10)[["gene", "mutated", "patients", "obs_exp", "padj"]].round(4).to_dict("records")},
        "cis_dosage": {"genes": int(len(cis)), "protein_cis_sig": int((cis["protein_padj"] < 0.05).sum())},
        "dependency": json.loads((PROCESSED / "layer_dependency_lines.json").read_text()) | {
            "selective_genes_padj05": int(((dep["padj"] < 0.05) & (dep["selectivity"] > 0)).sum())},
        "tcga": {"genes": int(tcga["log2fc"].notna().sum()), "pairs": int(tcga["n_pairs"].max()),
                 "survival_patients": int(tcga["n_survival"].max())},
        "phosphosites": {"sites": int(len(sites)), "significant": int((sites["padj"] < 0.05).sum())},
        "kinases": {"tested": int(len(kin)), "activated_padj05": int(((kin["padj"] < 0.05) & (kin["z"] > 0)).sum()),
                    "top": kin.head(8)[["gene", "n_substrate_sites", "z", "padj"]].round(4).to_dict("records")},
        "purity": {"tumors": int(pur["n_tumors_with_purity"].iloc[0])},
    }
    (PROCESSED / "layers_summary.json").write_text(json.dumps(summary, indent=1, default=float))
    return summary


if __name__ == "__main__":
    print(json.dumps(run(), indent=1, default=float))
