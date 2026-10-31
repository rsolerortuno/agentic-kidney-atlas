"""Convergence v3 components, target / biomarker profiles and weight-sensitivity analysis (deterministic, no LLM).

Changes from v2 (see evals/history.md for why):

1. Omics is continuous. Each layer contributes min(1, |log2FC| / cap) when it passes padj, instead of a 0/1 flag,
   so the ~550 genes that all had a "perfect" binary omics component are ordered by effect size.
2. Copy number counts only when it moves the gene product: loss/gain frequency x cis effect (Spearman of copy number
   with protein, or RNA when protein is missing). The 3p arm loss no longer gives the same bonus to ~490 genes.
3. Genetics uses somatic recurrence measured here (CPTAC + TCGA-KIRC, length-corrected) instead of the Open Targets
   somatic datatype, whose uniform cancer-census score (0.456) was attached to many unrelated genes.
4. Tractability comes from the Open Targets bulk table for every gene, not only for the genes chosen for enrichment.
5. New streams: DepMap selective dependency and TCGA-KIRC replication.

Two profiles answer two different questions and leave out literature and clinical-drug information on purpose, so they
can be evaluated against clinical outcomes without circularity:

- target    = tumour-up expression, purity-extrapolated tumour-intrinsic signal, selective dependency, driver status,
              kinase activity and predicted (not clinical) tractability.
- biomarker = expression change in either direction, external replication, prognosis (TCGA overall survival), grade or
              stage association, RNA-protein agreement and cis copy-number effect.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd

from .config import load
from .paths import PROCESSED, RAW

EXTRA = RAW / "extra"
LAYERS = ("rna", "protein", "phospho")


def _clip(x) -> np.ndarray:
    """Clip to [0, 1] as a plain array (missing -> 0)."""
    return np.nan_to_num(np.clip(np.asarray(x, dtype=float), 0, 1), nan=0.0)


def _map(index: pd.Index, series: pd.Series) -> np.ndarray:
    return np.asarray(index.map(series), dtype=float)


def _by_symbol(path, col="gene") -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    d = pd.read_parquet(path)
    d = d.dropna(subset=[col])
    return d.drop_duplicates(col).set_index(col)


def tractability_table() -> pd.DataFrame:
    """Open Targets 26.09 tractability for every target: clinical precedence vs predicted SM / AB tractability."""
    parts = sorted(EXTRA.glob("ot_target_tractability_*.parquet"))
    if not parts:
        return pd.DataFrame(columns=["clinical", "predicted_sm_ab"])
    t = pd.concat([pd.read_parquet(p) for p in parts], ignore_index=True)
    t = t[t["value"]]
    clinical = t["category"].isin(["Approved Drug", "Advanced Clinical", "Phase 1 Clinical"])
    # strong predicted buckets only: SM with a high-quality ligand or pocket, AB with high-confidence surface /
    # secreted location (broad buckets such as "UniProt SigP or TMHMM" cover a third of the genome)
    predicted = (((t["modality"] == "SM") & t["category"].isin(["High-Quality Ligand", "High-Quality Pocket"]))
                 | ((t["modality"] == "AB") & (t["category"] == "UniProt loc high conf")))
    hg = pd.read_csv(RAW / "hgnc_complete_set.txt", sep="\t", usecols=["symbol", "ensembl_gene_id"], dtype=str).dropna()
    sym = dict(zip(hg["ensembl_gene_id"], hg["symbol"]))
    t["gene"] = t["targetId"].map(sym)
    t = t.dropna(subset=["gene"])
    genes = sorted(set(t["gene"]))
    out = pd.DataFrame(index=pd.Index(genes, name="gene"))
    out["clinical"] = out.index.isin(t.loc[clinical, "gene"])
    out["predicted_sm_ab"] = out.index.isin(t.loc[predicted, "gene"])
    mods = t[predicted].groupby("gene")["modality"].agg(lambda s: ",".join(sorted(set(s))))
    out["predicted_modalities"] = mods.reindex(out.index)
    return out


def components(t: pd.DataFrame) -> pd.DataFrame:
    """Add every v3 component (all in [0, 1]) to the gene universe table `t` (one row per HGNC symbol)."""
    cfg = load()["convergence"]
    caps, lw, thr = cfg["omics_caps"], cfg["omics_layer_weights"], load()["omics"]
    t = t.copy().set_index("gene", drop=False)

    # ---- omics: continuous, gated by significance; "up" version for the target profile
    s, s_up = {}, {}
    for layer in LAYERS:
        lfc, padj = t[f"{layer}_log2fc"], t[f"{layer}_padj"]
        mag = _clip(lfc.abs() / caps[layer]) * (padj < thr["padj"]).to_numpy(float)
        s[layer] = mag
        s_up[layer] = mag * (lfc > 0).to_numpy(float)
    conc = t["rna_sig"] & t["protein_sig"] & (np.sign(t["rna_log2fc"]) == np.sign(t["protein_log2fc"]))
    s["concordance"] = np.minimum(s["rna"], s["protein"]) * conc.to_numpy(float)
    s_up["concordance"] = np.minimum(s_up["rna"], s_up["protein"]) * conc.to_numpy(float)
    tot = sum(lw.values())
    t["c_omics"] = sum(lw[k] * s[k] for k in lw) / tot
    t["c_omics_up"] = sum(lw[k] * s_up[k] for k in lw) / tot

    # ---- copy number x cis effect
    cis = _by_symbol(PROCESSED / "layer_cis_dosage.parquet")
    if len(cis):
        rho_p = cis["protein_rho"].where(cis["protein_padj"] < 0.05)
        rho_r = cis["rna_rho"].where(cis["rna_padj"] < 0.05)
        rho = rho_p.fillna(rho_r)
        t["cis_rho"] = _map(t.index, rho)
        t["c_cnv"] = t["cnv_freq"].fillna(0).to_numpy(float) * _clip(t["cis_rho"] / cfg["cis_rho_saturation"])
    else:
        t["cis_rho"], t["c_cnv"] = np.nan, t["cnv_freq"].fillna(0)

    # ---- somatic recurrence (genetics); Open Targets germline genetic association kept
    mut = _by_symbol(PROCESSED / "layer_mutations.parquet")
    t["mutation_padj"] = _map(t.index, mut["padj"]) if len(mut) else np.nan
    t["mutation_freq"] = _map(t.index, mut["freq"]) if len(mut) else np.nan
    q = t["mutation_padj"].clip(lower=1e-300)
    t["c_mutation"] = _clip(-np.log10(q) / cfg["mutation_log10q_saturation"]) * (q < 0.05).to_numpy(float)
    t["c_genetics"] = np.maximum(t["ot_genetic"].fillna(0).to_numpy(float), t["c_mutation"].to_numpy(float))

    # ---- tractability (all genes)
    tr = tractability_table()
    clin = np.asarray(t.index.map(tr["clinical"]).fillna(False), bool) if len(tr) else np.zeros(len(t), bool)
    pred = np.asarray(t.index.map(tr["predicted_sm_ab"]).fillna(False), bool) if len(tr) else np.zeros(len(t), bool)
    t["has_clinical_drug_any_indication"] = clin
    t["predicted_tractable_sm_ab"] = pred
    t["c_tractability"] = np.where(clin, 1.0, np.where(pred, 0.5, 0.0))
    # target profile: predicted modality evidence only. Clinical precedence is drug history and would leak into the
    # evaluation against clinical targets (a first version used pred | clin; see evals/history.md, M13).
    t["c_tractability_predicted"] = pred.astype(float)

    # ---- DepMap selective dependency
    dep = _by_symbol(PROCESSED / "layer_dependency.parquet")
    dc = cfg["dependency"]
    if len(dep):
        ok = (dep["n_ccrcc_lines"] >= dc["min_lines"]) & (dep["padj"] < dc["padj"]) & ~dep["common_essential"]
        c = pd.Series(_clip(dep["selectivity"] / dc["selectivity_saturation"]) * _clip(
            dep["frac_ccrcc_dependent"] / dc["fraction_saturation"]) * ok.to_numpy(float), index=dep.index)
        t["c_dependency"] = np.nan_to_num(_map(t.index, c))
        t["depmap_mean_effect_ccrcc"] = _map(t.index, dep["mean_effect_ccrcc"])
        t["depmap_selectivity"] = _map(t.index, dep["selectivity"])
    else:
        t["c_dependency"] = 0.0

    # ---- TCGA-KIRC replication and prognosis
    tc = _by_symbol(PROCESSED / "layer_tcga.parquet")
    if len(tc):
        lfc = _map(t.index, tc["log2fc"])
        padj = _map(t.index, tc["padj"])
        ref = t["rna_log2fc"].where(t["rna_sig"], t["protein_log2fc"].where(t["protein_sig"])).to_numpy(float)
        same = np.sign(lfc) == np.sign(ref)          # False when the gene is not significant in CPTAC (ref NaN)
        t["tcga_log2fc"], t["tcga_padj"] = lfc, padj
        t["c_replication"] = _clip(np.abs(lfc) / caps["rna"]) * ((padj < 0.05) & same)
        hr = _map(t.index, tc["cox_log_hr_per_sd"])
        hp = _map(t.index, tc["cox_padj"])
        t["tcga_cox_log_hr"], t["tcga_cox_padj"] = hr, hp
        t["c_prognosis"] = _clip(np.abs(hr) / cfg["prognosis_log_hr_saturation"]) * (hp < 0.05)
    else:
        t["c_replication"], t["c_prognosis"] = 0.0, 0.0

    # ---- kinase activity (protein-adjusted KSEA)
    kin = _by_symbol(PROCESSED / "layer_kinase_activity.parquet")
    if len(kin):
        c = pd.Series(_clip(kin["z"] / cfg["kinase_z_saturation"]) * ((kin["padj"] < 0.1) & (kin["z"] > 0)).to_numpy(float),
                      index=kin.index)
        t["c_kinase"] = np.nan_to_num(_map(t.index, c))
        t["kinase_z"] = _map(t.index, kin["z"])
    else:
        t["c_kinase"] = 0.0

    # ---- purity: up-regulation that survives extrapolation to a pure tumour
    pur = _by_symbol(PROCESSED / "layer_purity.parquet", col="hgnc_symbol")
    if len(pur):
        a = _map(t.index, pur["rna_lfc_at_purity1"])
        se = _map(t.index, pur["rna_lfc_at_purity1_se"])
        t["rna_lfc_at_purity1"], t["rna_lfc_at_purity1_se"] = a, se
        t["rna_purity_rho"] = _map(t.index, pur["rna_purity_rho"])
        robust = (a - 1.96 * se) > 0
        t["c_tumor_intrinsic"] = _clip(a / caps["rna"]) * robust * (t["rna_log2fc"] > 0).to_numpy()
    else:
        t["c_tumor_intrinsic"] = 0.0

    # ---- grade / stage association and RNA-protein agreement (biomarker profile)
    clinp = PROCESSED / "clinical_associations.parquet"
    if clinp.exists():
        ca = pd.read_parquet(clinp)
        ca = ca[ca["padj"] < 0.05]
        best = ca.assign(a=ca["rho"].abs()).groupby("hgnc_symbol")["a"].max()
        t["c_grade"] = np.nan_to_num(_map(t.index, pd.Series(_clip(best / cfg["grade_rho_saturation"]), index=best.index)))
    else:
        t["c_grade"] = 0.0
    corr = _by_symbol(PROCESSED / "rna_protein_correlation.parquet", col="hgnc_symbol")
    if len(corr):
        rho = corr["spearman_rho"].where(corr["padj"] < 0.05)
        t["c_measurable"] = np.nan_to_num(_map(t.index, pd.Series(_clip(rho / cfg["rna_protein_rho_saturation"]),
                                                                  index=rho.index)))
    else:
        t["c_measurable"] = 0.0
    return t.reset_index(drop=True)


def profile_scores(t: pd.DataFrame) -> pd.DataFrame:
    prof = load()["profiles"]
    for name, w in prof.items():
        t[f"score_{name}"] = sum(wt * t[f"c_{k}"] for k, wt in w.items())
    return t


def weighted(t: pd.DataFrame, weights: dict) -> np.ndarray:
    return sum(w * t[f"c_{k}"].to_numpy(float) for k, w in weights.items())


def sensitivity(t: pd.DataFrame, weights: dict, n: int | None = None, top: int = 200, seed: int | None = None) -> pd.DataFrame:
    """Rank stability under random re-weighting: weights ~ Dirichlet(kappa * w / sum(w)) * sum(w).

    Returns, for the `top` genes of the base ranking, the median rank, 5th-95th percentile and the probability of
    staying in the top 10 / 25 / 100 across draws."""
    sc = load()["sensitivity"]
    n = n or sc["draws"]
    rng = np.random.default_rng(sc["seed"] if seed is None else seed)
    keys = list(weights)
    w = np.array([weights[k] for k in keys], float)
    C = np.column_stack([t[f"c_{k}"].to_numpy(float) for k in keys])
    base = C @ w
    order = np.argsort(-base, kind="stable")
    focus = order[:top]
    draws = rng.dirichlet(sc["kappa"] * w / w.sum(), size=n) * w.sum()
    ranks = np.empty((n, len(focus)), dtype=np.int32)
    for i, wd in enumerate(draws):
        s = C @ wd
        # rank of each focus gene = 1 + number of genes with a strictly higher score
        srt = np.sort(s)
        ranks[i] = len(s) - np.searchsorted(srt, s[focus], side="right") + 1
    genes = t["gene"].to_numpy()[focus]
    return pd.DataFrame({
        "gene": genes, "base_rank": np.arange(1, len(focus) + 1), "base_score": base[focus],
        "rank_median": np.median(ranks, axis=0), "rank_p05": np.percentile(ranks, 5, axis=0),
        "rank_p95": np.percentile(ranks, 95, axis=0), "p_top10": (ranks <= 10).mean(axis=0),
        "p_top25": (ranks <= 25).mean(axis=0), "p_top100": (ranks <= 100).mean(axis=0)})


def write_profiles(t: pd.DataFrame) -> dict:
    cfg = load()
    out = {}
    tables = {"convergence": cfg["convergence"]["stream_weights"], **cfg["profiles"]}
    for name, w in tables.items():
        sens = sensitivity(t, w)
        sens.to_parquet(PROCESSED / f"sensitivity_{name}.parquet", index=False)
        out[name] = {"robust_top10": sens[sens["p_top10"] >= 0.8]["gene"].tolist(),
                     "top10": sens.head(10)[["gene", "rank_p05", "rank_p95", "p_top10"]].to_dict("records")}
    cols = ["gene"] + [c for c in t.columns if c.startswith("c_") or c.startswith("score_")]
    t[cols].to_parquet(PROCESSED / "scores_v3.parquet", index=False)
    (PROCESSED / "sensitivity_summary.json").write_text(json.dumps(out, indent=1, default=float))
    return out
