"""Within-tumour stratification: genetic subtypes (VHL, PBRM1, BAP1, SETD2) and the ccA / ccB expression subtype.

Deterministic, no LLM. Tumour-versus-normal analyses find genes that separate ccRCC from kidney; these find genes that
separate ccRCC tumours from each other, which is what a subgroup-specific target or a prognostic biomarker needs.

1. Genetic subtypes. A tumour is "mutant" for a driver when it carries a nonsynonymous mutation (the same calls as the
   recurrence layer: cBioPortal, CPTAC discovery patients and TCGA-KIRC PanCancer Atlas). Mutant vs wild-type tumours
   are compared with Welch t-tests on CPTAC RNA and protein and on TCGA-KIRC RNA, BH within driver x layer. A
   subgroup-specific target candidate must be higher in mutant tumours in CPTAC RNA and protein (q < 0.1), replicate in
   TCGA (same direction, q < 0.1), and be up in tumour vs normal overall (so it is a tumour antigen in that subgroup,
   not a gene that is merely less lost).
2. ccA / ccB. ClearCode34 (Brooks et al. 2014) splits ccRCC into ccA (better prognosis) and ccB. Each tumour gets
   score = mean z(ccA genes) - mean z(ccB genes) and is called ccA when the score is positive (a signature-score
   approximation of the published centroid classifier). Check: ccB must have worse overall survival in TCGA. Then every
   gene is compared ccB vs ccA in CPTAC RNA and TCGA RNA. ClearCode34 genes themselves are flagged, and they get no
   ccA/ccB component in the scores, because they define the split.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
from scipy import stats

from .layers import EXTRA, NONSYN, bh, cox_univariate
from .omics import LAYERS, read_matrix
from .outcomes import CLEARCODE34, clearcode34
from .paths import PROCESSED
from .safety_clinical import GRADE

DRIVERS = ["VHL", "PBRM1", "BAP1", "SETD2"]
Q = 0.1
MIN_GROUP = 8


def _tcga_tumour() -> pd.DataFrame:
    x = pd.read_csv(EXTRA / "tcga_kirc_hiseqv2.gz", sep="\t", index_col=0)
    t = x[[c for c in x.columns if c.endswith("-01")]]
    t.columns = [c[:12] for c in t.columns]
    return t.loc[:, ~t.columns.duplicated()]


def _mutant_sets() -> dict[str, dict[str, tuple[set, set]]]:
    """{cohort: {driver: (mutant patients, sequenced patients)}}."""
    out = {}
    for key in ("cptac", "tcga"):
        m = pd.read_parquet(EXTRA / f"mutations_{key}.parquet")
        samples = json.loads((EXTRA / f"mutations_{key}.samples.json").read_text())
        seq = {s.rsplit("-", 1)[0] if key == "cptac" else s[:12] for s in samples}
        m = m[m["mutation_type"].isin(NONSYN)]
        out[key] = {g: (set(m.loc[m["gene"] == g, "patient"]) & seq, seq) for g in DRIVERS}
    return out


def welch(x: pd.DataFrame, a: list[str], b: list[str]) -> pd.DataFrame:
    """Rows = genes. Mean difference (a - b) and Welch t-test; genes with < MIN_GROUP values per group are dropped."""
    A, B = x[a].to_numpy(float), x[b].to_numpy(float)
    na, nb = np.isfinite(A).sum(axis=1), np.isfinite(B).sum(axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        r = stats.ttest_ind(A, B, axis=1, equal_var=False, nan_policy="omit")
        diff = np.nanmean(A, axis=1) - np.nanmean(B, axis=1)
    out = pd.DataFrame({"gene": x.index, "diff": diff, "pvalue": np.asarray(r.pvalue, float), "n_a": na, "n_b": nb})
    out = out[(out["n_a"] >= MIN_GROUP) & (out["n_b"] >= MIN_GROUP) & np.isfinite(out["pvalue"])].copy()
    out["padj"] = bh(out["pvalue"])
    return out


def genetic_subtypes() -> pd.DataFrame:
    muts = _mutant_sets()
    data = {"cptac_rna": read_matrix(LAYERS["rna"][0]), "cptac_protein": read_matrix(LAYERS["protein"][0]),
            "tcga_rna": _tcga_tumour()}
    rows = []
    for g in DRIVERS:
        for layer, x in data.items():
            mut, seq = muts["tcga" if layer.startswith("tcga") else "cptac"][g]
            cols = [c for c in x.columns if c in seq]
            a = [c for c in cols if c in mut]
            b = [c for c in cols if c not in mut]
            if len(a) < MIN_GROUP or len(b) < MIN_GROUP:
                continue
            r = welch(x, a, b)
            r["driver"], r["layer"] = g, layer
            rows.append(r)
    df = pd.concat(rows, ignore_index=True)
    df.to_parquet(PROCESSED / "subtype_genetic_de.parquet", index=False)
    return df


def subgroup_targets(df: pd.DataFrame) -> pd.DataFrame:
    """Genes higher in mutant tumours in CPTAC RNA + protein, replicated in TCGA, and tumour-up overall."""
    from .scores import tractability_table
    u = pd.read_parquet(PROCESSED / "gene_universe.parquet").set_index("gene")
    tr = tractability_table()
    out = []
    for g, d in df.groupby("driver"):
        w = d.pivot_table(index="gene", columns="layer", values=["diff", "padj"])
        w.columns = [f"{b}_{a}" for a, b in w.columns]
        need = [f"{l}_{k}" for l in ("cptac_rna", "cptac_protein", "tcga_rna") for k in ("diff", "padj")]
        if not all(c in w for c in need):
            continue
        ok = ((w["cptac_rna_diff"] > 0) & (w["cptac_rna_padj"] < Q) & (w["cptac_protein_diff"] > 0)
              & (w["cptac_protein_padj"] < Q) & (w["tcga_rna_diff"] > 0) & (w["tcga_rna_padj"] < Q))
        hit = w[ok].copy()
        hit["driver"] = g
        hit["tumour_up"] = hit.index.map(lambda s: bool((u.get("rna_log2fc", pd.Series()).get(s, 0) or 0) > 0
                                                        and (u.get("rna_padj", pd.Series()).get(s, 1) or 1) < 0.05))
        hit["predicted_tractable"] = hit.index.map(lambda s: bool(tr["predicted_sm_ab"].get(s, False)) if len(tr) else False)
        hit["clinical_drug"] = hit.index.map(lambda s: bool(tr["clinical"].get(s, False)) if len(tr) else False)
        out.append(hit.reset_index())
    t = pd.concat(out, ignore_index=True) if out else pd.DataFrame()
    if len(t):
        t["score"] = (t["cptac_protein_diff"].clip(0, 1.5) / 1.5 + t["tcga_rna_diff"].clip(0, 2) / 2) / 2
        t = t.sort_values(["driver", "score"], ascending=[True, False])
    t.to_parquet(PROCESSED / "subtype_targets.parquet", index=False)
    return t


# ---------------------------------------------------------------- ccA / ccB
def ccab_calls(x: pd.DataFrame) -> pd.Series:
    """Signature score mean z(ccA) - mean z(ccB) per tumour (columns of x); > 0 = ccA."""
    from .literature_entities import resolve_gene
    z = (x.T - x.mean(axis=1)) / x.std(axis=1).replace(0, np.nan)
    def genes(lst):
        return [g for g in ((resolve_gene(s) or {"UNG2": "UNG"}.get(s)) for s in lst) if g and g in z.columns]
    return z[genes(CLEARCODE34["ccA"])].mean(axis=1) - z[genes(CLEARCODE34["ccB"])].mean(axis=1)


def ccab() -> tuple[pd.DataFrame, dict]:
    cp, tc = read_matrix(LAYERS["rna"][0]), _tcga_tumour()
    cpp = read_matrix(LAYERS["protein"][0])
    s_cp, s_tc = ccab_calls(cp), ccab_calls(tc)
    rows = []
    for layer, x, s in (("cptac_rna", cp, s_cp), ("cptac_protein", cpp, s_cp), ("tcga_rna", tc, s_tc)):
        b = [c for c in x.columns if c in s.index and s[c] < 0]
        a = [c for c in x.columns if c in s.index and s[c] > 0]
        r = welch(x, b, a)        # ccB - ccA
        r["layer"] = layer
        rows.append(r)
    df = pd.concat(rows, ignore_index=True)
    cc = set(clearcode34())
    df["clearcode34_gene"] = df["gene"].isin(cc)
    df.to_parquet(PROCESSED / "subtype_ccab_de.parquet", index=False)
    # positive control: ccB should have worse overall survival in TCGA
    surv = pd.read_csv(EXTRA / "tcga_kirc_survival.txt", sep="\t")
    surv = surv[surv["sample"].str.endswith("-01")].dropna(subset=["OS", "OS.time"])
    surv["patient"] = surv["sample"].str[:12]
    surv = surv.drop_duplicates("patient").set_index("patient")
    common = s_tc.index.intersection(surv.index)
    ccb = (s_tc[common] < 0).astype(float).to_numpy()[:, None]
    beta, p = cox_univariate(ccb, surv.loc[common, "OS.time"].to_numpy(float), surv.loc[common, "OS"].to_numpy(int))
    cli = pd.read_csv(LAYERS_CLINICAL, sep="\t", index_col=0).iloc[1:]
    grade = pd.to_numeric(cli["Histologic_Grade"].astype(str).str.extract(r"(\d)")[0], errors="coerce")
    g = grade.reindex(s_cp.index)
    rho, rp = stats.spearmanr(s_cp[g.notna()], g.dropna())
    pd.DataFrame({"tumour": list(s_cp.index) + list(s_tc.index),
                  "cohort": ["CPTAC"] * len(s_cp) + ["TCGA-KIRC"] * len(s_tc),
                  "ccab_score": list(s_cp.values) + list(s_tc.values)}).assign(
        subtype=lambda d: np.where(d["ccab_score"] > 0, "ccA", "ccB")).to_parquet(PROCESSED / "subtype_ccab_calls.parquet",
                                                                                  index=False)
    meta = {"cptac": {"ccA": int((s_cp > 0).sum()), "ccB": int((s_cp < 0).sum())},
            "tcga": {"ccA": int((s_tc > 0).sum()), "ccB": int((s_tc < 0).sum())},
            "control_tcga_os_ccB_vs_ccA": {"hr": round(float(np.exp(beta[0])), 3), "pvalue": float(f"{p[0]:.3g}"),
                                           "n": int(len(common))},
            "control_cptac_grade_vs_score": {"spearman_rho": round(float(rho), 3), "pvalue": float(f"{rp:.3g}")}}
    return df, meta


LAYERS_CLINICAL = EXTRA.parent / "HS_CPTAC_CCRCC_CLI.tsi"


def run() -> dict:
    gen = genetic_subtypes()
    tg = subgroup_targets(gen)
    cab, meta = ccab()
    sig = gen[gen["padj"] < Q]
    w = cab.pivot_table(index="gene", columns="layer", values=["diff", "padj"])
    w.columns = [f"{b}_{a}" for a, b in w.columns]
    repl = w[(w["cptac_rna_padj"] < 0.05) & (w["tcga_rna_padj"] < 0.05)
             & (np.sign(w["cptac_rna_diff"]) == np.sign(w["tcga_rna_diff"]))]
    summary = {
        "drivers": {g: {"cptac_mutant": int(gen[(gen.driver == g) & (gen.layer == "cptac_rna")]["n_a"].max() or 0),
                        "tcga_mutant": int(gen[(gen.driver == g) & (gen.layer == "tcga_rna")]["n_a"].max() or 0),
                        "significant_by_layer": sig[sig.driver == g].groupby("layer").size().to_dict(),
                        "subgroup_targets": tg[tg.driver == g]["gene"].head(15).tolist() if len(tg) else [],
                        "subgroup_targets_tractable": tg[(tg.driver == g) & tg.predicted_tractable & tg.tumour_up]["gene"]
                        .head(10).tolist() if len(tg) else []} for g in DRIVERS},
        "ccab": meta | {"replicated_genes": int(len(repl)),
                        "top_ccB_up": repl[~repl.index.isin(clearcode34())].nlargest(10, "tcga_rna_diff").index.tolist(),
                        "top_ccA_up": repl[~repl.index.isin(clearcode34())].nsmallest(10, "tcga_rna_diff").index.tolist()},
    }
    (PROCESSED / "subtypes_summary.json").write_text(json.dumps(summary, indent=1, default=float))
    return summary


if __name__ == "__main__":
    print(json.dumps(run(), indent=1, default=float))
