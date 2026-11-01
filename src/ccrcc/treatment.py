"""Treatment-response biomarkers: which genes predict benefit from immunotherapy or from anti-angiogenic drugs?

Deterministic, no LLM. Two trials with pre-treatment tumour RNA-seq, both from the cBioPortal datahub (iAtlas
harmonized; see extra_download.py):

- IMmotion150 (McDermott et al. 2018, Nat Med): metastatic RCC randomised to atezolizumab + bevacizumab (n = 88),
  atezolizumab (86) or sunitinib (89). Outcomes: progression-free survival (PFS) and objective response (mRECIST).
- Choueiri et al. 2016 (Clin Cancer Res): nivolumab biomarker trial, 16 tumours, 3 responders. Too small for
  genome-wide claims; it is used only to check the direction of IMmotion150 immunotherapy hits.

Per gene and arm: univariate Cox PFS hazard ratio per SD of log2(TPM + 1) and the response AUC (responders vs the rest,
Mann-Whitney). The predictive question is an interaction: does the gene separate immunotherapy (atezolizumab +
bevacizumab) from the anti-angiogenic arm (sunitinib)? It is tested as the difference of the two arm-specific log
hazard ratios, z = (b_IO - b_TKI) / sqrt(se_IO^2 + se_TKI^2). All p-values are BH-adjusted within arm or contrast.

Positive controls (published in IMmotion150): the angiogenesis signature predicts longer PFS on sunitinib, and the
T-effector signature favours atezolizumab + bevacizumab over sunitinib. If the controls do not reproduce, the gene
calls are not trusted.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
from scipy import stats

from .layers import bh, cox_univariate
from .paths import PROCESSED, RAW

EXTRA = RAW / "extra"
SIGNATURES = {  # McDermott et al. 2018, Fig. 2
    "T_effector": ["CD8A", "EOMES", "PRF1", "IFNG", "CD274"],
    "angiogenesis": ["VEGFA", "KDR", "ESM1", "PECAM1", "ANGPTL4", "CD34"],
    "myeloid_inflammation": ["IL6", "CXCL1", "CXCL2", "CXCL3", "CXCL8", "PTGS2"],
}
ARMS = {"atezo_bev": "atezolizumab + bevacizumab", "atezo": "atezolizumab", "sunitinib": "sunitinib",
        "nivolumab": "nivolumab (Choueiri 2016)"}
MIN_MEAN_LOG2 = 1.0     # genes with mean log2(TPM + 1) below this are not tested
MAX_LOG_HR = 3.0        # |log HR per SD| above this means the Newton fit did not converge (separation)


def _read_clinical(key: str) -> pd.DataFrame:
    c = pd.read_csv(EXTRA / f"{key}_clinical_patient.txt", sep="\t", comment="#")
    s = pd.read_csv(EXTRA / f"{key}_clinical_sample.txt", sep="\t", comment="#")
    m = c.merge(s, on="PATIENT_ID", how="inner")
    m = m[m["SAMPLE_TREATMENT"].eq("Pre")].drop_duplicates("PATIENT_ID").set_index("SAMPLE_ID")
    m["pfs_time"] = pd.to_numeric(m["PFS_MONTHS"], errors="coerce")
    m["pfs_event"] = m["PFS_STATUS"].astype(str).str.startswith("1").astype(int)
    m["responder"] = m["RESPONSE"].map(lambda r: np.nan if pd.isna(r) else float(r in ("Complete Response",
                                                                                       "Partial Response")))
    if key == "immotion150":
        ici, other = m["ICI_RX"].fillna(""), m["NON_ICI_RX"].fillna("")
        m["arm"] = np.select([ici.eq("Atezolizumab") & other.eq("Bevacizumab"), ici.eq("Atezolizumab"),
                              other.eq("Sunitinib")], ["atezo_bev", "atezo", "sunitinib"], default="")
    else:
        m["arm"] = "nivolumab"
    return m[m["arm"] != ""]


def _read_expression(key: str, samples: pd.Index) -> pd.DataFrame:
    x = pd.read_csv(EXTRA / f"{key}_expression.txt", sep="\t")
    x = x.dropna(subset=["Hugo_Symbol"]).drop(columns=[c for c in ("Entrez_Gene_Id",) if c in x.columns])
    x = x.groupby("Hugo_Symbol").mean(numeric_only=True)
    x = x[[c for c in x.columns if c in samples]]
    if key == "immotion150":            # TPM -> log2(TPM + 1); the Choueiri file is already log2 upper-quartile
        x = np.log2(x.clip(lower=0) + 1)
    return x[x.mean(axis=1) >= MIN_MEAN_LOG2]


def _auc(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Per-row AUC of responders (y = 1) vs non-responders, with the Mann-Whitney normal-approximation p-value."""
    pos, neg = x[:, y == 1], x[:, y == 0]
    n1, n0 = pos.shape[1], neg.shape[1]
    ranks = stats.rankdata(np.hstack([pos, neg]), axis=1)
    u = ranks[:, :n1].sum(axis=1) - n1 * (n1 + 1) / 2
    auc = u / (n1 * n0)
    z = (u - n1 * n0 / 2) / np.sqrt(n1 * n0 * (n1 + n0 + 1) / 12)
    return auc, 2 * stats.norm.sf(np.abs(z))


def _arm_stats(x: pd.DataFrame, clin: pd.DataFrame, arm: str, study: str) -> pd.DataFrame:
    c = clin[(clin["arm"] == arm) & clin["pfs_time"].notna()]
    xs = x[c.index.intersection(x.columns)]
    # genes with little spread or mostly zero in this arm make the Cox fit diverge (|log HR| in the hundreds)
    xs = xs[(xs.std(axis=1) >= 0.25) & ((xs > 0).mean(axis=1) >= 0.5)]
    c = c.loc[xs.columns]
    Z = ((xs.T - xs.mean(axis=1)) / xs.std(axis=1).replace(0, np.nan)).fillna(0).to_numpy(float)
    beta, p, se = cox_univariate(Z, c["pfs_time"].to_numpy(float), c["pfs_event"].to_numpy(int), return_se=True)
    bad = ~np.isfinite(beta) | (np.abs(beta) > MAX_LOG_HR)     # not converged
    beta, se, p = np.where(bad, np.nan, beta), np.where(bad, np.nan, se), np.where(bad, np.nan, p)
    out = pd.DataFrame({"gene": xs.index, "study": study, "arm": arm, "n": len(c), "n_events": int(c["pfs_event"].sum()),
                        "cox_log_hr": beta, "cox_se": se, "cox_pvalue": p})
    resp = c["responder"].dropna()
    if resp.nunique() == 2 and min(resp.value_counts()) >= 3:
        auc, ap = _auc(xs[resp.index].to_numpy(float), resp.to_numpy(int))
        out["auc_response"], out["auc_pvalue"] = auc, ap
        out["n_responders"] = int(resp.sum())
    out["cox_padj"] = bh(out["cox_pvalue"])
    if "auc_pvalue" in out:
        out["auc_padj"] = bh(out["auc_pvalue"])
    return out


def signature_scores(x: pd.DataFrame) -> pd.DataFrame:
    z = (x.T - x.mean(axis=1)) / x.std(axis=1).replace(0, np.nan)
    return pd.DataFrame({k: z[[g for g in v if g in z.columns]].mean(axis=1) for k, v in SIGNATURES.items()})


def run() -> dict:
    rows, controls = [], {}
    clin = _read_clinical("immotion150")
    x = _read_expression("immotion150", clin.index)
    sig = signature_scores(x)
    for arm in ("atezo_bev", "atezo", "sunitinib"):
        rows.append(_arm_stats(x, clin, arm, "IMmotion150"))
        # signature positive controls: Cox on the standardized signature score in the arm
        c = clin[(clin["arm"] == arm)].loc[lambda d: d.index.isin(sig.index) & d["pfs_time"].notna()]
        s = sig.loc[c.index]
        s = (s - s.mean()) / s.std()
        b, p, se = cox_univariate(s.to_numpy(float), c["pfs_time"].to_numpy(float), c["pfs_event"].to_numpy(int),
                                  return_se=True)
        controls[arm] = {k: {"hr_per_sd": round(float(np.exp(bb)), 3), "pvalue": float(f"{pp:.3g}"), "se": float(ss)}
                         for k, bb, pp, ss in zip(s.columns, b, p, se)}
    arm = pd.concat(rows, ignore_index=True)

    # predictive contrast: immunotherapy combination vs anti-angiogenic arm
    io = arm[arm["arm"] == "atezo_bev"].set_index("gene")
    tki = arm[arm["arm"] == "sunitinib"].set_index("gene")
    g = io.index.intersection(tki.index)
    z = (io.loc[g, "cox_log_hr"] - tki.loc[g, "cox_log_hr"]) / np.sqrt(io.loc[g, "cox_se"] ** 2 + tki.loc[g, "cox_se"] ** 2)
    inter = pd.DataFrame({"gene": g, "study": "IMmotion150", "arm": "atezo_bev_vs_sunitinib",
                          "n": int(io["n"].iloc[0] + tki["n"].iloc[0]),
                          "cox_log_hr": (io.loc[g, "cox_log_hr"] - tki.loc[g, "cox_log_hr"]).to_numpy(),
                          "cox_se": np.sqrt(io.loc[g, "cox_se"] ** 2 + tki.loc[g, "cox_se"] ** 2).to_numpy(),
                          "cox_pvalue": 2 * stats.norm.sf(np.abs(z.to_numpy()))}).dropna(subset=["cox_pvalue"])
    inter["cox_padj"] = bh(inter["cox_pvalue"])
    for k in SIGNATURES:
        a, t = controls["atezo_bev"][k], controls["sunitinib"][k]
        d = np.log(a["hr_per_sd"]) - np.log(t["hr_per_sd"])
        zz = d / np.sqrt(a["se"] ** 2 + t["se"] ** 2)
        controls.setdefault("atezo_bev_vs_sunitinib", {})[k] = {"ratio_of_hr": round(float(np.exp(d)), 3),
                                                                 "pvalue": float(f"{2 * stats.norm.sf(abs(zz)):.3g}")}

    # nivolumab: direction check only
    c2 = _read_clinical("choueiri2016")
    x2 = _read_expression("choueiri2016", c2.index)
    niv = _arm_stats(x2, c2, "nivolumab", "Choueiri2016")
    df = pd.concat([arm, inter, niv], ignore_index=True)
    df.to_parquet(PROCESSED / "layer_treatment.parquet", index=False)

    ok_angio = controls["sunitinib"]["angiogenesis"]["hr_per_sd"] < 1 and controls["sunitinib"]["angiogenesis"]["pvalue"] < 0.05
    ok_teff = controls["atezo_bev_vs_sunitinib"]["T_effector"]["ratio_of_hr"] < 1
    summary = {
        "studies": {"IMmotion150": {"patients": int(len(clin)), "arms": clin["arm"].value_counts().to_dict(),
                                    "responders": int(clin["responder"].sum()), "genes_tested": int(len(x))},
                    "Choueiri2016": {"patients": int(len(c2)), "responders": int(c2["responder"].sum()),
                                     "genes_tested": int(len(x2))}},
        "signature_controls": controls,
        "controls_reproduced": {"angiogenesis_longer_pfs_on_sunitinib": bool(ok_angio),
                                "t_effector_favours_atezo_bev_over_sunitinib": bool(ok_teff)},
        "significant": {a: int((df[(df["arm"] == a)]["cox_padj"] < 0.1).sum()) for a in df["arm"].unique()},
        "top_predictive": inter.nsmallest(10, "cox_pvalue")[["gene", "cox_log_hr", "cox_padj"]].round(4).to_dict("records"),
    }
    (PROCESSED / "treatment_summary.json").write_text(json.dumps(summary, indent=1, default=float))
    return summary


if __name__ == "__main__":
    print(json.dumps(run(), indent=1, default=float))
