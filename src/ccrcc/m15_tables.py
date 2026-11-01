"""Tables and summaries of the M15 analyses, shaped for the report template (no LLM, nothing computed here)."""
from __future__ import annotations

import json

import pandas as pd

from .paths import PROCESSED


def _json(name: str) -> dict:
    p = PROCESSED / name
    return json.loads(p.read_text(encoding="utf-8")) if p.exists() else {}


def _rows(df: pd.DataFrame, cols: list[str], n: int, nd: int = 3) -> list[dict]:
    out = []
    for r in df.head(n)[cols].to_dict("records"):
        out.append({k: (None if pd.isna(v) else (round(float(v), nd) if isinstance(v, float) else v)) for k, v in r.items()})
    return out


def build(candidates: list[str]) -> dict:
    out: dict = {"treatment": _json("treatment_summary.json"), "subtypes": _json("subtypes_summary.json"),
                 "adc_summary": _json("adc_summary.json"), "immune_summary": _json("immune_summary.json"),
                 "sc_rep": _json("sc_replication_summary.json"),
                 "confirmatory": _json("confirmatory_summary.json")}
    t = PROCESSED / "layer_treatment.parquet"
    if t.exists():
        d = pd.read_parquet(t)
        arms = {}
        for arm in ("atezo_bev", "sunitinib", "atezo_bev_vs_sunitinib"):
            a = d[d["arm"] == arm].dropna(subset=["cox_pvalue"]).sort_values("cox_pvalue")
            a["direction"] = a["cox_log_hr"].map(lambda x: "higher expression, longer benefit" if x < 0 else "higher expression, earlier progression")
            arms[arm] = _rows(a, ["gene", "cox_log_hr", "cox_padj", "direction"], 12)
        out["treatment_top"] = arms
        cand = d[(d["gene"].isin(candidates)) & d["arm"].isin(["atezo_bev", "sunitinib", "atezo_bev_vs_sunitinib"])]
        w = cand.pivot_table(index="gene", columns="arm", values=["cox_log_hr", "cox_padj"])
        w.columns = [f"{b}_{a}" for a, b in w.columns]
        out["treatment_candidates"] = _rows(w.reset_index(), ["gene"] + [c for c in w.columns], 40)
    s = PROCESSED / "subtype_targets.parquet"
    if s.exists():
        st = pd.read_parquet(s)
        out["subgroup_targets"] = {
            drv: _rows(g.sort_values(["tumour_up", "predicted_tractable", "score"], ascending=False),
                       ["gene", "cptac_rna_diff", "cptac_protein_diff", "tcga_rna_diff", "tumour_up", "predicted_tractable",
                        "clinical_drug"], 12)
            for drv, g in st.groupby("driver")}
    a = PROCESSED / "adc_suitability.parquet"
    if a.exists():
        d = pd.read_parquet(a).sort_values("rank")
        cols = ["rank", "gene", "suitability", "c_surface", "n_transmembrane", "extracellular_residues", "protein_log2fc",
                "ihc_vital_max_level", "ihc_vital_max_tissue", "ihc_kidney_max_level", "ihc_measured", "rcc_pct_donor_p25",
                "rcc_pct_donor_median", "control"]
        out["adc_top"] = _rows(d[d["c_surface"] > 0], cols, 25)
        out["adc_controls"] = _rows(d[d["control"]], cols, 10)
    i = PROCESSED / "immune_profile.parquet"
    if i.exists():
        d = pd.read_parquet(i).sort_values("rank")
        cols = ["rank", "gene", "immune_target_score", "immune_share", "c_coupling", "c_io_response", "rna_log2fc", "control"]
        out["immune_top"] = _rows(d, cols, 20)
        out["immune_controls"] = _rows(d[d["control"]], cols, 10)
    r = PROCESSED / "sc_replication.parquet"
    if r.exists():
        d = pd.read_parquet(r)
        out["sc_rep_targets"] = _rows(d[d["gene"].isin(candidates)].sort_values("gene"),
                                      ["gene", "li_top", "zhang_top", "li_rcc_share", "zhang_rcc_share", "li_tvp_log2fc",
                                       "zhang_tvp_log2fc"], 40)
    return out
