"""Load the treatment, subtype, ADC, immune, cell-replication and confirmatory-cohort results into the knowledge graph.

Every record goes through an Evidence record and is attached to its gene by a claim, so agents can cite it and
verify.py can check the gene link, the numbers and the freshness of every sentence that uses it.

Predicates (subject is a Gene node):
  response_assessed_in        -> dataset   IMmotion150 (and Choueiri 2016): progression and response per arm
  subtype_assessed_in         -> dataset   mutant vs wild type (VHL, PBRM1, BAP1, SETD2) and ccB vs ccA
  adc_assessed_in             -> dataset   surface topology, tumour increase, vital-organ staining, homogeneity
  immune_profile_assessed_in  -> dataset   link to CD8 / Treg / M2 / immune score, response, immune profile score
  cell_attribution_checked_in -> dataset   dominant cell source in a second single-cell dataset (Zhang 2021)
  confirmatory_replication_in -> dataset   protein and phosphosite change in the confirmatory CPTAC cohort

Only a focus set of genes is loaded (candidates of every ranking plus the discovery pool), which keeps the graph small
and covers every gene a report sentence is likely to name.
"""
from __future__ import annotations

import json
from datetime import date

import numpy as np
import pandas as pd

from .config import load
from .graph import GraphStore, evidence_id
from .paths import PROCESSED, RAW

TODAY = date.today().isoformat()
DS = {"treatment": ("DATASET:IMMOTION150", "IMmotion150 and Choueiri 2016 treatment-response trials (cBioPortal datahub)"),
      "subtype": ("DATASET:SUBTYPES", "CPTAC and TCGA-KIRC tumours split by driver mutation and ccA/ccB"),
      "adc": ("DATASET:ADC-SUITABILITY", "UniProt topology, Human Protein Atlas staining and scRNA homogeneity"),
      "immune": ("DATASET:IMMUNE-DECONVOLUTION", "CPTAC CIBERSORT and xCell immune deconvolution"),
      "zhang": ("DATASET:ZHANG2021-SCRNA", "Zhang et al. 2021 ccRCC single-cell RNA-seq (CZ CELLxGENE)"),
      "confirmatory": ("DATASET:CPTAC-CONFIRMATORY", "CPTAC ccRCC confirmatory cohort, DIA proteome and phosphoproteome (PDC000411)")}


def _r(x, nd=3):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return None
    x = float(x)
    return float(f"{x:.3g}") if abs(x) < 1e-3 and x != 0 else round(x, nd)


def _read(name):
    p = PROCESSED / name
    return pd.read_parquet(p) if p.exists() else None


def focus_genes(n: int = 300) -> set[str]:
    genes = set(load()["cell_context"]["targets"])
    s = _read("scores_v3.parquet")
    if s is not None:
        for col in ("score", "score_target", "score_biomarker"):
            if col in s:
                genes |= set(s.nlargest(n, col)["gene"])
    for name, col, key in (("discovery_candidates.parquet", None, "gene"), ("adc_suitability.parquet", "suitability", "gene"),
                           ("immune_profile.parquet", "immune_target_score", "gene")):
        d = _read(name)
        if d is not None:
            genes |= set(d[key] if col is None else d.nlargest(n, col)[key])
    t = _read("layer_treatment.parquet")
    if t is not None:
        for arm in ("atezo_bev", "sunitinib", "atezo_bev_vs_sunitinib"):
            genes |= set(t[t["arm"] == arm].nsmallest(60, "cox_pvalue")["gene"])
    st = _read("subtype_targets.parquet")
    if st is not None:
        genes |= set(st["gene"])
    return genes


def load_applied(store: GraphStore | None = None) -> dict:
    store = store or GraphStore()
    hg = pd.read_csv(RAW / "hgnc_complete_set.txt", sep="\t", usecols=["hgnc_id", "symbol"], dtype=str)
    gid = dict(zip(hg["symbol"], hg["hgnc_id"]))
    existing = {r[0] for r in store.conn.execute("SELECT id FROM nodes WHERE type = 'Gene'")}
    counts: dict[str, int] = {}
    genes = focus_genes()

    def claim(sym, pred, dsk, ev_type, props):
        h = gid.get(sym)
        if not h:
            return
        if h not in existing:
            store.add_node(h, "Gene", sym, hgnc_id=h)
            existing.add(h)
        ev = store.add_evidence(evidence_id(ev_type, sym), ev_type, DS[dsk][1], TODAY, gene=sym,
                                **{k: v for k, v in props.items() if v is not None})
        store.add_claim(h, pred, DS[dsk][0], ev)
        counts[ev_type] = counts.get(ev_type, 0) + 1

    with store.batch():
        for _, (did, name) in DS.items():
            store.add_node(did, "Dataset", name)

        t = _read("layer_treatment.parquet")
        if t is not None:
            by = {k: g.set_index("gene") for k, g in t.groupby("arm")}
            sig = json.loads((PROCESSED / "treatment_summary.json").read_text())
            for g in genes:
                rec, any_arm = {}, False
                for arm, label in (("atezo_bev", "atezolizumab_bevacizumab"), ("sunitinib", "sunitinib"), ("atezo", "atezolizumab"),
                                   ("atezo_bev_vs_sunitinib", "difference_between_arms"), ("nivolumab", "nivolumab_choueiri2016"),
                                   ("nivolumab_cm025", "nivolumab_checkmate025"), ("everolimus_cm025", "everolimus_checkmate025"),
                                   ("nivolumab_vs_everolimus", "difference_nivolumab_vs_everolimus")):
                    d = by.get(arm)
                    if d is not None and g in d.index:
                        r = d.loc[g]
                        rec[label] = {"log_hr_per_sd": _r(r["cox_log_hr"]), "padj": _r(r["cox_padj"]), "n": int(r["n"]),
                                      "direction": "longer benefit with higher expression" if r["cox_log_hr"] < 0 else
                                                   "earlier progression with higher expression"}
                        if "auc_response" in r and pd.notna(r.get("auc_response")):
                            rec[label]["response_auc"] = _r(r["auc_response"])
                        any_arm = True
                if any_arm:
                    claim(g, "response_assessed_in", "treatment", "treatment_response",
                          {"arms": rec, "tested_genes": sig["studies"]["IMmotion150"]["genes_tested"],
                           "caution": "no gene differs between arms after correction" if sig["significant"].get(
                               "atezo_bev_vs_sunitinib", 0) == 0 else None,
                           "caution_checkmate025": "no gene differs between nivolumab and everolimus after correction" if
                               sig["significant"].get("nivolumab_vs_everolimus", 0) == 0 else None})

        gen = _read("subtype_genetic_de.parquet")
        cc = _read("subtype_ccab_de.parquet")
        if gen is not None and cc is not None:
            gen_i = {k: g.set_index("gene") for k, g in gen.groupby(["driver", "layer"])}
            cc_i = {k: g.set_index("gene") for k, g in cc.groupby("layer")}
            st = _read("subtype_targets.parquet")
            sub_t = {(r.gene, r.driver) for r in st.itertuples()} if st is not None else set()
            for g in genes:
                rec = {}
                for drv in ("VHL", "PBRM1", "BAP1", "SETD2"):
                    item = {}
                    for lay in ("cptac_rna", "cptac_protein", "tcga_rna"):
                        d = gen_i.get((drv, lay))
                        if d is not None and g in d.index:
                            item[lay] = {"mutant_minus_wildtype_log2": _r(d.loc[g, "diff"]), "padj": _r(d.loc[g, "padj"])}
                    if item:
                        item["subgroup_target_candidate"] = (g, drv) in sub_t
                        rec[drv] = item
                ccab = {}
                for lay, d in cc_i.items():
                    if g in d.index:
                        ccab[lay] = {"ccB_minus_ccA_log2": _r(d.loc[g, "diff"]), "padj": _r(d.loc[g, "padj"])}
                if ccab:
                    ccab["defines_the_split"] = bool(cc.loc[cc["gene"] == g, "clearcode34_gene"].iloc[0])
                    rec["ccA_ccB"] = ccab
                if rec:
                    claim(g, "subtype_assessed_in", "subtype", "subtype_profile", rec)

        adc = _read("adc_suitability.parquet")
        if adc is not None:
            for r in adc[adc["gene"].isin(genes) | adc["control"]].itertuples():
                claim(r.gene, "adc_assessed_in", "adc", "adc_profile", {
                    "suitability": _r(r.suitability), "rank": int(r.rank), "of": int(len(adc)),
                    "surface_score": _r(r.c_surface), "transmembrane_segments": None if pd.isna(r.n_transmembrane) else int(r.n_transmembrane),
                    "extracellular_residues": None if pd.isna(r.extracellular_residues) else int(r.extracellular_residues),
                    "protein_log2fc": _r(r.protein_log2fc), "rna_log2fc": _r(r.rna_log2fc),
                    "vital_organ_staining_level": _r(r.ihc_vital_max_level), "vital_organ_tissue": None if pd.isna(r.ihc_vital_max_tissue) else r.ihc_vital_max_tissue,
                    "kidney_staining_level": _r(r.ihc_kidney_max_level), "staining_data_available": bool(r.ihc_measured),
                    "tumour_cells_expressing_pct_weakest_quarter_of_donors": _r(r.rcc_pct_donor_p25, 1),
                    "tumour_cells_expressing_pct_median_donor": _r(r.rcc_pct_donor_median, 1), "known_adc_or_radioligand_target": bool(r.control),
                    "surface_call": getattr(r, "surface_call", None), "surface_evidence": getattr(r, "surface_evidence", None),
                    "observed_at_cell_surface_cspa": bool(getattr(r, "cspa_high_confidence", False)),
                    "vital_organ_protein_ts_max": _r(getattr(r, "vital_protein_ts_max", np.nan)),
                    "vital_organ_protein_ts_tissue": None if pd.isna(getattr(r, "vital_protein_ts_tissue", None)) else r.vital_protein_ts_tissue,
                    "vital_organ_protein_enriched": None if pd.isna(getattr(r, "vital_protein_enriched", None)) else bool(r.vital_protein_enriched)})

        imm = _read("immune_profile.parquet")
        if imm is not None:
            for r in imm[imm["gene"].isin(genes) | imm["control"]].itertuples():
                claim(r.gene, "immune_profile_assessed_in", "immune", "immune_profile", {
                    "immune_target_score": _r(r.immune_target_score), "rank": int(r.rank), "of": int(len(imm)),
                    "immune_cell_share_of_tumour_transcripts": _r(r.immune_share),
                    "rho_with_cd8_t_cells": _r(r.cd8_rho), "rho_with_regulatory_t_cells": _r(r.treg_rho),
                    "rho_with_m2_macrophages": _r(r.m2_rho), "rho_with_immune_score": _r(r.immune_score_rho),
                    "log_hr_atezolizumab_bevacizumab": _r(getattr(r, "io_log_hr", None)), "rna_log2fc": _r(r.rna_log2fc),
                    "checkpoint_control": bool(r.control)})

        sc = _read("sc_replication.parquet")
        if sc is not None:
            for r in sc[sc["gene"].isin(genes)].itertuples():
                claim(r.gene, "cell_attribution_checked_in", "zhang", "sc_replication", {
                    "li2022_dominant_group": r.li_top, "zhang2021_dominant_group": r.zhang_top,
                    "same_dominant_group": bool(r.agree_top), "li2022_tumour_cell_share": _r(r.li_rcc_share),
                    "zhang2021_tumour_cell_share": _r(r.zhang_rcc_share), "share_profile_spearman": _r(r.share_profile_rho),
                    "li2022_tumour_vs_tubule_log2fc": _r(r.li_tvp_log2fc), "zhang2021_tumour_vs_tubule_log2fc": _r(r.zhang_tvp_log2fc),
                    "zhang2021_paired_donors": None if pd.isna(r.zhang_tvp_donors) else int(r.zhang_tvp_donors)})

        cp, cph = _read("confirmatory_protein.parquet"), _read("confirmatory_phospho.parquet")
        if cp is not None:
            ph = cph.groupby("gene") if cph is not None else None
            for r in cp[cp["gene"].isin(genes)].itertuples():
                props = {"protein_conf_log2fc": _r(r.conf_log2fc), "protein_conf_padj": _r(r.conf_padj),
                         "protein_conf_pairs": int(r.conf_n_pairs), "protein_discovery_log2fc": _r(r.disc_log2fc),
                         "protein_discovery_padj": _r(r.disc_padj), "protein_replicates": bool(r.replicates)}
                if ph is not None and r.gene in ph.groups:
                    g = ph.get_group(r.gene)
                    sg = g[g["disc_sig"]]
                    props.update({"phosphosites_measured_in_both": int(len(g)), "phosphosites_significant_in_discovery": int(len(sg)),
                                  "phosphosites_replicating": int(sg["replicates"].sum())})
                claim(r.gene, "confirmatory_replication_in", "confirmatory", "confirmatory_replication", props)
    return {"claims_added": counts, "focus_genes": len(genes)}


if __name__ == "__main__":
    print(json.dumps(load_applied(), indent=1))
