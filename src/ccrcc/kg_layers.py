"""M13: load the extra evidence layers into the knowledge graph (every claim through an Evidence record).

Predicates added (subject is always a Gene node):
  somatically_mutated_in   -> disease   recurrence above the length-corrected background (CPTAC + TCGA-KIRC)
  mutation_assessed_in     -> dataset   tested, not recurrent
  copy_number_drives       -> disease   copy number correlates with the gene's protein (or RNA) in tumours
  selective_dependency_in  -> disease   DepMap: ccRCC lines depend on the gene more than other lines
  dependency_assessed_in   -> dataset   measured, not selective
  upregulated_in / downregulated_in / not_changed_in -> DATASET:TCGA-KIRC-RNA   external replication
  prognostic_in            -> disease   TCGA-KIRC overall survival (Cox, BH < 0.05)
  activated_in / inhibited_in -> disease   kinase activity (protein-adjusted KSEA, BH < 0.1)
  phosphosite_changed_in   -> disease   site-level tumour vs normal (BH < 0.05, |log2FC| >= 1), one record per site
  purity_assessed_in       -> dataset   correlation with ESTIMATE purity and log2FC extrapolated to purity 1
"""
from __future__ import annotations

import json
from datetime import date

import numpy as np
import pandas as pd

from .config import load
from .graph import GraphStore, evidence_id
from .paths import PROCESSED, RAW

DISEASE = "MONDO:0005005"
TODAY = date.today().isoformat()
DATASETS = {
    "mutations": ("DATASET:SOMATIC-CPTAC-TCGA", "Somatic mutations, CPTAC ccRCC discovery + TCGA-KIRC (cBioPortal)"),
    "depmap": ("DATASET:DEPMAP-26Q1", "DepMap Public 26Q1 CRISPR (Chronos) gene effect"),
    "tcga": ("DATASET:TCGA-KIRC-RNA", "TCGA-KIRC RNA-seq tumour vs adjacent normal + overall survival (UCSC Xena)"),
    "phosphosite": ("DATASET:CPTAC-CCRCC-PHOSPHOSITE", "CPTAC ccRCC phosphosites (BCM pipeline)"),
    "purity": ("DATASET:CPTAC-CCRCC-PURITY", "CPTAC ESTIMATE tumour purity (WashU)"),
    "cnv": ("DATASET:CPTAC-CCRCC-CNV", "CPTAC ccRCC discovery copy number (gene level)"),
}


def _f(x, nd=4):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return None
    return float(f"{x:.3g}") if abs(x) < 1e-3 and x != 0 else round(float(x), nd)


def load_layers(store: GraphStore | None = None) -> dict:
    store = store or GraphStore()
    cfg = load()["convergence"]
    hg = pd.read_csv(RAW / "hgnc_complete_set.txt", sep="\t", usecols=["hgnc_id", "symbol"], dtype=str)
    gid = dict(zip(hg["symbol"], hg["hgnc_id"]))
    existing = {r[0] for r in store.conn.execute("SELECT id FROM nodes WHERE type = 'Gene'")}
    counts: dict[str, int] = {}

    def gene(sym):
        h = gid.get(sym)
        if not h:
            return None
        if h not in existing:
            store.add_node(h, "Gene", sym, hgnc_id=h)
            existing.add(h)
        return h

    def claim(sym, pred, obj, ev_type, ev_key, source, **props):
        g = gene(sym)
        if g is None:
            return
        ev = store.add_evidence(evidence_id(ev_type, *ev_key), ev_type, source, TODAY,
                                **{k: v for k, v in props.items() if v is not None})
        store.add_claim(g, pred, obj, ev)
        counts[ev_type] = counts.get(ev_type, 0) + 1

    with store.batch():
        for _, (did, name) in DATASETS.items():
            store.add_node(did, "Dataset", name)

        p = PROCESSED / "layer_mutations.parquet"
        if p.exists():
            meta = json.loads((PROCESSED / "layer_mutations_meta.json").read_text())
            for r in pd.read_parquet(p).itertuples():
                if r.mutated < 2:
                    continue
                sig = r.padj < 0.05
                claim(r.gene, "somatically_mutated_in" if sig else "mutation_assessed_in",
                      DISEASE if sig else DATASETS["mutations"][0], "somatic_recurrence", (r.gene,),
                      DATASETS["mutations"][0], gene=r.gene, mutated_patients=int(r.mutated), patients=int(r.patients),
                      freq=_f(r.freq), cptac_mutated=int(r.cptac_mutated), tcga_mutated=int(r.tcga_mutated),
                      expected=_f(r.expected), obs_exp=_f(r.obs_exp), padj=_f(r.padj), cds_bp=int(r.cds_bp),
                      method="binomial vs coding-length background, nonsynonymous, hypermutators excluded",
                      cohorts={k: v["patients"] for k, v in meta.items()})

        p = PROCESSED / "layer_cis_dosage.parquet"
        if p.exists():
            cnvf = pd.read_parquet(PROCESSED / "cnv_frequencies.parquet").dropna(subset=["hgnc_symbol"])
            cnvf = cnvf.drop_duplicates("hgnc_symbol").set_index("hgnc_symbol")
            for r in pd.read_parquet(p).itertuples():
                prot_sig = pd.notna(r.protein_padj) and r.protein_padj < 0.05 and r.protein_rho > 0
                rna_sig = pd.notna(r.rna_padj) and r.rna_padj < 0.05 and r.rna_rho > 0
                if not (prot_sig or rna_sig) or r.gene not in cnvf.index:
                    continue
                c = cnvf.loc[r.gene]
                claim(r.gene, "copy_number_drives", DISEASE, "cis_dosage", (r.gene,), DATASETS["cnv"][0], gene=r.gene,
                      rna_rho=_f(r.rna_rho, 3), rna_padj=_f(r.rna_padj), protein_rho=_f(r.protein_rho, 3),
                      protein_padj=_f(r.protein_padj), loss_freq=_f(c.loss_freq, 3), gain_freq=_f(c.gain_freq, 3),
                      arm=c.arm)

        p = PROCESSED / "layer_dependency.parquet"
        if p.exists():
            dc = cfg["dependency"]
            lines = json.loads((PROCESSED / "layer_dependency_lines.json").read_text())
            for r in pd.read_parquet(p).itertuples():
                if r.n_ccrcc_lines < dc["min_lines"] or not np.isfinite(r.pvalue):
                    continue
                sel = r.padj < dc["padj"] and r.selectivity > 0 and r.frac_ccrcc_dependent > 0 and not r.common_essential
                if not sel and r.frac_ccrcc_dependent < 0.25 and not r.common_essential:
                    continue   # keep the graph small: non-dependencies are not stored
                claim(r.gene, "selective_dependency_in" if sel else "dependency_assessed_in",
                      DISEASE if sel else DATASETS["depmap"][0], "depmap_dependency", (r.gene,), DATASETS["depmap"][0],
                      gene=r.gene, mean_effect_ccrcc=_f(r.mean_effect_ccrcc, 3), mean_effect_other=_f(r.mean_effect_other, 3),
                      selectivity=_f(r.selectivity, 3), frac_ccrcc_dependent=_f(r.frac_ccrcc_dependent, 3),
                      frac_all_dependent=_f(r.frac_all_dependent, 3), n_ccrcc_lines=int(r.n_ccrcc_lines),
                      padj=_f(r.padj), common_essential=bool(r.common_essential), release=lines["release"])

        p = PROCESSED / "layer_tcga.parquet"
        if p.exists():
            for r in pd.read_parquet(p).itertuples():
                if pd.notna(r.log2fc):
                    sig = r.padj < 0.05 and abs(r.log2fc) >= 1
                    pred = ("upregulated_in" if r.log2fc > 0 else "downregulated_in") if sig else "not_changed_in"
                    claim(r.gene, pred, DATASETS["tcga"][0], "tcga_de", (r.gene,), DATASETS["tcga"][0], gene=r.gene,
                          log2fc=_f(r.log2fc), padj=_f(r.padj), n_pairs=int(r.n_pairs), significant=bool(sig),
                          unit="log2(RSEM+1)", method="paired t-test, BH")
                if pd.notna(r.cox_padj) and r.cox_padj < 0.05:
                    claim(r.gene, "prognostic_in", DISEASE, "tcga_survival", (r.gene,), DATASETS["tcga"][0], gene=r.gene,
                          log_hr_per_sd=_f(r.cox_log_hr_per_sd, 3), hr_per_sd=_f(float(np.exp(r.cox_log_hr_per_sd)), 3),
                          padj=_f(r.cox_padj), n=int(r.n_survival), endpoint="overall survival",
                          direction="worse with higher expression" if r.cox_log_hr_per_sd > 0 else
                          "better with higher expression")

        p = PROCESSED / "layer_kinase_activity.parquet"
        if p.exists():
            for r in pd.read_parquet(p).itertuples():
                sig = r.padj < 0.1
                pred = ("activated_in" if r.z > 0 else "inhibited_in") if sig else "kinase_activity_assessed_in"
                claim(r.gene, pred, DISEASE if sig else DATASETS["phosphosite"][0], "kinase_activity", (r.gene,),
                      DATASETS["phosphosite"][0], gene=r.gene, ksea_z=_f(r.z, 3), padj=_f(r.padj),
                      n_substrate_sites=int(r.n_substrate_sites), mean_substrate_log2fc=_f(r.mean_substrate_log2fc, 3),
                      method="KSEA on protein-adjusted site log2FC, OmniPath kinase-substrate")

        p = PROCESSED / "layer_phosphosite.parquet"
        if p.exists():
            s = pd.read_parquet(p)
            s = s[(s["padj"] < 0.05) & (s["log2fc"].abs() >= 1)].sort_values("padj").drop_duplicates("site_id")
            for r in s.itertuples():
                claim(r.gene, "phosphosite_changed_in", DISEASE, "phosphosite_de", (r.site_id,),
                      DATASETS["phosphosite"][0], gene=r.gene, site=r.site, log2fc=_f(r.log2fc),
                      protein_adjusted_log2fc=_f(getattr(r, "protein_adjusted_log2fc", np.nan)), padj=_f(r.padj),
                      n_pairs=int(r.n_pairs))

        p = PROCESSED / "layer_purity.parquet"
        if p.exists():
            for r in pd.read_parquet(p).dropna(subset=["hgnc_symbol"]).drop_duplicates("hgnc_symbol").itertuples():
                if pd.isna(r.rna_lfc_at_purity1) and pd.isna(r.protein_lfc_at_purity1):
                    continue
                claim(r.hgnc_symbol, "purity_assessed_in", DATASETS["purity"][0], "purity", (r.hgnc_symbol,),
                      DATASETS["purity"][0], gene=r.hgnc_symbol, rna_purity_rho=_f(r.rna_purity_rho, 3),
                      rna_purity_padj=_f(r.rna_purity_padj), rna_lfc_at_purity1=_f(r.rna_lfc_at_purity1, 3),
                      rna_lfc_at_purity1_se=_f(r.rna_lfc_at_purity1_se, 3),
                      protein_purity_rho=_f(r.protein_purity_rho, 3),
                      protein_lfc_at_purity1=_f(r.protein_lfc_at_purity1, 3),
                      protein_lfc_at_purity1_se=_f(r.protein_lfc_at_purity1_se, 3),
                      n_tumors=int(r.n_tumors_with_purity), purity_source="ESTIMATE (RNA-based)")
    stats = store.stats()
    (PROCESSED / "graph_stats.json").write_text(json.dumps(stats, indent=1))
    return {"claims_added": counts, "claims_total": stats["claims_total"],
            "claims_without_resolvable_evidence": stats["claims_without_resolvable_evidence"]}


if __name__ == "__main__":
    print(json.dumps(load_layers(), indent=1))
