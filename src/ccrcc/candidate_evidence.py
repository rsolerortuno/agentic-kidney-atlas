"""Evidence coverage and modality rationale, without changing ranking weights.

Missing observations are never converted to negative biological evidence here.
These are transparent descriptive rules, not a clinical treatment recommendation.
"""
from __future__ import annotations

import math
from pathlib import Path

import pandas as pd

from .config import load
from .layers import selective_dependency
from .paths import PROCESSED


def finite(value) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def observation(name, status, detail, source, fields=None):
    return dict(name=name, status=status, detail=detail, source=source, fields=fields or {})


def expression(name, row, lfc, q, source, reference=None, min_effect=0, replication=False):
    effect, padj = row.get(lfc), row.get(q)
    if not finite(effect) or not finite(padj):
        return observation(name, "missing", "No usable effect and adjusted p-value", source)
    fields = {lfc: effect, q: padj}
    if replication and (not finite(reference) or reference == 0):
        return observation(name, "missing", "Replication comparison unavailable: discovery direction is missing or zero; replication-cohort measurements are shown", source, fields)
    status = "non_favorable"
    detail = "Measured; does not meet the stated effect/significance criterion"
    if padj < .05 and abs(effect) >= min_effect:
        status, detail = "favorable", "Significant expression difference (not proof of dependency)"
        if finite(reference) and effect * reference < 0:
            status, detail = "contradictory", "Significant change in the opposite direction to discovery"
    elif finite(reference) and effect * reference < 0:
        detail = "Opposite point-estimate direction, but not significant; replication not established"
    return observation(name, status, detail, source, fields)


def build(genes: list[str], processed: Path = PROCESSED) -> list[dict]:
    names = ["gene_universe", "layer_tcga", "confirmatory_protein", "sc_replication",
             "layer_dependency", "adc_suitability"]
    tables = {}
    for name in names:
        path = processed / (name + ".parquet")
        tables[name] = pd.read_parquet(path).drop_duplicates("gene").set_index("gene") if path.exists() else pd.DataFrame()
    def row(name, gene):
        table = tables[name]
        return table.loc[gene].to_dict() if gene in table.index else {}
    out = []
    thresholds = load()["omics"]["min_abs_log2fc"]
    for gene in dict.fromkeys(genes):
        u, tc, conf, sc, dep, adc = [row(name, gene) for name in names]
        obs = [expression("Discovery RNA", u, "rna_log2fc", "rna_padj", "gene_universe.parquet", min_effect=thresholds["rna"]),
               expression("Discovery protein", u, "protein_log2fc", "protein_padj", "gene_universe.parquet", min_effect=thresholds["protein"]),
               expression("TCGA RNA replication", tc, "log2fc", "padj", "layer_tcga.parquet", reference=u.get("rna_log2fc"), replication=True)]
        if not finite(conf.get("conf_n_pairs")) or conf["conf_n_pairs"] < 10:
            obs.append(observation("Confirmatory protein", "missing", "Not tested or fewer than 10 matched pairs", "confirmatory_protein.parquet", {"n_pairs": conf.get("conf_n_pairs")}))
        else:
            obs.append(expression("Confirmatory protein", conf, "conf_log2fc", "conf_padj", "confirmatory_protein.parquet", reference=u.get("protein_log2fc"), replication=True))
        if not sc or not isinstance(sc.get("li_top"), str) or not isinstance(sc.get("zhang_top"), str):
            obs.append(observation("Cell attribution", "missing", "Two usable cell-type assignments unavailable", "sc_replication.parquet"))
        else:
            agree = sc["li_top"] == sc["zhang_top"]
            obs.append(observation("Cell attribution", "favorable" if agree else "contradictory",
                "Dominant cell source agrees" if agree else "Dominant cell source differs across cohorts",
                "sc_replication.parquet", {k:sc.get(k) for k in ("li_top", "zhang_top", "li_rcc_share", "zhang_rcc_share")}))
        dep_ok = False
        if dep and finite(dep.get("mean_effect_ccrcc")) and finite(dep.get("padj")):
            dep_ok = bool(selective_dependency(pd.DataFrame([dep])).iloc[0])
            obs.append(observation("Selective cell-line dependency", "favorable" if dep_ok else "non_favorable",
                "Meets the configured selective-dependency rule" if dep_ok else "Measured; does not meet the selective-dependency rule",
                "layer_dependency.parquet", {k:dep.get(k) for k in ("mean_effect_ccrcc", "n_ccrcc_dependent", "n_ccrcc_lines", "padj")}))
        else:
            obs.append(observation("Selective cell-line dependency", "missing", "Usable dependency results unavailable", "layer_dependency.parquet"))
        if not adc or not adc.get("ihc_measured") or not finite(adc.get("ihc_vital_max_level")):
            obs.append(observation("Vital-organ staining", "missing", "Staining evidence unavailable; this is not evidence of safety", "adc_suitability.parquet"))
        else:
            level = adc["ihc_vital_max_level"]
            obs.append(observation("Vital-organ staining", "non_favorable" if level > 0 else "favorable",
                "Staining detected in a vital organ; potential off-tumour exposure" if level > 0 else "No staining in the assayed vital-organ panel; safety remains unproven",
                "adc_suitability.parquet", {"maximum_level":level,"tissue":adc.get("ihc_vital_max_tissue")}))
        up = finite(u.get("protein_log2fc")) and u["protein_log2fc"] > 0 and u.get("protein_sig",False)
        surface = finite(adc.get("c_surface")) and adc["c_surface"] > 0
        expression_signal = any(o["status"] == "favorable" for o in obs[:2])
        modes = [dict(name="Biomarker", status="exploratory" if expression_signal else "not_established", rationale="Expression differences motivate measurement studies; diagnostic, prognostic and predictive utility require separate validation." if expression_signal else "No discovery expression difference meeting the configured criteria is established here; biomarker utility is not established.")]
        if up and surface:
            modes.append(dict(name="Antibody / ADC antigen", status="hypothesis", rationale="Protein increase plus surface-topology annotation motivates accessibility studies. Staining and cell attribution above remain limiting evidence."))
        modes.append(dict(name="Inhibition", status="hypothesis" if dep_ok else "not_established",
            rationale="Selective dependency motivates functional intervention, but does not establish chemical inhibition." if dep_ok else "Expression alone does not establish benefit from inhibiting this protein; functional intervention evidence is insufficient here."))
        modes.append(dict(name="Degradation / MGD", status="not_established", rationale="No candidate-specific productive ternary complex or degradation experiment is established by this analysis. A pocket or exposed loop does not establish an MGD mechanism."))
        measured = sum(x["status"] != "missing" for x in obs)
        out.append(dict(gene=gene, observations=obs, modalities=modes, measured=measured, total=len(obs)))
    return out
