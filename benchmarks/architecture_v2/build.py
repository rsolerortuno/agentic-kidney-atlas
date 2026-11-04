"""Build the v2 benchmark: integrated, masked evidence dossiers with rubric-derived answers.

v1 (PR #10) gave each case two records and wrote the decision rule into every question, so a modern model, a single
agent or a team would all score near 100 % and the comparison could not discriminate. v2 tests what decomposition is
supposed to help with: many heterogeneous records per case, look-alike distractors from other entities, decisive
details hidden among irrelevant ones, and one integrated verdict that depends on several layers.

Every dossier is built from the processed atlas tables (read-only). Gene names are replaced by ENTITY codes. Half of
the dossiers are native; the other half are counterfactual copies of a native dossier with one decisive record
changed (recorded in provenance.json). Answers come from the rubric in PROTOCOL, which every arm also receives, so
the task is integration of supplied evidence under stated definitions, not guessing hidden criteria.
"""
from __future__ import annotations

import copy
import hashlib
import json
import random
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
PROCESSED = ROOT / "data" / "processed"
HERE = Path(__file__).parent
FROZEN = HERE / "frozen"

GENES = ["CA9", "CD70", "NNMT", "SLC16A3", "ENPP3", "HAVCR1", "COL23A1", "RAB42", "PSAT1", "PCK1", "CD74", "HK2"]
DISTRACTOR_POOL = ["ANGPTL4", "VWF", "LOX", "ALDOB", "UMOD", "PLIN2", "ESM1", "SLC6A19", "TIMP1", "FABP6", "KDR", "MT1G"]
QUESTIONS = ["protein_up", "replicated", "tumour_intrinsic", "surface_accessible", "vital_safety_concern",
             "selective_dependency"]
VERDICTS = ["surface_target_candidate", "loss_biomarker", "infiltration_signal", "needs_more_data"]

PROTOCOL = """Definitions (apply them to the records of the TARGET entity only; records of other entities are not evidence
about it). Each sub-question is answered yes, no or insufficient. Missing data are insufficient, never no.
1 protein_up: CPTAC discovery protein record. yes if padj<0.05 and log2fc>0; no if padj<0.05 and log2fc<0;
  insufficient if padj>=0.05 or no protein record.
2 replicated: reference direction = CPTAC discovery RNA, which must have padj<0.05 (otherwise insufficient).
  Independent cohorts: TCGA RNA and the confirmatory protein cohort. yes if at least one independent cohort has
  padj<0.05 in the reference direction and none has padj<0.05 in the opposite direction; no if any independent
  cohort has padj<0.05 opposite to the reference; insufficient otherwise.
3 tumour_intrinsic: single-cell tumour record (share_of_tumour_transcripts per cell group) and purity record
  (rna_lfc_at_purity1 with its SE). no if a non-RCC group has a larger transcript share than RCC; yes if RCC has the
  largest share AND rna_lfc_at_purity1 - 1.96*SE > 0; insufficient otherwise (missing record, or the purity interval
  includes or is below 0). pct_expressing is the share of cells with any expression and is NOT the transcript share.
4 surface_accessible: surface record. yes if surface_call = surface; no if internal or secreted; insufficient if
  unknown or no record.
5 vital_safety_concern: normal-tissue staining (IHC) and normal proteome records. yes if measured vital-organ staining
  level >= 2 or vital_protein_enriched is true; no if staining was measured with level <= 1 and the proteome does not
  report enrichment; insufficient if staining was not measured and the proteome does not report enrichment. Kidney
  staining is not a vital organ here.
6 selective_dependency: DepMap record. yes if n_ccrcc_lines>=10, padj<0.1, selectivity>0, n_ccrcc_dependent>=3 and
  not common_essential; no if a record exists and any condition fails; insufficient if no record.
verdict (first rule that applies):
  surface_target_candidate: protein_up=yes, replicated=yes, tumour_intrinsic=yes, surface_accessible=yes and
    vital_safety_concern=no;
  loss_biomarker: protein_up=no and replicated=yes;
  infiltration_signal: tumour_intrinsic=no and protein_up is not no;
  needs_more_data: otherwise.
Treatment, survival and within-arm associations are context: none of the definitions uses them."""


def _f(x, nd=4):
    return None if x is None or pd.isna(x) else (float(f"{x:.3g}") if abs(x) < 1e-3 and x != 0 else round(float(x), nd))


def tables() -> dict:
    t = {n: pd.read_parquet(PROCESSED / f"{n}.parquet") for n in
         ("de_results", "layer_tcga", "confirmatory_protein", "layer_purity", "sc_celltype_expression",
          "adc_suitability", "normal_proteome", "layer_dependency", "layer_treatment")}
    t["sc"] = t.pop("sc_celltype_expression").query("tissue == 'tumor'")
    return t


def records(gene: str, t: dict) -> list[dict]:
    """All dossier records for one gene (type + fields), without ids or entity codes."""
    out = []
    de = t["de_results"].query("hgnc_symbol == @gene")
    for r in de.itertuples():
        out.append({"type": "cptac_expression", "layer": r.layer, "log2fc": _f(r.log2fc), "padj": _f(r.padj),
                    "n_pairs": int(r.n_pairs)})
    tc = t["layer_tcga"].query("gene == @gene")
    for r in tc.itertuples():
        if pd.notna(r.log2fc):
            out.append({"type": "tcga_expression", "layer": "rna", "log2fc": _f(r.log2fc), "padj": _f(r.padj),
                        "n_pairs": int(r.n_pairs)})
        if pd.notna(r.cox_log_hr_per_sd):
            out.append({"type": "tcga_survival", "log_hr_per_sd": _f(r.cox_log_hr_per_sd), "padj": _f(r.cox_padj)})
    cf = t["confirmatory_protein"].query("gene == @gene and conf_n_pairs >= 10")
    for r in cf.itertuples():
        out.append({"type": "confirmatory_protein", "log2fc": _f(r.conf_log2fc), "padj": _f(r.conf_padj),
                    "n_pairs": int(r.conf_n_pairs)})
    pu = t["layer_purity"].query("hgnc_symbol == @gene")
    for r in pu.itertuples():
        if pd.notna(r.rna_lfc_at_purity1):
            out.append({"type": "purity_adjustment", "rna_lfc_at_purity1": _f(r.rna_lfc_at_purity1, 3),
                        "rna_lfc_at_purity1_se": _f(r.rna_lfc_at_purity1_se, 3), "rna_purity_rho": _f(r.rna_purity_rho, 3)})
    sc = t["sc"].query("hgnc_symbol == @gene").sort_values("share_of_tissue_transcripts", ascending=False)
    if len(sc):
        top = sc.head(5)
        if "RCC" not in set(top.cell_group):
            top = pd.concat([top, sc[sc.cell_group == "RCC"]])
        out.append({"type": "single_cell_tumour", "groups": [
            {"cell_group": r.cell_group, "share_of_tumour_transcripts": _f(r.share_of_tissue_transcripts, 3),
             "pct_expressing": _f(r.pct_expressing, 1)} for r in top.itertuples()]})
    ad = t["adc_suitability"].query("gene == @gene")
    for r in ad.itertuples():
        out.append({"type": "surface_annotation", "surface_call": r.surface_call, "evidence": r.surface_evidence,
                    "transmembrane_segments": None if pd.isna(r.n_transmembrane) else int(r.n_transmembrane)})
        out.append({"type": "normal_tissue_ihc", "measured": bool(pd.notna(r.ihc_vital_max_level)),
                    "vital_max_level": None if pd.isna(r.ihc_vital_max_level) else _f(r.ihc_vital_max_level, 1),
                    "vital_max_tissue": None if pd.isna(r.ihc_vital_max_tissue) else r.ihc_vital_max_tissue,
                    "kidney_max_level": None if pd.isna(r.ihc_kidney_max_level) else _f(r.ihc_kidney_max_level, 1)})
    npr = t["normal_proteome"].query("gene == @gene")
    for r in npr.itertuples():
        out.append({"type": "normal_proteome", "vital_tissues_detected": int(r.vital_protein_detected),
                    "vital_protein_ts_max": _f(r.vital_protein_ts_max, 2), "vital_protein_ts_tissue": r.vital_protein_ts_tissue,
                    "vital_protein_enriched": bool(r.vital_protein_enriched)})
    dp = t["layer_dependency"].query("gene == @gene")
    for r in dp.itertuples():
        out.append({"type": "depmap", "n_ccrcc_lines": int(r.n_ccrcc_lines), "n_ccrcc_dependent": int(r.n_ccrcc_dependent),
                    "selectivity": _f(r.selectivity, 3), "padj": _f(r.padj), "common_essential": bool(r.common_essential),
                    "mean_effect_ccrcc": _f(r.mean_effect_ccrcc, 3)})
    tr = t["layer_treatment"].query("gene == @gene and study in ['IMmotion150', 'CheckMate025']")
    for r in tr.itertuples():
        kind = "treatment_interaction" if "_vs_" in r.arm else "treatment_within_arm"
        out.append({"type": kind, "study": r.study, "arm": r.arm, "log_hr_per_sd": _f(r.cox_log_hr, 3), "padj": _f(r.cox_padj)})
    return out


# ------------------------------------------------------------------ rubric
def _one(recs, typ, **kw):
    for r in recs:
        if r["type"] == typ and all(r.get(k) == v for k, v in kw.items()):
            return r
    return None


def gold(recs: list[dict]) -> dict:
    sig = lambda r: r is not None and r.get("padj") is not None and r["padj"] < 0.05
    prot = _one(recs, "cptac_expression", layer="protein")
    q1 = "insufficient" if not sig(prot) else ("yes" if prot["log2fc"] > 0 else "no")
    rna = _one(recs, "cptac_expression", layer="rna")
    if not sig(rna):
        q2 = "insufficient"
    else:
        ref = 1 if rna["log2fc"] > 0 else -1
        ind = [r for r in (_one(recs, "tcga_expression"), _one(recs, "confirmatory_protein")) if sig(r)]
        signs = [1 if r["log2fc"] > 0 else -1 for r in ind]
        q2 = "no" if -ref in signs else "yes" if ref in signs else "insufficient"
    sc, pu = _one(recs, "single_cell_tumour"), _one(recs, "purity_adjustment")
    q3 = "insufficient"
    if sc:
        shares = {g["cell_group"]: g["share_of_tumour_transcripts"] for g in sc["groups"]
                  if g.get("share_of_tumour_transcripts") is not None}
        rcc = shares.get("RCC")
        if rcc is not None and any(v > rcc for k, v in shares.items() if k != "RCC"):
            q3 = "no"
        elif rcc is not None and pu and pu["rna_lfc_at_purity1"] is not None and pu["rna_lfc_at_purity1_se"] is not None \
                and pu["rna_lfc_at_purity1"] - 1.96 * pu["rna_lfc_at_purity1_se"] > 0:
            q3 = "yes"
    su = _one(recs, "surface_annotation")
    q4 = "insufficient" if su is None else {"surface": "yes", "internal": "no", "secreted": "no"}.get(su["surface_call"], "insufficient")
    ihc, npr = _one(recs, "normal_tissue_ihc"), _one(recs, "normal_proteome")
    enriched = bool(npr and npr["vital_protein_enriched"])
    if (ihc and ihc["measured"] and ihc["vital_max_level"] >= 2) or enriched:
        q5 = "yes"
    elif ihc and ihc["measured"]:
        q5 = "no"
    else:
        q5 = "insufficient"
    dp = _one(recs, "depmap")
    q6 = "insufficient" if dp is None else ("yes" if dp["n_ccrcc_lines"] >= 10 and dp["padj"] is not None and dp["padj"] < 0.1
                                           and dp["selectivity"] > 0 and dp["n_ccrcc_dependent"] >= 3
                                           and not dp["common_essential"] else "no")
    a = dict(zip(QUESTIONS, [q1, q2, q3, q4, q5, q6]))
    a["verdict"] = verdict(a)
    return a


def verdict(a: dict) -> str:
    if (a.get("protein_up"), a.get("replicated"), a.get("tumour_intrinsic"), a.get("surface_accessible"),
            a.get("vital_safety_concern")) == ("yes", "yes", "yes", "yes", "no"):
        return "surface_target_candidate"
    if a.get("protein_up") == "no" and a.get("replicated") == "yes":
        return "loss_biomarker"
    if a.get("tumour_intrinsic") == "no" and a.get("protein_up") != "no":
        return "infiltration_signal"
    return "needs_more_data"


# ------------------------------------------------------------------ counterfactuals
def _perturbations():
    """Each changes one decisive record of a native dossier; the new answer is recomputed by the rubric."""
    def flip_tcga(rs):
        r = _one(rs, "tcga_expression"); r["log2fc"] = -abs(r["log2fc"]) if r["log2fc"] > 0 else abs(r["log2fc"]); r["padj"] = 1e-08
        rs[:] = [x for x in rs if x["type"] != "confirmatory_protein"]
        return "TCGA RNA set significant in the opposite direction; confirmatory record removed"
    def rcc_loses_share(rs):
        g = _one(rs, "single_cell_tumour")["groups"]
        rcc = next(x for x in g if x["cell_group"] == "RCC"); other = next(x for x in g if x["cell_group"] != "RCC")
        rcc["share_of_tumour_transcripts"], other["share_of_tumour_transcripts"] = 0.2, 0.6
        rcc["pct_expressing"] = 70.0
        return "RCC transcript share lowered below another group while RCC pct_expressing raised to 70"
    def internal(rs):
        _one(rs, "surface_annotation")["surface_call"] = "internal"
        _one(rs, "surface_annotation")["evidence"] = "UniProt: Endoplasmic reticulum membrane | SURFY: not classified | CSPA: not observed"
        return "surface call changed to internal (ER membrane)"
    def stain(rs):
        r = _one(rs, "normal_tissue_ihc"); r.update(measured=True, vital_max_level=3.0, vital_max_tissue="liver")
        return "vital-organ staining set to level 3 in liver"
    def protein_gone(rs):
        rs[:] = [x for x in rs if not (x["type"] == "cptac_expression" and x["layer"] == "protein")]
        return "discovery protein record withheld"
    def purity_wide(rs):
        r = _one(rs, "purity_adjustment"); r["rna_lfc_at_purity1_se"] = round(abs(r["rna_lfc_at_purity1"]) / 1.2, 3)
        return "purity SE widened so the interval includes 0"
    def protein_down(rs):
        r = _one(rs, "cptac_expression", layer="protein"); r["log2fc"] = -abs(r["log2fc"]); r["padj"] = 1e-10
        rn = _one(rs, "cptac_expression", layer="rna"); rn["log2fc"] = -abs(rn["log2fc"]); rn["padj"] = 1e-10
        t = _one(rs, "tcga_expression"); t["log2fc"] = -abs(t["log2fc"]); t["padj"] = 1e-08
        rs[:] = [x for x in rs if x["type"] != "confirmatory_protein"]
        return "discovery RNA, protein and TCGA set significantly down; confirmatory removed"
    def clear_safety(rs):
        r = _one(rs, "normal_tissue_ihc"); r.update(measured=True, vital_max_level=1.0, vital_max_tissue="lung")
        p = _one(rs, "normal_proteome")
        if p:
            p.update(vital_protein_enriched=False, vital_protein_ts_max=1.1)
        return "vital-organ staining set to level 1 and proteome enrichment removed"
    def dependency(rs):
        r = _one(rs, "depmap"); r.update(padj=0.01, selectivity=0.4, n_ccrcc_dependent=8, common_essential=False)
        return "DepMap record set to a selective dependency"
    return [("flip_tcga", flip_tcga), ("rcc_loses_share", rcc_loses_share), ("internal", internal), ("stain", stain),
            ("protein_gone", protein_gone), ("purity_wide", purity_wide), ("protein_down", protein_down),
            ("dependency", dependency), ("clear_safety", clear_safety)]


def _applicable(name, rs):
    need = {"flip_tcga": ["tcga_expression"], "rcc_loses_share": ["single_cell_tumour"], "internal": ["surface_annotation"],
            "stain": ["normal_tissue_ihc"], "protein_gone": [], "purity_wide": ["purity_adjustment"],
            "protein_down": ["tcga_expression"], "dependency": ["depmap"], "clear_safety": ["normal_tissue_ihc"]}[name]
    return all(_one(rs, t) for t in need)


def build(out: Path = FROZEN, seed: int = 20261104) -> dict:
    if (out / "manifest.json").exists():
        raise FileExistsError("frozen suite exists; build into a new folder instead of editing it")
    t = tables()
    rng = random.Random(seed)
    codes = [f"ENTITY_{i:02d}" for i in range(1, 100)]
    rng.shuffle(codes)
    native = {g: records(g, t) for g in GENES}
    distract = {g: records(g, t) for g in DISTRACTOR_POOL}
    perts = _perturbations()
    cases, golds, prov = [], {}, {}
    seen = {v: 0 for v in VERDICTS}
    for g in GENES:
        seen[gold(native[g])["verdict"]] += 1
    for i, g in enumerate(GENES):
        variants = [("native", "none", copy.deepcopy(native[g]))]
        base_gold = gold(native[g])
        # first perturbation (gene-specific rotation) that changes the verdict; failing that, one that changes any answer
        options = []
        for k in range(len(perts)):
            name, fn = perts[(i + k) % len(perts)]
            if _applicable(name, native[g]):
                rs = copy.deepcopy(native[g]); note = fn(rs)
                if gold(rs) != base_gold:
                    options.append((gold(rs)["verdict"] != base_gold["verdict"], name, note, rs))
        if options:
            # balance the verdict classes: among verdict-changing options take the rarest resulting verdict
            changing = [o for o in options if o[0]] or options
            _, name, note, rs = min(changing, key=lambda o: seen[gold(o[3])["verdict"]])
            seen[gold(rs)["verdict"]] += 1
            variants.append((name, note, rs))
        for v, (kind, note, rs) in enumerate(variants):
            cid = f"g{i + 1:02d}v{v}"
            target = codes[2 * i + v]
            others = rng.sample([d for d in DISTRACTOR_POOL], 2)
            pool = [dict(r, entity=target) for r in rs]
            for j, o in enumerate(others):
                # look-alike records of other entities: same record types, different values, no "unrelated" label
                pick = rng.sample(distract[o], min(6, len(distract[o])))
                pool += [dict(r, entity=codes[60 + 2 * i + j]) for r in pick]
            rng.shuffle(pool)
            for n, r in enumerate(pool):
                r["id"] = f"{cid}.r{n + 1:02d}"
            cases.append({"id": cid, "target_entity": target, "records": pool})
            golds[cid] = {**gold(rs), "gene_index": i, "variant": v, "kind": kind, "is_synthetic": v > 0,
                          "target_record_ids": [r["id"] for r in pool if r["entity"] == target]}
            prov[cid] = {"gene": g, "perturbation": kind, "note": note, "distractor_genes": others}
    out.mkdir(parents=True, exist_ok=True)
    files = {"cases.json": cases, "gold.json": golds, "provenance.json": prov, "protocol.txt": PROTOCOL}
    for name, obj in files.items():
        text = obj if isinstance(obj, str) else json.dumps(obj, indent=1, ensure_ascii=False)
        (out / name).write_text(text + "\n", encoding="utf-8", newline="\n")
    manifest = {"version": 2, "seed": seed, "n_cases": len(cases), "n_genes": len(GENES),
                "files": {n: hashlib.sha256((out / n).read_bytes()).hexdigest() for n in files}}
    (out / "manifest.json").write_text(json.dumps(manifest, indent=1) + "\n", encoding="utf-8", newline="\n")
    return manifest


if __name__ == "__main__":
    print(json.dumps(build(), indent=1))
