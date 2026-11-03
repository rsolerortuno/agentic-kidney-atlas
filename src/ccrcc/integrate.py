"""Deterministic convergence scoring (no LLM). Writes runs/<date>/convergence.json.

Every score component carries the evidence ids that justify it; only ids present in the graph are kept.
"""
from __future__ import annotations

import json
import os
from datetime import date
from pathlib import Path

import pandas as pd

from .config import load
from .graph import GraphStore, evidence_id
from .literature_entities import resolve_gene
from .paths import LITERATURE, PROCESSED, RUNS
from .scores import components, profile_scores


CURRENT = RUNS / ".current"


def current_run():
    """The run the agent step is working on: CCRCC_RUN_DIR, else the folder scripts/integrate.py started (runs/.current),
    else the latest dated run. Steps repeated on a later day must not silently switch to a new, empty folder."""
    if os.environ.get("CCRCC_RUN_DIR"):
        return run_dir()
    if CURRENT.exists():
        d = RUNS / CURRENT.read_text(encoding="utf-8").strip()
        if d.is_dir():
            return d
    runs = sorted(p for p in RUNS.glob("[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]") if p.is_dir())
    return runs[-1] if runs else run_dir()


def run_dir(day: str | None = None):
    """runs/<date>/ by default; CCRCC_RUN_DIR (relative to runs/ or absolute) redirects a whole /informe run, which
    is how the replicate runs used to measure run-to-run variability are kept apart (src/ccrcc/baselines.py)."""
    override = os.environ.get("CCRCC_RUN_DIR")
    if override and not day:
        d = Path(override) if Path(override).is_absolute() else RUNS / override
        d.mkdir(parents=True, exist_ok=True)
        return d
    d = RUNS / (day or date.today().isoformat())
    d.mkdir(parents=True, exist_ok=True)
    return d


def compute() -> pd.DataFrame:
    cfg = load()["convergence"]
    lw, sw = cfg["omics_layer_weights"], cfg["stream_weights"]
    t = pd.read_parquet(PROCESSED / "gene_universe.parquet")
    ot = pd.read_parquet(PROCESSED / "external_ot_associations.parquet").set_index("symbol")
    details_path = PROCESSED / "external_targets.json"
    details = json.loads(details_path.read_text()) if details_path.exists() else {}
    rel = pd.read_parquet(LITERATURE / "relations.parquet")

    lit_ids: dict[str, list[str]] = {}
    for r in rel.itertuples():
        for txt in (r.subject, r.object):
            g = resolve_gene(txt)
            if g:
                lit_ids.setdefault(g, []).append(evidence_id("literature", r.relation_id))

    t["n_literature"] = t["gene"].map(lambda g: len(set(lit_ids.get(g, []))))
    t["ot_literature"] = t["gene"].map(ot["dt_literature"]).fillna(0)
    t["c_literature"] = pd.concat([(t["n_literature"] / 3).clip(upper=1), t["ot_literature"]], axis=1).max(axis=1)
    t["ot_score"] = t["gene"].map(ot["score"]).fillna(0)
    t["ot_genetic"] = t["gene"].map(ot["dt_genetic_association"]).fillna(0)
    t["ot_somatic"] = t["gene"].map(ot["dt_somatic_mutation"]).fillna(0)
    if cfg.get("version", 2) >= 3:
        t = profile_scores(components(t))
    else:  # v2, kept so evals/history.md can be reproduced
        omics_raw = (lw["rna"] * t["rna_sig"] + lw["protein"] * t["protein_sig"] + lw["phospho"] * t["phospho_sig"]
                     + lw["concordance"] * t["concordant"])
        t["c_omics"] = omics_raw / sum(lw.values())
        t["c_cnv"] = t["cnv_freq"].fillna(0)
        t["c_genetics"] = t[["ot_genetic", "ot_somatic"]].max(axis=1)

        def tract(g):
            d = (details.get(g) or {}).get("open_targets")
            if not d:
                return 0.0
            return 1.0 if d["has_approved_or_clinical_drug"] else (0.5 if d["tractable_modalities"] else 0.0)

        t["c_tractability"] = t["gene"].map(tract)
    t["score"] = sum(sw[k] * t[f"c_{k}"] for k in sw)
    t["n_streams"] = sum((t[f"c_{k}"] > 0).astype(int) for k in sw)
    t["_lit_ids"] = t["gene"].map(lambda g: sorted(set(lit_ids.get(g, []))))
    return t.sort_values(["score", "n_streams"], ascending=False).reset_index(drop=True)


def _evidence(row, ot: pd.DataFrame, details: dict, store: GraphStore) -> dict:
    g = row["gene"]
    ot_ev = [evidence_id("open_targets_assoc", ot.loc[g, "ensembl"])] if g in ot.index else []
    ev = {"omics": [evidence_id("omics_de", layer, g) for layer in ("rna", "protein", "phospho")
                    if pd.notna(row.get(f"{layer}_padj"))],
          "cnv": [evidence_id("omics_cnv", g)],
          "literature": row["_lit_ids"] + (ot_ev if row["ot_literature"] > 0 else []),
          "genetics": (ot_ev if row.get("ot_genetic", 0) > 0 else [])
                      + ([evidence_id("somatic_recurrence", g)] if row.get("c_mutation", 0) > 0 else []),
          "tractability": [],
          "dependency": [evidence_id("depmap_dependency", g)] if row.get("c_dependency", 0) > 0 else [],
          "replication": [evidence_id("tcga_de", g)] if row.get("c_replication", 0) > 0 else [],
          "prognosis": [evidence_id("tcga_survival", g)] if row.get("c_prognosis", 0) > 0 else [],
          "kinase": [evidence_id("kinase_activity", g)] if row.get("c_kinase", 0) > 0 else [],
          "tumor_intrinsic": [evidence_id("purity", g)] if pd.notna(row.get("rna_lfc_at_purity1")) else [],
          "cis_dosage": [evidence_id("cis_dosage", g)] if row.get("c_cnv", 0) > 0 else []}
    d = (details.get(g) or {}).get("open_targets")
    if d:
        ev["tractability"] = [evidence_id("open_targets_tract", g)] + [
            evidence_id("open_targets_drug", x["chembl"], g) for x in d["drugs"]][:5]
    valid = set(store.check_evidence_ids([i for v in ev.values() for i in v])["valid"])
    return {k: [i for i in v if i in valid] for k, v in ev.items()}


def _r(x, nd=3):
    return None if x is None or pd.isna(x) else round(float(x), nd)


EXTRA_FIELDS = {  # M13 per-gene values carried into convergence.json (agents and the report quote them)
    "mutation_freq": 3, "mutation_padj": None, "cis_rho": 3, "depmap_mean_effect_ccrcc": 3, "depmap_selectivity": 3,
    "tcga_log2fc": 3, "tcga_padj": None, "tcga_cox_log_hr": 3, "tcga_cox_padj": None, "kinase_z": 3,
    "rna_lfc_at_purity1": 3, "rna_lfc_at_purity1_se": 3, "rna_purity_rho": 3}


def _row(i, row, cfg, details, ot, store, sens: dict) -> dict:
    d = (details.get(row["gene"]) or {}).get("open_targets") or {}
    comp_keys = list(cfg["convergence"]["stream_weights"])
    for prof in cfg.get("profiles", {}).values():
        comp_keys += [k for k in prof if k not in comp_keys]
    extra = {}
    for k, nd in EXTRA_FIELDS.items():
        v = row.get(k)
        extra[k] = None if v is None or pd.isna(v) else (float(f"{v:.3g}") if nd is None else round(float(v), nd))
    return {
        "rank": i + 1, "gene": row["gene"], "score": round(float(row["score"]), 4), "n_streams": int(row["n_streams"]),
        "components": {k: round(float(row[f"c_{k}"]), 3) for k in comp_keys if f"c_{k}" in row},
        "profile_scores": {p: round(float(row[f"score_{p}"]), 4) for p in cfg.get("profiles", {}) if f"score_{p}" in row},
        "rank_interval_90": sens.get(row["gene"]),
        **extra,
        "rna_log2fc": None if pd.isna(row.get("rna_log2fc")) else round(float(row["rna_log2fc"]), 3),
        "protein_log2fc": None if pd.isna(row.get("protein_log2fc")) else round(float(row["protein_log2fc"]), 3),
        "cnv_loss_freq": None if pd.isna(row.get("loss_freq")) else round(float(row["loss_freq"]), 3),
        "cnv_gain_freq": None if pd.isna(row.get("gain_freq")) else round(float(row["gain_freq"]), 3),
        "open_targets_score": round(float(row["ot_score"]), 4),
        "ot_genetic_association": round(float(row["ot_genetic"]), 4),
        "ot_somatic_mutation": round(float(row["ot_somatic"]), 4),
        "ot_literature": round(float(row["ot_literature"]), 4),
        "n_literature_relations": int(row["n_literature"]),
        "drugs": [{"name": x["name"], "max_stage": x["max_stage"]} for x in d.get("drugs", [])][:5],
        "evidence_ids": _evidence(row, ot, details, store),
    }


def write(n: int | None = None, day: str | None = None) -> dict:
    from .scores import write_profiles
    cfg = load()
    n = n or max(cfg["convergence"]["top_n_eval"], cfg["convergence"]["top_n_report"])
    t = compute()
    store = GraphStore()
    ot = pd.read_parquet(PROCESSED / "external_ot_associations.parquet").set_index("symbol")
    details = json.loads((PROCESSED / "external_targets.json").read_text())
    sens_summary = write_profiles(t) if cfg["convergence"].get("version", 2) >= 3 else {}

    def intervals(name):
        p = PROCESSED / f"sensitivity_{name}.parquet"
        if not p.exists():
            return {}
        s = pd.read_parquet(p)
        out = {r.gene: {"rank_p05": int(r.rank_p05), "rank_median": float(r.rank_median), "rank_p95": int(r.rank_p95),
                        "p_top10": round(float(r.p_top10), 3), "p_top25": round(float(r.p_top25), 3),
                        "p_top100": round(float(r.p_top100), 3)} for r in s.itertuples()}
        lp = PROCESSED / f"loso_{name}.parquet"
        if lp.exists():
            for r in pd.read_parquet(lp).itertuples():
                if r.gene in out:
                    out[r.gene]["leaves_top25_without"] = list(r.leaves_top_without)
        return out

    ranking = [_row(i, row, cfg, details, ot, store, intervals("convergence")) for i, row in t.head(n).iterrows()]
    profiles = {}
    for name in cfg.get("profiles", {}):
        tp = t.sort_values(f"score_{name}", ascending=False).reset_index(drop=True)
        rows = []
        for i, row in tp.head(cfg["convergence"]["top_n_report"]).iterrows():
            r = _row(i, row, cfg, details, ot, store, intervals(name))
            r["score"] = round(float(row[f"score_{name}"]), 4)
            rows.append(r)
        profiles[name] = rows
    out = {"generated": date.today().isoformat(), "config": cfg["convergence"], "profile_weights": cfg.get("profiles", {}),
           "sensitivity_config": cfg.get("sensitivity", {}), "sensitivity_summary": sens_summary,
           "universe_size": int(len(t)), "ranking": ranking, "profiles": profiles}
    path = run_dir(day) / "convergence.json"
    if path.exists():
        # Rebuilding a frozen run must not dirty it with a new date when nothing else changed.
        old = json.loads(path.read_text())
        if {k: v for k, v in old.items() if k != "generated"} == json.loads(json.dumps({k: v for k, v in out.items() if k != "generated"})):
            out["generated"] = old.get("generated", out["generated"])
    path.write_text(json.dumps(out, indent=1))
    from .ranking_evidence import load as load_ranking_evidence
    ranking_ids = load_ranking_evidence(out, store)
    return {"path": str(path), "top10": [r["gene"] for r in ranking[:10]], "universe": len(t),
            "ranking_evidence_ids": ranking_ids,
            **{f"{k}_top10": [r["gene"] for r in v[:10]] for k, v in profiles.items()}}

