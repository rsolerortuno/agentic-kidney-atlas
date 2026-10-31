"""Deterministic convergence scoring (no LLM). Writes runs/<date>/convergence.json.

Every score component carries the evidence ids that justify it; only ids present in the graph are kept.
"""
from __future__ import annotations

import json
from datetime import date

import pandas as pd

from .config import load
from .graph import GraphStore, evidence_id
from .literature_entities import resolve_gene
from .paths import LITERATURE, PROCESSED, RUNS


def run_dir(day: str | None = None):
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

    omics_raw = (lw["rna"] * t["rna_sig"] + lw["protein"] * t["protein_sig"] + lw["phospho"] * t["phospho_sig"]
                 + lw["concordance"] * t["concordant"])
    t["c_omics"] = omics_raw / sum(lw.values())
    t["c_cnv"] = t["cnv_freq"].fillna(0)
    t["n_literature"] = t["gene"].map(lambda g: len(set(lit_ids.get(g, []))))
    t["ot_literature"] = t["gene"].map(ot["dt_literature"]).fillna(0)
    t["c_literature"] = pd.concat([(t["n_literature"] / 3).clip(upper=1), t["ot_literature"]], axis=1).max(axis=1)
    t["ot_score"] = t["gene"].map(ot["score"]).fillna(0)
    t["ot_genetic"] = t["gene"].map(ot["dt_genetic_association"]).fillna(0)
    t["ot_somatic"] = t["gene"].map(ot["dt_somatic_mutation"]).fillna(0)
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
          "genetics": ot_ev if row["c_genetics"] > 0 else [],
          "tractability": []}
    d = (details.get(g) or {}).get("open_targets")
    if d:
        ev["tractability"] = [evidence_id("open_targets_tract", g)] + [
            evidence_id("open_targets_drug", x["chembl"], g) for x in d["drugs"]][:5]
    valid = set(store.check_evidence_ids([i for v in ev.values() for i in v])["valid"])
    return {k: [i for i in v if i in valid] for k, v in ev.items()}


def write(n: int | None = None, day: str | None = None) -> dict:
    cfg = load()
    n = n or max(cfg["convergence"]["top_n_eval"], cfg["convergence"]["top_n_report"])
    t = compute()
    store = GraphStore()
    ot = pd.read_parquet(PROCESSED / "external_ot_associations.parquet").set_index("symbol")
    details = json.loads((PROCESSED / "external_targets.json").read_text())
    ranking = []
    for i, row in t.head(n).iterrows():
        d = (details.get(row["gene"]) or {}).get("open_targets") or {}
        ranking.append({
            "rank": i + 1, "gene": row["gene"], "score": round(float(row["score"]), 4), "n_streams": int(row["n_streams"]),
            "components": {k: round(float(row[f"c_{k}"]), 3) for k in cfg["convergence"]["stream_weights"]},
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
        })
    out = {"generated": date.today().isoformat(), "config": cfg["convergence"],
           "universe_size": int(len(t)), "ranking": ranking}
    path = run_dir(day) / "convergence.json"
    path.write_text(json.dumps(out, indent=1))
    return {"path": str(path), "top10": [r["gene"] for r in ranking[:10]], "universe": len(t)}
