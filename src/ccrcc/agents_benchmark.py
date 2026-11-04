"""Numbers for the report's "Agent benchmark" tab, read from the benchmark result files (nothing typed by hand).

Sources: benchmarks/architecture_v2/runs*/results.json (interpretation benchmark, one folder per condition) and
benchmarks/architecture_pipeline/results.json (full pipeline runs). Missing files give an empty dict and the tab
is not shown.
"""
from __future__ import annotations

import json
from pathlib import Path

from .paths import ROOT

V2 = ROOT / "benchmarks" / "architecture_v2"
PIPE = ROOT / "benchmarks" / "architecture_pipeline"
ARMS = {"B": "One agent", "C": "Specialists only", "D": "Specialists + orchestrator"}
CONDITIONS = {  # folder -> (label, structured output?, first attempt only?)
    "runs_first_attempt": ("4 per call, first attempt", False, True),
    "runs": ("4 per call, after one retry", False, False),
    "runs_load12": ("12 per call, after one retry", False, False),
    "runs_structured": ("4 per call, structured output", True, False),
    "runs_structured_load12": ("12 per call, structured output", True, False),
}
PIPE_CONFIGS = {  # arm field in session.json -> label
    "multi": "Multi-agent, before",
    "multi_improved": "Multi-agent, improved",
    "multi_improved_sonnet": "Improved, cheaper orchestrator (Sonnet)",
    "single": "One agent (Opus)",
}


def _load(p: Path) -> dict | None:
    return json.loads(p.read_text(encoding="utf-8")) if p.exists() else None


def _arm(v: dict) -> dict:
    return {"accuracy": round(100 * v["verdict_accuracy"]), "answered": round(100 * v["answered"]),
            "confident_errors": v["confident_errors_per_case"], "cost_case": v["cost_usd_per_case"],
            "cost_total": v["cost_usd_total"], "wall_case": v["wall_s_per_case"],
            "failed": v["failed_calls"], "calls": v["calls"]}


def build() -> dict:
    conds = []
    for folder, (label, structured, first) in CONDITIONS.items():
        r = _load(V2 / folder / "results.json")
        if r:
            conds.append({"key": folder, "label": label, "structured": structured, "first_attempt": first,
                          "arms": {a: _arm(r["arms"][a]) for a in ARMS},
                          "majority": round(100 * r.get("majority_class_verdict_accuracy", 0))})
    p = _load(PIPE / "results.json")
    if not conds or not p:
        return {}
    by = {c["key"]: c for c in conds}
    runs = p["runs"]
    pipe = []
    for arm, label in PIPE_CONFIGS.items():
        rr = [r for r in runs if r.get("arm") == arm and not r.get("incomplete")]
        if not rr:
            continue
        mean = lambda k: sum(r[k] for r in rr) / len(rr)
        pipe.append({"arm": arm, "label": label, "runs": rr,
                     "cost": mean("cost_usd"), "minutes": mean("duration_min"), "checklist": mean("checklist_score"),
                     "unsupported": sum(r["unsupported"] for r in rr), "share": mean("finding_evidence_cited_share")})
    # total model spend of the benchmark calls themselves (retries included), for the "what it cost" note
    spend_v2 = sum(c["arms"]["B"]["cost_total"] + c["arms"]["D"]["cost_total"]
                   for c in conds if not c["first_attempt"])  # D already includes the specialist (C) calls
    spend_pipe = sum(r["cost_usd"] or 0 for r in runs)
    return {"conditions": conds, "by": by, "arms": ARMS, "pipeline": pipe,
            "checks": list(runs[0]["checklist"].keys()), "pipe_runs": runs,
            "checklist_max": runs[0]["checklist_max"], "spend_v2": spend_v2, "spend_pipe": spend_pipe,
            "check_labels": CHECK_LABELS}


CHECK_LABELS = {
    "ca9_in_targets": "CA9 is in the target list",
    "cd70_partly_immune_derived": "CD70 is said to be made partly by T cells",
    "hk2_not_replicated": "HK2 is said not to replicate in the confirmatory cohort",
    "epas1_selective_dependency": "EPAS1 is named as the selective dependency",
    "no_gene_separates_arms": "No gene separates the arms of the trials",
    "checkmate025_used": "The CheckMate 025 trial data are used",
    "col23a1_prioritised": "COL23A1 (a new candidate) is prioritised",
    "rab42_prioritised": "RAB42 (a new candidate) is prioritised",
    "psat1_pck1_biomarkers": "PSAT1 and PCK1 are listed as biomarkers",
    "coxfa4l2_discussed": "COXFA4L2 (NDUFA4L2) is discussed",
    "normal_proteome_used": "The normal-tissue proteome (GTEx) is used",
    "surface_location_reasoned": "Where the protein sits in the cell is reasoned about",
    "prism_drug_response_used": "The PRISM drug screen is used",
}
