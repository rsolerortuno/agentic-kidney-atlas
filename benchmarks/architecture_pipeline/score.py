"""Score full pipeline runs: multi-agent (/informe) vs single agent (/informe_single), both with the project MCP tools.

There is no gold answer for an open analysis, so quality is measured four ways, all fixed before the runs:
1. automatic verification (statements checked, unsupported, stale) and how many agent findings the verifier dropped
   on first pass;
2. coverage: evidence types and genes cited, sections filled, discovery genes in the target list;
3. a checklist of facts the data establish and a complete analysis should state (CHECKS, regexes on the narrative);
4. agreement of the final target list with the published multi-agent run.
Cost, tokens, turns and wall time come from Claude Code's JSON result.
"""
from __future__ import annotations

import json
import re
import sqlite3
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).parent
SECTIONS = ["omics", "literature", "external", "cells", "spatial", "structure", "discovery", "validation", "robustness",
            "treatment", "subtypes", "adc", "immune", "integration"]


def _para(text: str, *patterns: str) -> bool:
    """All patterns occur in one paragraph."""
    return any(all(re.search(p, para, re.I) for p in patterns) for para in text.split("\n\n"))


CHECKS = {
    "ca9_in_targets": lambda t, tg, ids: "CA9" in tg,
    "cd70_partly_immune_derived": lambda t, tg, ids: _para(t, r"\bCD70\b", r"\bT[ -]cells?\b|\bCD8\b|lymphocyt|immune"),
    "hk2_not_replicated": lambda t, tg, ids: _para(t, r"\bHK2\b", r"not replicat|does not replicate|did not replicate|fails? to replicate"),
    "epas1_selective_dependency": lambda t, tg, ids: _para(t, r"\bEPAS1\b", r"dependen"),
    "no_gene_separates_arms": lambda t, tg, ids: bool(re.search(r"no gene[^.]{0,80}(separat|differ|distinguish)", t, re.I)),
    # CheckMate 025 is the nivolumab-vs-everolimus trial; reports may name the arms without the trial name.
    "checkmate025_used": lambda t, tg, ids: bool(re.search(r"CheckMate", t, re.I)) or _para(t, r"\bnivolumab\b", r"\beverolimus\b"),
    "col23a1_prioritised": lambda t, tg, ids: "COL23A1" in tg,
    "rab42_prioritised": lambda t, tg, ids: "RAB42" in tg or _para(t, r"\bRAB42\b", r"discover|unexplored|exploratory"),
    "psat1_pck1_biomarkers": lambda t, tg, ids: "PSAT1" in tg and "PCK1" in tg,
    "coxfa4l2_discussed": lambda t, tg, ids: bool(re.search(r"COXFA4L2|NDUFA4L2", t)),
    "normal_proteome_used": lambda t, tg, ids: bool(re.search(r"GTEx|normal[- ]tissue prote|proteome", t, re.I)),
    "surface_location_reasoned": lambda t, tg, ids: bool(re.search(r"\binternal\b|intracellular|inside the cell|not (on|at) the (cell )?surface", t, re.I)),
    # PRISM counts when the narrative cites a PRISM record or names the screen; "compound" alone matches drug names.
    "prism_drug_response_used": lambda t, tg, ids: any(i.startswith("EV:prism_response:") for i in ids)
        or bool(re.search(r"\bPRISM\b|drug[- ]response|drug sensitivity", t, re.I)),
}


def narrative_text(n: dict) -> str:
    parts = [p.get("text", "") for p in n.get("summary", [])]
    for paras in (n.get("sections") or {}).values():
        parts += [p.get("text", "") for p in paras]
    parts += [t.get("rationale", "") for t in n.get("targets", [])]
    return "\n\n".join(parts)


def evidence_types(ids: set[str]) -> set[str]:
    return {i.split(":")[1] for i in ids if i.count(":") >= 2}


def score_run(run: Path, official_targets: list[str]) -> dict:
    n = json.loads((run / "narrative.json").read_text(encoding="utf-8"))
    vf = json.loads((run / "verification_final.json").read_text(encoding="utf-8"))
    va = json.loads((run / "verification_agents.json").read_text(encoding="utf-8"))["agents"]
    f = json.loads((run / "findings.json").read_text(encoding="utf-8"))
    targets = [t["gene"] for t in n.get("targets", [])]
    text = narrative_text(n)
    ids = set()
    for block in f.get("agent_findings", {}).values():
        for x in block.get("findings", []):
            ids |= set(x.get("evidence_ids", []))
    narr_ids = set()
    for p in n.get("summary", []) + [p for v in (n.get("sections") or {}).values() for p in v] + n.get("targets", []):
        narr_ids |= set(p.get("evidence_ids", []))
    genes = {x.get("gene") for b in f.get("agent_findings", {}).values() for x in b.get("findings", []) if x.get("gene")}
    checks = {k: bool(fn(text, targets, narr_ids)) for k, fn in CHECKS.items()}
    jac = len(set(targets) & set(official_targets)) / len(set(targets) | set(official_targets)) if targets else 0
    return {
        "run": run.name, "statements_checked": vf["statements"], "unsupported": vf["unsupported"], "stale": vf["stale"],
        "agent_findings_submitted": sum(a.get("submitted", 0) for a in va.values()),
        "agent_findings_kept_first_pass": sum(a.get("kept", 0) for a in va.values()),
        "agent_findings_final": sum(len(b.get("findings", [])) for b in f.get("agent_findings", {}).values()),
        "genes_with_findings": len(genes), "evidence_types_cited": len(evidence_types(ids | narr_ids)),
        "distinct_evidence_ids": len(ids | narr_ids),
        # integration loss: share of the specialists' verified evidence that the final narrative cites
        "finding_evidence_cited_share": round(len(ids & narr_ids) / len(ids), 3) if ids else None,
        "sections_filled": sum(1 for s in SECTIONS if (n.get("sections") or {}).get(s)),
        "targets": targets, "jaccard_with_published_run": round(jac, 3),
        "checklist": checks, "checklist_score": sum(checks.values()), "checklist_max": len(CHECKS),
        "narrative_words": len(text.split()),
    }


def summarize() -> dict:
    official = [t["gene"] for t in json.loads((ROOT / "runs/2026-11-03/narrative.json").read_text(encoding="utf-8"))["targets"]]
    out = {"published_run_targets": official, "runs": []}
    for meta in sorted((HERE / "runs").glob("*/session.json")):
        run = meta.parent
        s = json.loads(meta.read_text(encoding="utf-8"))
        rec = score_run(run, official) if (run / "verification_final.json").exists() else {"run": run.name, "incomplete": True}
        rec.update(arm=s["arm"], cost_usd=s.get("total_cost_usd"), duration_min=round((s.get("duration_ms") or 0) / 60000, 1),
                   turns=s.get("num_turns"), output_tokens=(s.get("usage") or {}).get("output_tokens"),
                   models=sorted((s.get("modelUsage") or {}).keys()))
        out["runs"].append(rec)
    (HERE / "results.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    return out


if __name__ == "__main__":
    print(json.dumps(summarize(), indent=1))
