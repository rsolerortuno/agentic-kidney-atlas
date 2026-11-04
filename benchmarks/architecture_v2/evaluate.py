"""Score the v2 runs against the rubric answers and compare quality with cost.

Primary: verdict accuracy and sub-question accuracy. Secondary: all seven answers right, confident errors (yes/no
where the rubric says otherwise), over-abstention (insufficient where the rubric is definite), citations of other
entities' records, failures, cost, output tokens and wall time. Differences between arms are paired and bootstrapped
over the 12 genes (native and counterfactual dossiers and the replicates of a gene are averaged first), so repeated
measurements of one gene are not counted as independent evidence.
"""
from __future__ import annotations

import json
import random
from collections import defaultdict
from pathlib import Path

from .build import FROZEN, QUESTIONS, VERDICTS
from .run import RUNS, SPECIALISTS

ARMS = ["B", "C", "D"]


def _calls(tag_dir: Path):
    """Cost, output tokens and wall time per arm for every batch."""
    out = defaultdict(lambda: {"cost": 0.0, "out_tokens": 0, "wall_s": 0.0, "calls": 0, "failed_calls": 0})
    for bd in sorted(tag_dir.glob("rep*/batch*")):
        rec = {p.stem: json.loads(p.read_text(encoding="utf-8")) for p in bd.glob("*.json")}
        def add(arm, names, wall):
            o = out[(bd.parent.name, arm)]
            for n in names:
                r = rec.get(n, {})
                o["cost"] += r.get("cost_usd") or 0
                o["out_tokens"] += (r.get("usage") or {}).get("output_tokens") or 0
                o["calls"] += 1
                o["failed_calls"] += r.get("status") != "ok"
            o["wall_s"] += wall
        spec = [f"C_{n}" for n in SPECIALISTS]
        spec_wall = max((rec.get(n, {}).get("elapsed_s") or 0) for n in spec)
        add("B", ["B_single"], rec.get("B_single", {}).get("elapsed_s") or 0)
        add("C", spec, spec_wall)
        add("D", spec + ["D_orchestrator"], spec_wall + (rec.get("D_orchestrator", {}).get("elapsed_s") or 0))
    return out


def score_rows(tag_dir: Path, gold: dict) -> list[dict]:
    rows = []
    for af in sorted(tag_dir.glob("rep*/answers.json")):
        ans = json.loads(af.read_text(encoding="utf-8"))["answers"]
        for arm in ARMS:
            for cid, g in gold.items():
                a = ans.get(arm, {}).get(cid) or {}
                sub = {q: a.get(q) == g[q] for q in QUESTIONS}
                cited = [e for e in a.get("evidence_ids", []) if isinstance(e, str)]
                target = set(g["target_record_ids"])
                rows.append({
                    "rep": af.parent.name, "arm": arm, "case": cid, "gene": g["gene_index"], "synthetic": g["is_synthetic"],
                    "answered": bool(a), "verdict_ok": a.get("verdict") == g["verdict"],
                    "subq_acc": sum(sub.values()) / len(QUESTIONS), "all_ok": all(sub.values()) and a.get("verdict") == g["verdict"],
                    **{f"q_{q}": ok for q, ok in sub.items()},
                    "confident_errors": sum(1 for q in QUESTIONS if a.get(q) in ("yes", "no") and a.get(q) != g[q]),
                    "over_abstentions": sum(1 for q in QUESTIONS if a.get(q) == "insufficient" and g[q] != "insufficient"),
                    "false_candidate": a.get("verdict") == "surface_target_candidate" and g["verdict"] != "surface_target_candidate",
                    "cited": len(cited), "cited_other_entity": sum(1 for e in cited if e not in target)})
    return rows


def bootstrap(rows, a, b, field, n=5000, seed=7):
    per = defaultdict(lambda: defaultdict(list))
    for r in rows:
        per[r["gene"]][r["arm"]].append(float(r[field]))
    d = [sum(v[a]) / len(v[a]) - sum(v[b]) / len(v[b]) for v in per.values() if v.get(a) and v.get(b)]
    if not d:
        return None
    rng = random.Random(seed)
    draws = sorted(sum(rng.choices(d, k=len(d))) / len(d) for _ in range(n))
    return {"diff": round(sum(d) / len(d), 4), "ci95": [round(draws[int(.025 * n)], 4), round(draws[int(.975 * n)], 4)],
            "genes": len(d)}


def summarize(tag: str = "main") -> dict:
    tag_dir = RUNS if tag == "main" else RUNS.parent / f"runs_{tag}"
    gold = json.loads((FROZEN / "gold.json").read_text(encoding="utf-8"))
    rows = score_rows(tag_dir, gold)
    calls = _calls(tag_dir)
    reps = sorted({r["rep"] for r in rows})
    arms = {}
    for arm in ARMS:
        rr = [r for r in rows if r["arm"] == arm]
        if not rr:
            continue
        mean = lambda k: round(sum(float(r[k]) for r in rr) / len(rr), 4)
        cost = sum(calls[(rep, arm)]["cost"] for rep in reps)
        correct = sum(r["verdict_ok"] for r in rr)
        arms[arm] = {
            "case_answers": len(rr), "answered": mean("answered"), "verdict_accuracy": mean("verdict_ok"),
            "subquestion_accuracy": mean("subq_acc"), "all_seven_correct": mean("all_ok"),
            "per_question": {q: mean(f"q_{q}") for q in QUESTIONS},
            "native_verdict_accuracy": round(sum(r["verdict_ok"] for r in rr if not r["synthetic"]) /
                                             max(1, sum(not r["synthetic"] for r in rr)), 4),
            "counterfactual_verdict_accuracy": round(sum(r["verdict_ok"] for r in rr if r["synthetic"]) /
                                                     max(1, sum(r["synthetic"] for r in rr)), 4),
            "confident_errors_per_case": mean("confident_errors"), "over_abstentions_per_case": mean("over_abstentions"),
            "false_target_candidates": sum(r["false_candidate"] for r in rr),
            "other_entity_citation_rate": round(sum(r["cited_other_entity"] for r in rr) / max(1, sum(r["cited"] for r in rr)), 4),
            "cost_usd_total": round(cost, 3), "cost_usd_per_case": round(cost / len(rr), 4),
            "cost_usd_per_correct_verdict": round(cost / correct, 4) if correct else None,
            "output_tokens_per_case": round(sum(calls[(rep, arm)]["out_tokens"] for rep in reps) / len(rr), 1),
            "wall_s_per_case": round(sum(calls[(rep, arm)]["wall_s"] for rep in reps) / len(rr), 2),
            "calls": sum(calls[(rep, arm)]["calls"] for rep in reps),
            "failed_calls": sum(calls[(rep, arm)]["failed_calls"] for rep in reps)}
    counts = defaultdict(int)
    for g in gold.values():
        counts[g["verdict"]] += 1
    out = {"tag": tag, "replicates": reps, "cases": len(gold), "arms": arms,
           "majority_class_verdict_accuracy": round(max(counts.values()) / len(gold), 4),
           "paired_differences": {f"{a}-{b}": {f: bootstrap(rows, a, b, f) for f in ("verdict_ok", "subq_acc", "all_ok")}
                                  for a, b in (("C", "B"), ("D", "B"), ("D", "C"))}}
    (tag_dir / "results.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    return out


if __name__ == "__main__":
    import sys
    print(json.dumps(summarize(sys.argv[1] if len(sys.argv) > 1 else "main"), indent=1))
