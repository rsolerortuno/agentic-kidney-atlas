"""M8: evaluation.

(a) Recall of known ccRCC genes (curated list + Open Targets top 20) in the top-N convergence ranking.
    Because the Open Targets association is itself one scoring stream, recall is also reported for an
    ablated score without the Open Targets stream (to expose circularity).
(b) Known-answer questions (evals/qa.yaml) answered by Claude Code headless (`claude -p`) with the project
    MCP servers; graded by regex.
"""
from __future__ import annotations

import json
import re
import shutil
import subprocess
from pathlib import Path

import pandas as pd
import yaml

from .config import load
from .integrate import compute
from .paths import PROCESSED, ROOT

CURATED = ["VHL", "PBRM1", "SETD2", "BAP1", "KDM5C", "EPAS1", "CA9", "VEGFA", "MTOR", "HIF1A"]
EVALS = ROOT / "evals"


def known_genes() -> list[str]:
    ot = pd.read_parquet(PROCESSED / "external_ot_associations.parquet")
    return list(dict.fromkeys(CURATED + ot["symbol"].head(20).tolist()))


def recall(top_n: int | None = None) -> dict:
    cfg = load()["convergence"]
    top_n = top_n or cfg["top_n_eval"]
    t = compute()
    known = known_genes()
    full_top = t["gene"].head(top_n).tolist()
    sw = cfg["stream_weights"]
    ablated = t.assign(s=sum(sw[k] * t[f"c_{k}"] for k in sw if k != "genetics"))
    abl_top = ablated.sort_values("s", ascending=False)["gene"].head(top_n).tolist()
    rank = {g: i + 1 for i, g in enumerate(t["gene"])}
    rows = [{"gene": g, "source": "curated" if g in CURATED else "open_targets_top20", "rank": rank.get(g),
             "in_top": g in full_top, "in_top_without_ot": g in abl_top} for g in known]
    ranks = [r["rank"] for r in rows if r["rank"]]
    universe = len(t)
    rec = sum(r["in_top"] for r in rows) / len(known)
    return {"top_n": top_n, "known": len(known), "universe": universe,
            "recall": rec,
            "recall_without_open_targets_stream": sum(r["in_top_without_ot"] for r in rows) / len(known),
            "recall_at": {k: sum(1 for x in ranks if x <= k) / len(known) for k in (50, 100, 200, 500, 1000)},
            "median_rank": float(pd.Series(ranks).median()),
            "random_expected_recall": top_n / universe,
            "fold_enrichment": rec / (top_n / universe),
            "genes": rows}


def find_claude() -> str | None:
    exe = shutil.which("claude")
    if exe:
        return exe
    import os
    base = Path(os.environ.get("LOCALAPPDATA", "")) / "Claude-Data"
    found = sorted(base.glob("*/claude-code/*/*/claude.exe"), key=lambda p: p.stat().st_mtime, reverse=True)
    return str(found[0]) if found else None


QA_PROMPT = ("You are answering a factual question about this project's ccRCC analysis. Use ONLY the MCP tools "
             "(omics-mcp, kg-mcp, external-mcp) to look it up; do not answer from memory. "
             "Reply with one short sentence.\nQuestion: {q}")


def run_qa(model: str = "claude-sonnet-5-5", timeout: int = 300) -> dict:
    qa = yaml.safe_load((EVALS / "qa.yaml").read_text(encoding="utf-8"))
    claude = find_claude()
    allowed = "mcp__omics-mcp,mcp__kg-mcp,mcp__external-mcp"
    results = []
    for item in qa:
        answer, error = "", None
        if not claude:
            error = "claude CLI not found"
        else:
            try:
                p = subprocess.run([claude, "-p", QA_PROMPT.format(q=item["question"]), "--model", model,
                                    "--mcp-config", str(ROOT / ".mcp.json"), "--allowedTools", allowed],
                                   capture_output=True, text=True, timeout=timeout, cwd=ROOT, encoding="utf-8")
                answer = (p.stdout or "").strip()
                if p.returncode != 0 or "Not logged in" in answer:
                    error = answer[:200] or p.stderr[:200]
            except subprocess.TimeoutExpired:
                error = "timeout"
        correct = bool(not error and re.search(item["pattern"], answer, re.I))
        results.append({"id": item["id"], "question": item["question"], "expected": item["expected"],
                        "answer": answer[:500], "correct": correct, "error": error})
    return {"model": model, "correct": sum(r["correct"] for r in results), "total": len(results), "results": results}


def grade_answers(path: Path, model: str) -> dict:
    """Grade a saved answers file ([{"id", "answer"}]) with the same regexes as run_qa."""
    qa = {q["id"]: q for q in yaml.safe_load((EVALS / "qa.yaml").read_text(encoding="utf-8"))}
    answers = {a["id"]: a["answer"] for a in json.loads(Path(path).read_text(encoding="utf-8"))}
    results = []
    for qid, item in qa.items():
        ans = answers.get(qid, "")
        results.append({"id": qid, "question": item["question"], "expected": item["expected"], "answer": ans[:500],
                        "correct": bool(re.search(item["pattern"], ans, re.I)), "error": None if ans else "no answer"})
    return {"model": model, "correct": sum(r["correct"] for r in results), "total": len(results), "results": results}


def write_results(rec: dict, qa: dict | None) -> Path:
    lines = ["# Evaluation results", "",
             f"## (a) Recall of known ccRCC genes in the top {rec['top_n']} convergence genes", "",
             f"- Known genes: {rec['known']} (curated 10 + Open Targets top 20, deduplicated)",
             f"- **Recall: {rec['recall']:.0%}**",
             f"- Recall with the Open Targets genetics stream removed (circularity check): "
             f"{rec['recall_without_open_targets_stream']:.0%}",
             f"- Supplementary: recall@50/100/200/500/1000 = "
             + " / ".join(f"{v:.0%}" for v in rec["recall_at"].values())
             + f"; median rank {rec['median_rank']:.0f} of {rec['universe']:,} genes; "
             f"random expectation {rec['random_expected_recall']:.2%} → {rec['fold_enrichment']:.0f}× enrichment",
             "", "| Gene | Source | Rank | In top | In top without genetics stream |", "|---|---|---|---|---|"]
    for r in rec["genes"]:
        lines.append(f"| {r['gene']} | {r['source']} | {r['rank']} | {'✅' if r['in_top'] else '❌'} | "
                     f"{'✅' if r['in_top_without_ot'] else '❌'} |")
    lines += ["", "## (b) Known-answer questions (`claude -p` + project MCP servers)", ""]
    if qa:
        lines += [f"- Model: {qa['model']}", f"- **Correct: {qa['correct']}/{qa['total']}**", "",
                  "| ID | Question | Expected | Correct | Answer / error |", "|---|---|---|---|---|"]
        for r in qa["results"]:
            ans = (r["error"] and f"ERROR: {r['error']}") or r["answer"]
            ans = ans.replace("|", "\\|").replace("\n", " ")[:220]
            lines.append(f"| {r['id']} | {r['question']} | {r['expected']} | {'✅' if r['correct'] else '❌'} | {ans} |")
    else:
        lines.append("Not run.")
    lines += ["", "## Analysis", "",
              "Recall@100 is below the 60 % target. The cause is structural, not a bug: see "
              "[history.md](history.md) for the evaluated versions, the single methodological change, and why "
              "mutational drivers with unchanged bulk expression (VHL, PBRM1, SETD2, BAP1, KDM5C…) cannot reach "
              "the top 100 of an expression-convergence ranking without a somatic-mutation layer."]
    out = EVALS / "results.md"
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    (EVALS / "results.json").write_text(json.dumps({"recall": rec, "qa": qa}, indent=1), encoding="utf-8")
    return out
