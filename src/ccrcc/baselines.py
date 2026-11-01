"""M13: what does the agentic layer add? Three comparisons (outputs in evals/baselines/).

1. Closed-book LLM. The same question is put to Claude Sonnet 5.5 (the subagent model) through `claude -p` with no
   tools, no MCP servers and no project files, three times. Its candidates are compared with the pipeline's on
   reproducibility (Jaccard between runs), on overlap with RCC clinical targets (what memorized knowledge gives),
   on support by the data-only profiles (top 500) and on traceability (statements that cite evidence: none can).
2. Agents vs deterministic ranking. Which candidates the orchestrator promoted or left out relative to the numeric
   convergence rank, next to the M13 values behind those calls (purity-extrapolated log2FC, selective dependency,
   replication, target and biomarker profile scores).
3. Run-to-run variability of the full agentic pipeline: the official run plus replicate runs of `/informe`
   (CCRCC_RUN_DIR=replicates/repN), compared by Jaccard and top-k agreement of the final target lists.
"""
from __future__ import annotations

import itertools
import json
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

import pandas as pd

from .config import load
from .literature_entities import resolve_gene
from .paths import PROCESSED, ROOT, RUNS

OUT = ROOT / "evals" / "baselines"
CLOSED_PROMPT = (
    "{question}\nAnswer from your own knowledge only: you have no tools, files or databases. Return ONLY a JSON object, "
    "no prose and no code fences: {{\"candidates\": [{{\"gene\": \"<HGNC symbol>\", \"role\": \"target|biomarker|"
    "target_and_biomarker\", \"rationale\": \"one sentence, including any quantitative facts you know\"}}]}} "
    "with exactly 10 candidates ordered by strength of evidence.")
NO_TOOLS = "Bash,Read,Write,Edit,Glob,Grep,WebFetch,WebSearch,Agent,Task,NotebookEdit,TodoWrite"


def _claude() -> str | None:
    return shutil.which("claude")


def closed_book(runs: int = 3, model: str = "claude-sonnet-5-5", timeout: int = 600, cached_only: bool = False) -> list[dict]:
    OUT.mkdir(parents=True, exist_ok=True)
    q = load()["question"]
    claude = _claude()
    results = []
    for i in range(1, runs + 1):
        path = OUT / f"closed_book_{i}.json"
        if path.exists():
            results.append(json.loads(path.read_text(encoding="utf-8")))
            continue
        if cached_only:
            raise FileNotFoundError(f"Frozen model response missing: {path}")
        if not claude:
            raise RuntimeError("claude CLI not found")
        with tempfile.TemporaryDirectory() as tmp:   # no .mcp.json, no CLAUDE.md, no project files
            p = subprocess.run([claude, "-p", CLOSED_PROMPT.format(question=q), "--model", model,
                                "--disallowedTools", NO_TOOLS, "--strict-mcp-config"],
                               capture_output=True, text=True, timeout=timeout, cwd=tmp, encoding="utf-8")
        raw = (p.stdout or "").strip()
        m = re.search(r"\{.*\}", raw, re.S)
        try:
            data = json.loads(m.group(0)) if m else {"candidates": []}
        except json.JSONDecodeError:
            data = {"candidates": [], "parse_error": True}
        rec = {"run": i, "model": model, "raw": raw[:6000], **data}
        path.write_text(json.dumps(rec, indent=1, ensure_ascii=False), encoding="utf-8")
        results.append(rec)
    return results


def _genes(cands) -> list[str]:
    out = []
    for c in cands:
        g = resolve_gene(str(c.get("gene", ""))) or str(c.get("gene", "")).upper()
        if g and g not in out:
            out.append(g)
    return out


def jaccard(a, b) -> float:
    a, b = set(a), set(b)
    return len(a & b) / len(a | b) if a | b else float("nan")


def _mean_pairwise(lists) -> float | None:
    pairs = list(itertools.combinations(lists, 2))
    return round(sum(jaccard(a, b) for a, b in pairs) / len(pairs), 3) if pairs else None


def _support(genes: list[str], scores: pd.DataFrame, clinical: set[str]) -> dict:
    top_t = set(scores.sort_values("score_target", ascending=False)["gene"].head(500))
    top_b = set(scores.sort_values("score_biomarker", ascending=False)["gene"].head(500))
    n = len(genes) or 1
    return {"n": len(genes), "in_rcc_clinical_targets": round(sum(g in clinical for g in genes) / n, 3),
            "in_target_or_biomarker_profile_top500": round(sum(g in top_t or g in top_b for g in genes) / n, 3),
            "not_clinical_but_data_supported": sorted(g for g in genes if g not in clinical and (g in top_t or g in top_b))}


def _numbers(text: str) -> int:
    return len(re.findall(r"(?<![\w.])\d+(?:\.\d+)?(?:%|-fold)?", text or ""))


def run_dirs(official: Path | None = None) -> list[Path]:
    dated = sorted(p.parent for p in RUNS.glob("[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]/findings.json"))
    reps = sorted(p.parent for p in (RUNS / "replicates").glob("*/findings.json"))
    return ([official] if official is not None else dated[-1:]) + reps


def _final_targets(run: Path) -> list[str]:
    f = json.loads((run / "findings.json").read_text(encoding="utf-8"))
    return [t["gene"] for t in f.get("targets", [])]


def agents_vs_ranking(run: Path, scores: pd.DataFrame) -> list[dict]:
    """Promotions / demotions by the orchestrator relative to the numeric convergence rank, with the M13 values."""
    f = json.loads((run / "findings.json").read_text(encoding="utf-8"))
    conv = f["convergence"]["ranking"]
    rank = {r["gene"]: r["rank"] for r in conv}
    s = scores.set_index("gene")
    out = []
    targets = [t["gene"] for t in f.get("targets", [])]
    top10 = [r["gene"] for r in conv[:10]]
    for g in list(dict.fromkeys(targets + top10)):
        row = s.loc[g] if g in s.index else None
        def val(c):
            return None if row is None or pd.isna(row.get(c)) else round(float(row[c]), 3)
        out.append({"gene": g, "agent_position": targets.index(g) + 1 if g in targets else None,
                    "convergence_rank": rank.get(g), "decision": ("kept" if g in targets and g in top10 else
                                                                  "promoted" if g in targets else "left out"),
                    "rna_log2fc": val("rna_log2fc"), "rna_lfc_at_purity1": val("rna_lfc_at_purity1"),
                    "c_tumor_intrinsic": val("c_tumor_intrinsic"), "c_dependency": val("c_dependency"),
                    "c_replication": val("c_replication"), "c_mutation": val("c_mutation"),
                    "target_profile_score": val("score_target"), "biomarker_profile_score": val("score_biomarker")})
    return out


def run(closed_runs: int = 3, cached_only: bool = False, official: Path | None = None) -> dict:
    from .integrate import compute
    from .outcomes import clinical_targets
    OUT.mkdir(parents=True, exist_ok=True)
    scores = compute()
    clinical = set(clinical_targets()["gene"])
    closed = closed_book(closed_runs, cached_only=cached_only)
    closed_lists = [_genes(r.get("candidates", [])) for r in closed]
    runs = run_dirs(official)
    pipe_lists = [_final_targets(r) for r in runs]
    official = runs[0] if runs else None
    res = {
        "closed_book": {
            "runs": len(closed_lists), "candidates": closed_lists,
            "mean_pairwise_jaccard": _mean_pairwise(closed_lists),
            "support": [_support(g, scores, clinical) for g in closed_lists],
            "numbers_per_run": [sum(_numbers(c.get("rationale", "")) for c in r.get("candidates", [])) for r in closed],
            "statements_citing_evidence": 0,
        },
        "pipeline": {
            "runs": [str(r.relative_to(ROOT)) for r in runs], "targets": pipe_lists,
            "mean_pairwise_jaccard": _mean_pairwise(pipe_lists),
            "top3_identical_across_runs": len({tuple(x[:3]) for x in pipe_lists}) == 1 if pipe_lists else None,
            "support": [_support(g, scores, clinical) for g in pipe_lists],
        },
        "closed_vs_pipeline_jaccard": (round(jaccard(set().union(*closed_lists), pipe_lists[0]), 3)
                                       if closed_lists and pipe_lists else None),
        "agents_vs_ranking": agents_vs_ranking(official, scores) if official else [],
    }
    if official:
        vf = json.loads((official / "verification_final.json").read_text(encoding="utf-8"))
        res["pipeline"]["statements_checked"] = vf["statements"]
        res["pipeline"]["unsupported"] = vf["unsupported"]
    (OUT / "baselines.json").write_text(json.dumps(res, indent=1, default=str), encoding="utf-8")
    return res
