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


def _family_note(genes: list[str]) -> str:
    """Gold genes from one gene family usually share one drug class; the bootstrap treats them as independent."""
    fam: dict[str, list[str]] = {}
    for g in genes:
        fam.setdefault(re.sub(r"\d+[A-Z]?$", "", g), []).append(g)
    shared = [v for v in fam.values() if len(v) > 1]
    if not shared:
        return "All gold genes come from different gene families."
    return (f"{', '.join('/'.join(v) for v in shared)} share one gene family and drug class, so the set holds about "
            f"{len(fam)} independent targets, not {len(genes)}; the bootstrap interval treats every gene as independent "
            "and is too narrow.")


def _fmt_metrics(block: dict) -> list[str]:
    lines = ["| Score | Gold genes | AUROC (95 % CI) | Recall@100 | Recall@500 | Median rank | Enrichment@100 | Gold genes in top 100 |",
             "|---|---|---|---|---|---|---|---|"]
    for label, m in block.items():
        if not m["n_gold"]:
            lines.append(f"| {label} | 0 | — | — | — | — | — | — |")
            continue
        ci = f"{m['auroc']:.3f} ({m['auroc_ci95'][0]:.2f}–{m['auroc_ci95'][1]:.2f})"
        lines.append(f"| {label} | {m['n_gold']} | {ci} | {m['recall_at_100']:.0%} | {m['recall_at_500']:.0%} | "
                     f"{m['median_rank']:,} | {m['fold_enrichment_at_100']}× | {', '.join(m['top_hits']) or '—'} |")
    return lines


def write_results(rec: dict, qa: dict | None, outcomes: dict | None = None, baselines: dict | None = None,
                  sensitivity: dict | None = None, sabotage: dict | None = None) -> Path:
    lines = ["# Evaluation results", "",
             "Convergence **v3** (M13). Earlier versions and why each change was made: [history.md](history.md).", "",
             f"## (a) Recall of known ccRCC genes in the top {rec['top_n']} convergence genes", "",
             f"- Known genes: {rec['known']} (curated 10 + Open Targets top 20, deduplicated)",
             f"- **Recall: {rec['recall']:.0%}** (v2: 22 %; target 60 %)",
             f"- Recall with the genetics stream removed: {rec['recall_without_open_targets_stream']:.0%}",
             f"- recall@50/100/200/500/1000 = "
             + " / ".join(f"{v:.0%}" for v in rec["recall_at"].values())
             + f"; median rank {rec['median_rank']:.0f} of {rec['universe']:,} genes; "
             f"random expectation {rec['random_expected_recall']:.2%} → {rec['fold_enrichment']:.0f}× enrichment",
             "", "| Gene | Source | Rank | In top | In top without genetics stream |", "|---|---|---|---|---|"]
    for r in rec["genes"]:
        lines.append(f"| {r['gene']} | {r['source']} | {r['rank']} | {'✅' if r['in_top'] else '❌'} | "
                     f"{'✅' if r['in_top_without_ot'] else '❌'} |")
    if outcomes:
        g = outcomes["gold_sets"]
        lines += ["", "## (b) Outcomes the scores do not see (M13)", "",
                  "The target and biomarker profiles use measurements only (no literature, no drug history), so they can "
                  "be scored against clinical outcomes. The overall convergence score includes literature and drug "
                  "history and is shown only to make that circularity visible. Gold-set construction: "
                  "`src/ccrcc/outcomes.py`.", "",
                  f"### RCC clinical targets ({len(g['rcc_clinical_targets'])} genes: phase ≥ 2 in renal cell carcinoma, "
                  "≥ 2 trials or approval, drugs with ≤ 5 targets)", ""] + _fmt_metrics(outcomes["clinical_targets"])
        lines += ["", f"Predicted tractability (strong SM / AB buckets) flags {outcomes.get('n_predicted_tractable', '—'):,} genes. "
                  "It alone gives most of the AUROC — clinical targets are, almost by definition, druggable proteins — so "
                  "the data layers should be judged by early enrichment: the share of gold genes in the top 100, which within "
                  "the tractable set would be ~"
                  f"{100 / max(outcomes.get('n_predicted_tractable', 1), 1):.1%} per tractable gold gene by chance, and by the "
                  "\"without tractability\" row."]
        if "target profile, tractability without ligand history" in outcomes["clinical_targets"]:
            a = outcomes["clinical_targets"]["target profile"]["auroc"]
            b = outcomes["clinical_targets"]["target profile, tractability without ligand history"]["auroc"]
            lines += ["", "The strong small-molecule bucket includes Open Targets' \"High-Quality Ligand\", which comes from "
                      "ChEMBL medicinal chemistry: drug-discovery history, not a property of the tumour. Without it (pocket "
                      f"and location buckets only, {outcomes.get('n_predicted_tractable_no_ligand', 0):,} genes) the target "
                      f"profile's AUROC is {b} instead of {a}. Much of the gap between the data-only rows and the full "
                      "profile is therefore known chemistry, and the headline AUROC should be read with that in mind."]
        lines += ["", f"### Retrospective post-2016 evaluation: first RCC trial on or after {g['cutoff']} "
                  f"({len(g['temporal_holdout_targets'])} genes: {', '.join(g['temporal_holdout_targets'])})", "",
                  "Other eligible clinical targets, including those with incomplete dates, are removed from the background. First observed RCC trial dates use all phases. Features are not frozen before 2016: this is retrospective, not prospective validation.",
                  "",
                  _family_note(g["temporal_holdout_targets"]), ""] + _fmt_metrics(outcomes["temporal_holdout"])
        lines += ["", f"### Prognostic biomarkers: ClearCode34 ({len(g['clearcode34'])} of 34 genes mapped; Brooks et al. 2014)", "",
                  "ClearCode34 was validated in TCGA, which the biomarker profile's prognosis component also uses; the last row "
                  "removes that component. A literature-based gold set (verbatim `biomarker_of` relations from the disease-level "
                  f"query) was empty ({len(g['literature_biomarkers_disease_query'])} genes): all such relations came from "
                  "gene-specific queries chosen from the omics ranking, so they were not used.", ""] + _fmt_metrics(outcomes["clearcode34"])
    if sensitivity:
        lines += ["", "## (c) Robustness to the weights", "",
                  "1,000 random re-weightings (Dirichlet, κ = 20, around the configured weights). Genes in the top 10 in "
                  "≥ 80 % of draws:", ""]
        for name, v in sensitivity.items():
            lines.append(f"- **{name}**: {', '.join(v['robust_top10']) or 'none'}")
        if any("top25_leaving_without_one_stream" in v for v in sensitivity.values()):
            lines += ["", "Re-weighting never removes a stream. Leaving one stream out at a time (its weight set to 0), "
                      "these top-25 genes drop out of the top 25 (the stream they depend on in brackets):", ""]
            for name, v in sensitivity.items():
                dep = v.get("top25_leaving_without_one_stream", {})
                lines.append(f"- **{name}**: " + (", ".join(f"{g} ({'/'.join(s)})" for g, s in dep.items()) or "none"))
    if baselines:
        cb, pl = baselines["closed_book"], baselines["pipeline"]
        lines += ["", "## (d) What the agentic layer adds", "",
                  "| | Closed-book LLM (Sonnet 5.5, no tools) | GraphRAG pipeline |", "|---|---|---|",
                  f"| Runs | {cb['runs']} | {len(pl['runs'])} |",
                  f"| Mean pairwise Jaccard of candidate lists | {cb['mean_pairwise_jaccard']} | {pl['mean_pairwise_jaccard']} |",
                  f"| Candidates that are RCC clinical targets | {', '.join(f'{x['in_rcc_clinical_targets']:.0%}' for x in cb['support'])} | "
                  f"{', '.join(f'{x['in_rcc_clinical_targets']:.0%}' for x in pl['support'])} |",
                  f"| Candidates in the data-only target or biomarker top 500 | "
                  f"{', '.join(f'{x['in_target_or_biomarker_profile_top500']:.0%}' for x in cb['support'])} | "
                  f"{', '.join(f'{x['in_target_or_biomarker_profile_top500']:.0%}' for x in pl['support'])} "
                  "(by construction: the agents pick candidates from these profiles) |",
                  f"| Statements citing checkable evidence | 0 | {pl.get('statements_checked', '—')} checked, "
                  f"{pl.get('unsupported', '—')} unsupported |", "",
                  f"Closed-book candidates: {'; '.join(', '.join(c) for c in cb['candidates'])}.", "",
                  f"Pipeline final targets: {'; '.join(', '.join(c) for c in pl['targets'])}.", "",
                  "The pipeline's Jaccard mostly reflects the deterministic ranking the agents start from, not agreement "
                  "reached independently; it shows the agent step is stable, not that its choices are right."]
        nb = baselines.get("no_agent_baseline")
        if nb:
            lines += ["", "### No-agent baseline", "",
                      f"Rule: {nb['rule']}. Overlap with the agents' final list (Jaccard): {nb['jaccard_with_agents']}. "
                      f"RCC clinical targets: {nb['support']['in_rcc_clinical_targets']:.0%} of the baseline list versus "
                      f"{pl['support'][0]['in_rcc_clinical_targets']:.0%} of the agents' list.", "",
                      f"Baseline list: {', '.join(nb['genes'])}.", "",
                      f"Added by the agents: {', '.join(nb['added_by_agents']) or 'none'}. "
                      f"Left out by the agents: {', '.join(nb['dropped_by_agents']) or 'none'}.", "",
                      "This is what the agent step changes. Whether the changes are improvements needs the reasons in "
                      "the report and, ultimately, experiments; the overlap alone does not say."]
    if sabotage:
        c = sabotage["corruptions"]
        pct = lambda x: "—" if x is None else f"{x:.0%}"
        lines += ["", "## (e) How often the verifier catches a wrong statement", "",
                  f"Each of the {sabotage['statements_passing']} statements of `runs/{sabotage['run']}` that pass "
                  "verification was copied and corrupted one way at a time; a corruption counts as caught when the "
                  "statement is rejected or loses evidence ids (the rule `verify.py final` applies).", "",
                  "| Corruption | Applied | Caught | Rate |", "|---|---|---|---|",
                  f"| Number: last digit of the first printed number changed | {c['number']['applied']} | {c['number']['caught']} | {pct(c['number']['rate'])} |",
                  f"| Direction word flipped (higher/lower, up/down, ...) | {c['direction']['applied']} | {c['direction']['caught']} | {pct(c['direction']['rate'])} |",
                  f"| ... of which the statement cites significant omics_de evidence | {c['direction']['applied_with_de_evidence']} | "
                  f"{c['direction']['caught_with_de_evidence']} | {pct(c['direction']['rate_with_de_evidence'])} |",
                  f"| Statement attributed to another gene | {c['gene']['applied']} | {c['gene']['caught']} | {pct(c['gene']['rate'])} |",
                  f"| Evidence of another gene cited | {c['evidence']['applied']} | {c['evidence']['caught']} | {pct(c['evidence']['rate'])} |", "",
                  "The direction check only reads tumour-versus-normal wording against cited omics_de records and printed "
                  "log2FC signs, so a flipped direction in a literature or treatment sentence is not caught. A changed last "
                  "digit is accepted when the new value still rounds to the stored one. Missed examples are listed in "
                  "`evals/verifier_sabotage.json`."]
    lines += ["", "## (f) Known-answer questions (`claude -p` + project MCP servers)", ""]
    if qa:
        lines += [f"- Model: {qa['model']}", f"- **Correct: {qa['correct']}/{qa['total']}**", "",
                  "| ID | Question | Expected | Correct | Answer / error |", "|---|---|---|---|---|"]
        for r in qa["results"]:
            ans = (r["error"] and f"ERROR: {r['error']}") or r["answer"]
            ans = ans.replace("|", "\\|").replace("\n", " ")[:220]
            lines.append(f"| {r['id']} | {r['question']} | {r['expected']} | {'✅' if r['correct'] else '❌'} | {ans} |")
    else:
        lines.append("Not run in this evaluation.")
    out = EVALS / "results.md"
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    (EVALS / "results.json").write_text(json.dumps({"recall": rec, "qa": qa, "outcomes": outcomes,
                                                    "sensitivity": sensitivity, "baselines": baselines},
                                                   indent=1, default=str), encoding="utf-8")
    return out


def rebuild_saved_evaluation(run: Path | None = None) -> Path:
    """Recompute deterministic evaluations; reuse only the frozen model answers."""
    from .baselines import run as compare_baselines
    from .outcomes import run as evaluate_outcomes

    qa = grade_answers(EVALS / "qa_answers.json", "saved claude-sonnet-5-5 answers (regex regraded)")
    sensitivity = json.loads((PROCESSED / "sensitivity_summary.json").read_text(encoding="utf-8"))
    from . import sabotage
    sab = sabotage.run(run) if run is not None and (Path(run) / "findings.json").exists() else None
    return write_results(recall(), qa, evaluate_outcomes(compute()),
                         compare_baselines(cached_only=True, official=run), sensitivity, sab)
