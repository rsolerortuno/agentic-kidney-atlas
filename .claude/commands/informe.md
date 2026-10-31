---
description: Orchestrate the ccRCC multi-agent analysis (4 subagents + deterministic integration + verification) and render the HTML report.
---

You are the orchestrator of the ccRCC GraphRAG analysis. Question:
"What are the candidate therapeutic targets and biomarkers in ccRCC with the strongest convergent evidence?"

Follow these steps exactly. All numbers come from tools/scripts; you never compute or invent values.

1. Run `uv run python scripts/integrate.py`. It prints the run directory path (`runs/<date>/convergence.json`).
   Call that directory RUN below.

2. Launch the four subagents **in parallel** (one message, four Agent/Task calls), each with this prompt
   (replace RUN with the absolute path):
   - `agente-omico`: "Question: <question>. Write your output to RUN/agent_omics.json."
   - `agente-literatura`: "Question: <question>. Write your output to RUN/agent_literature.json."
   - `agente-externo`: "Question: <question>. Write your output to RUN/agent_external.json."
   - `agente-celular`: "Question: <question>. Write your output to RUN/agent_cells.json."

3. Run `uv run python scripts/verify.py agents`. It removes evidence ids that do not exist or are not linked to
   the gene and drops unsupported findings, writing `*.verified.json`.

4. Read `RUN/convergence.json` (top 25) and the four `RUN/agent_*.verified.json` files. Integrate them: where do
   omics, literature and external knowledge converge, where do they disagree, which candidates are targets vs.
   biomarkers vs. mechanism. Then write `RUN/narrative.json`:

   ```json
   {"summary":  [{"text": "...", "evidence_ids": ["..."]}],
    "sections": {"omics": [...], "literature": [...], "external": [...], "cells": [...], "spatial": [...], "structure": [...], "discovery": [...], "integration": [...]},
    "targets":  [{"gene": "CA9", "category": "target|biomarker|target_and_biomarker",
                  "rationale": "...", "evidence_ids": ["..."]}],
    "limitations": [{"text": "..."}]}
   ```
   - `cells`, `spatial` and `structure` come from `agent_cells.verified.json`: which cell type expresses each
     candidate (tumor-intrinsic vs stromal / immune), what surrounds it in situ, and which candidates have a
     credible molecular-glue-degrader site. Say explicitly when a bulk-omics candidate is not tumor-intrinsic, and
     present MGD scores as a ranking of sites, never as evidence that a degrader exists.
   - `discovery` (2–3 paragraphs) comes from the `role: "discovery"` findings of `agent_cells.verified.json`: the most promising
     under-studied, undrugged candidates per track, what makes them novel, and the caveat that they are data-ranked hypotheses
     (no dependency or functional evidence).
   - 2–4 summary paragraphs; 1–3 paragraphs per section; 8–12 entries in `targets` ordered by your judgment of
     strength (you may depart from the numeric ranking, but say why in `rationale`).
   - Every summary/section paragraph and every target MUST cite evidence ids copied from the verified agent files
     or from `convergence.json`. A paragraph that discusses one specific gene should carry `"gene": "<symbol>"`
     and cite only that gene's evidence. Limitations need no ids.

5. Run `uv run python scripts/assemble.py` then `uv run python scripts/verify.py final`.
   If it reports any unsupported claim, fix `narrative.json` (remove or re-cite the statement) and repeat step 5
   until it prints `unsupported claims: 0`.

6. Run `uv run python scripts/render.py` and reply with the report path and a 5-line summary.
