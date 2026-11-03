---
description: Orchestrate the ccRCC multi-agent analysis (4 subagents + deterministic integration + verification) and render the HTML report.
---

You are the orchestrator of the ccRCC GraphRAG analysis. Question:
"What are the candidate therapeutic targets and biomarkers in ccRCC with the strongest convergent evidence?"

Follow these steps exactly. All numbers come from tools/scripts; you never compute or invent values.

1. Run `uv run python scripts/integrate.py`. It prints the run directory path (`runs/<date>/convergence.json`).
   Call that directory RUN below. It also records RUN in `runs/.current`, so the later scripts use RUN by default even
   if a step is repeated on another day. Never run this script again in the same analysis: it deletes the agent files.

2. Launch the five subagents **in parallel** (one message, five Agent/Task calls), each with this prompt
   (replace RUN with the absolute path):
   - `agente-omico`: "Question: <question>. Write your output to RUN/agent_omics.json."
   - `agente-literatura`: "Question: <question>. Write your output to RUN/agent_literature.json."
   - `agente-externo`: "Question: <question>. Write your output to RUN/agent_external.json."
   - `agente-celular`: "Question: <question>. Write your output to RUN/agent_cells.json."
   - `agente-aplicado`: "Question: <question>. Write your output to RUN/agent_applied.json."

3. Run `uv run python scripts/verify.py agents`. It removes evidence ids that do not exist or are not linked to
   the gene and drops unsupported findings, writing `*.verified.json`.

4. Read `RUN/convergence.json` (top 25) and the five `RUN/agent_*.verified.json` files. Integrate them: where do
   omics, literature and external knowledge converge, where do they disagree, which candidates are targets vs.
   biomarkers vs. mechanism. Then write `RUN/narrative.json`:

   ```json
   {"summary":  [{"text": "...", "evidence_ids": ["..."]}],
    "sections": {"omics": [...], "literature": [...], "external": [...], "cells": [...], "spatial": [...], "structure": [...], "discovery": [...], "validation": [...], "robustness": [...], "treatment": [...], "subtypes": [...], "adc": [...], "immune": [...], "integration": [...]},
    "targets":  [{"gene": "CA9", "category": "target|biomarker|target_and_biomarker",
                  "rationale": "...", "evidence_ids": ["..."]}],
    "limitations": [{"text": "..."}]}
   ```
   - `cells`, `spatial` and `structure` come from `agent_cells.verified.json`: which cell type expresses each
     candidate (tumor-intrinsic vs stromal / immune), what surrounds it in situ, and which candidates have a
     computationally detected cavity in a predicted monomer. Say explicitly when a bulk-omics candidate is not tumor-intrinsic, and
     present structural-context scores as descriptive heuristics, never as evidence of a degron or degrader. Use only curated experimental degron/ternary references and distinguish their original protein, ligase and compound from candidate-specific evidence.
   - `discovery` (2–4 paragraphs): from the `role: "discovery"` findings of `agent_cells.verified.json`. Treat each track
     separately (tumor, stroma, immune: the pool is mostly immune genes). For the leading candidates say what makes them
     novel and what the independent layers show — TCGA replication, whether the change survives extrapolation to a pure tumour
     (immune genes that vanish at purity 1 are infiltration-driven), survival and DepMap dependency — and keep the caveat
     that they are data-ranked hypotheses (no candidate is a selective dependency).
   - `validation` (2–3 paragraphs): what the independent layers say about the candidates — external replication
     in TCGA-KIRC and prognosis, somatic driver status, selective DepMap dependency, protein-adjusted kinase activity,
     and the purity check (does a bulk change survive extrapolation to a pure tumour or is it infiltration?). Name the
     candidates these layers strengthen and the ones they weaken.
   - `treatment`, `subtypes`, `adc` and `immune` (1-2 paragraphs each) come from `agent_applied.verified.json`. `treatment`: what
     IMmotion150 shows for the leading candidates, always as association with the adjusted p-value, and say that no gene
     separates the two main arms. `subtypes`: which candidates change with VHL, PBRM1, BAP1, SETD2 or ccA/ccB and which are
     subgroup candidates, with the infiltration caveat. `adc`: which candidates suit an antibody-drug conjugate or radioligand
     (surface, staining in vital organs, homogeneity) and which known targets act as controls. `immune`: which candidates are
     immune-cell genes, how they relate to cold and hot tumours, and the controls. Also use `confirmatory_replication` and
     `sc_replication_check` findings in `validation` and `cells`: say which protein calls replicate in the confirmatory cohort
     and which do not (HK2 is an example of a candidate that does not replicate there).
   - `robustness` (1–2 paragraphs): read `profiles.target` and `profiles.biomarker` and every row's
     `rank_interval_90` in `convergence.json`. Separate the candidates whose position survives random re-weighting
     (high `p_top10` / `p_top25`, narrow interval) from those that depend on the chosen weights.
   - **Unexplored candidates are part of the summary, not only of the `discovery` section.** The convergence ranking rewards
     literature and existing drugs, so it can never put a little-studied gene on top. The discovery screen is where those
     genes live. The summary must therefore contain one dedicated paragraph on the best unexplored candidates: two or three
     genes from `role: "discovery"` findings (prefer the tumour-cell track; say when the best ones belong to the stroma or
     immune track), and for each one say in plain sentences what supports it (cell source, replication in TCGA, protein or
     phospho evidence, tractability), what is still missing (accessible surface protein, staining in normal tissues,
     selective dependency, confirmatory-cohort data, heterogeneity across tumour cells) and which kind of study would test
     it. Frame them as "identified and prioritised for experimental validation", never as "discovered" or "proven":
     "novel" here means at most 15 ccRCC papers and no clinical drug, not that nothing was published, and the analysis does
     not show that treating the gene would help a patient. If one of them has earlier preclinical literature, say that the
     analysis reinforces an existing hypothesis; check it with the literature agent's findings and do not cite papers you
     did not see in the evidence.
   - Put at least two genes from the discovery screen in `targets`, ordered by your judgment, with `category` "target" and a
     `rationale` that begins with "Exploratory:" and names the main limitation.
   - Write the `discovery` section in prose, one paragraph per leading gene (not a list of scores): what the gene is, which
     cell makes it, what the independent layers show, what is missing. Put the numeric annotations at the end of each
     sentence, not in the middle of it.
   - 3–5 summary paragraphs; 1–3 paragraphs per section; 10–14 entries in `targets` ordered by your judgment of
     strength (you may depart from the numeric ranking, but say why in `rationale`). Use the target profile for
     target calls and the biomarker profile for biomarker calls; the overall convergence score mixes both.
   - Read `docs/numeric-claims.md`. Bind every number to its record with the `numeric_annotation` tool (kg-mcp): it
     returns the exact annotation to paste once into the text and the `numeric_claims` entry to add. Copying annotations
     from the verified agent files is fine when the paragraph cites the same evidence id. An annotation replaces its number; do not
     also print the number in plain text. Do not use ranking/global number
     pools. Preserve signs and units. Passing automatic checks is not complete scientific validation.
   - Write `narrative.json` in parts, not in one large write: first create it with `summary` and an empty `sections`
     object, `targets` list and `limitations` list; then add each section, the targets and the limitations with separate
     edits. A single very large write can be refused, and a partial file is easier to resume.
   - Every summary/section paragraph and every target MUST cite evidence ids copied from the verified agent files
     or from `convergence.json`. A paragraph that discusses one specific gene should carry `"gene": "<symbol>"`
     and cite only that gene's evidence. Limitations need no ids.

5. Run `uv run python scripts/assemble.py <RUN date>` then `uv run python scripts/verify.py final <RUN date>` (the date is
   RUN's folder name; without it the scripts use `runs/.current`).
   If it reports any unsupported claim, review and correct the original agent/narrative inputs (metric, citation, context or interpretation) and repeat step 5
   and re-run agent verification when its inputs change. Inspect all dropped statements; do not discard substantive results merely to obtain `unsupported claims: 0`.

6. Run `uv run python scripts/render.py <RUN date>` and reply with the report path and a 5-line summary.
