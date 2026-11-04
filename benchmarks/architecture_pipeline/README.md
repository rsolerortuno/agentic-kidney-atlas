# Single agent vs multi-agent on the real pipeline (with MCP tools)

The interpretation benchmark (`../architecture_v2`) handed the evidence to the agents. Here each arm must **find** the
evidence itself with the project's MCP tools and write the full analysis, as the production pipeline does.

- **Multi-agent (M):** `/informe` as published: an Opus 5.5 orchestrator launches the five Sonnet 5.5 specialists in
  parallel, verifies their findings and writes the narrative.
- **Single agent (S):** `/informe_single` (`.claude/commands/informe_single.md`), the same instructions except step 2:
  one Opus 5.5 agent reads the five specialist definitions and does all five jobs itself with the same MCP tools,
  writing one findings file. No subagents.
- Same question, data, graph, tools, verifier and narrative requirements; two runs per arm, interleaved
  (S, M, S, M); each run in its own folder. Cost, turns, tokens and wall time are Claude Code's own JSON result.
- Quality without a gold answer, fixed before the runs (`score.py`): automatic verification; coverage (evidence
  types, distinct evidence records, genes); a 13-item checklist of facts the data establish (e.g. CD70 partly made by
  T cells, HK2 not replicating, EPAS1 the only selective dependency, no gene separating treatment arms, use of the new
  proteome, CheckMate 025 and PRISM data); overlap with the published target list.

## Results

| Run | Statements checked / unsupported | Checklist (of 13) | Evidence types cited | Distinct evidence records | Cost | Wall time | Turns | Output tokens |
|---|---|---|---|---|---|---|---|---|
| M1 | 123 / 0 | 10 | 28 | 170 | $3.62 | 4.8 min | 40 | 30,095 |
| M2 | 112 / 0 | 11 | 30 | 172 | $3.59 | 5.0 min | 29 | 31,840 |
| S1 | 105 / 0 | 12 | 32 | 148 | $6.16 | 45.5 min | 256 | 96,795 |
| S2 | 118 / 0 | 11 | 33 | 148 | $5.75 | 44.9 min | 249 | 92,959 |
| **M mean** | 118 / 0 | 10.5 | 29 | 171 | **$3.60** | **4.9 min** | 35 | 31k |
| **S mean** | 112 / 0 | 11.5 | 32.5 | 148 | $5.95 | 45.2 min | 253 | 95k |

The three published multi-agent runs of 2026-11-03 (no cost data recorded) score 12, 11 and 11 on the checklist, so
the multi-agent arm's range is 10–12 and the single agent's 11–12. A first single-agent attempt was cut off by the
usage limit after 34 minutes and $2.52; it is not in the table.

What the numbers say:

1. **Quality is a tie.** Both arms pass verification with no unsupported statement and cover the same facts within
   one checklist item; with two runs per arm the difference is noise. The single agent cites slightly more evidence
   types, the multi-agent arm more distinct records.
2. **Multi-agent is cheaper and much faster.** About 40 % cheaper ($3.60 vs $5.95) and about 9 times faster
   (5 vs 45 minutes). The single agent needs some 250 sequential turns and three times the output tokens; the
   specialists run in parallel and on the cheaper model.
3. **Why:** the saving comes from delegating the tool-heavy retrieval to parallel workers on a cheaper model, and from
   each worker keeping a small context. Part of the cost gap is the model mix (all-Opus vs Opus + Sonnet); a
   single Sonnet agent was not tested. Long single sessions are also more exposed to usage limits, as the cut-off
   attempt showed.

## Conclusion, together with the interpretation benchmark

- **Interpreting evidence that is already gathered:** a single agent is better: as accurate or more, cheaper, faster
  and more reliable (`../architecture_v2`).
- **Gathering the evidence with tools and writing the analysis (this pipeline):** multi-agent is worth it, not for
  quality (a tie) but for cost (−40 %) and time (−89 %).
- For the atlas: keep the multi-agent pipeline; its value is parallel, cheaper retrieval, and the orchestrator's
  integration and correction step should stay with one capable model.

## Limits

Two runs per arm; one model pairing (Opus orchestrator / Sonnet specialists vs Opus alone); no gold answer, so
quality is measured by verification, coverage and a fact checklist, which cannot judge the depth of the reasoning;
the checklist is regex-based and can miss a fact stated in other words.

## Reproduce

`benchmarks/architecture_pipeline/run_arm.ps1 -Arm multi|single -Run 2026-11-04_<name>` runs one analysis through
`claude -p` (Opus 5.5, project MCP servers, JSON result; the multi-agent arm also allows the Agent/Task tools). Move
`runs/2026-11-04_<name>` to `runs/<name>` here, save its `session_<name>.json` there as `session.json` with an
`"arm"` field, and run `python -m benchmarks.architecture_pipeline.score`.
