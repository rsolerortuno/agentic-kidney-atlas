# Architecture comparison — protocol v1

**Status: prepared, not executed. Claude quota is unavailable; no model outcomes exist.**

## Scientific question and scope

Does specialist decomposition or centralized synthesis improve the interpretation of the same multi-omics evidence? This first pilot isolates interpretation with frozen evidence packets. It does **not** evaluate live MCP retrieval, novel target discovery, treatment efficacy, or the complete production agent workflow. Packet-level D reproduces its parallel-specialists/central-synthesis pattern, not its full prompts or report generation.

## Four conditions

- **A:** explicit three-valued domain rules applied to each packet. A is a strong rule baseline for this bounded benchmark, not the atlas's numerical target ranking. Because the tasks have explicit rules, A may reach the ceiling; that is expected and must not be presented as general superiority of deterministic systems.
- **B:** one generalist, all packet evidence, high reasoning effort and explicit self-review instruction.
- **C:** five specialist contexts, each reading all evidence with a different review lens. A fixed majority rule combines their labels; without three matching valid votes it abstains. Citations are the union of agreeing votes. This is an ensemble of independent reviews, not five mutually exclusive evidence partitions.
- **D:** exactly the same five worker outputs as C, followed by an LLM synthesis with the original packets available. This paired comparison isolates synthesis. D is charged the full worker cost even though workers are physically executed once and reused.

All model calls use the same configured Sonnet model and effort. Built-in tools, MCP servers, skills and repository access are disabled equally. Only model-visible packet JSON and instructions enter isolated contexts. Hidden labels and source provenance do not enter prompts.

## Cases and independence

The frozen suite has 12 seed cases in six families, with three variants each: **36 scenarios**. The first variant is a normalized projection of measured records in the local atlas graph; the other two are counterfactuals. Twelve native and 24 synthetic scenarios are reported separately. Synthetic positive treatment or degradation results are **not biological discoveries**. `provenance.json` records every source record, field projection and manipulation. Gene names are masked in model inputs; native and synthetic scenarios share neutral identifiers.

Families: independent protein replication, RNA/protein agreement, cellular origin, an explicitly defined surface-screen rule, treatment interaction, and experimental degradation evidence. Numerical thresholds are operational task definitions, not clinical recommendations. Examples deliberately distinguish significance from absence, transcript share from cell detection frequency, and pockets from experimental degradation.

Six cases run in a model context to make the subscription-based pilot practical. Variants of the same seed never share a context. Batch and record order are randomized deterministically across three replicates. There are **126 physical classification calls**: 18 generalist, 90 specialist, 18 synthesis. This produces 108 case answers per LLM architecture and 36 deterministic answers. Cases in a context are not independent replications.

These are small, selected cases from one disease, with repeated genes across families. There is no random sampling of all possible biological tasks. No prompts have been tuned on model results. Any later prompt changes require a new version and a new held-out suite; publishing these cases makes future reuse a public benchmark rather than a secret test set.

## Resource controls

Per six-case batch, B and D have a provider-reported API-equivalent ceiling of $3. Each worker is capped at $0.50; C allocates $2.50 and leaves $0.50 unused. D adds a $0.50 synthesis call. Output caps are B 12,000; five workers 2,000 each; synthesis 2,000. C and D's sharing enables a direct synthesis ablation but means C does not spend the whole common ceiling.

These are **common ceilings, not equal actual compute or tokens**. Claude Code controls provider budget per process; the final request may overshoot and cap behavior must be checked in the first smoke run. Reasoning-token accounting may not be exposed separately. Report actual costs, cached and uncached usage, elapsed times, failures and cap violations; do not claim exact compute matching. Subscription API-equivalent costs are not necessarily cash charges. The main classification execution has a maximum allocated estimate of $108 before overshoot, excluding judges. No automatic paid top-up is enabled by this runner.

## Primary outcomes

1. Exact structured decision accuracy (supported / contradicted / insufficient).
2. False confident conclusions, reported separately as critical classification errors.
3. Required-evidence recall and relevant-record citation precision. These test ID selection, **not semantic entailment of every prose sentence**.
4. Correct and unnecessary abstention, reported separately.
5. Both members of a native/counterfactual pair correct. Pair stability alone is not rewarded if both answers are wrong.
6. Parse, timeout and worker failures; model/provider usage and time; cost per correctly solved case.

The independent evaluator reads frozen labels; it never calls baseline A or the production numeric verifier. Metamorphic tests include sign reversal, withholding a measurement, and adding a significant but non-decisive within-arm association. Record-order randomization and unrelated-gene distractors occur in every arm. This first design does not isolate a separate clean-vs-distractor treatment effect; do not claim a measured distraction-resistance benefit.

## Analysis

Compare B−A, C−B, D−C and D−B on paired case means. Average variants and repeated outputs within each biological seed; use a seeded percentile bootstrap over the 12 seeds. Report per-family and native/synthetic results, not just a pooled number. The pilot intervals do not fully account for shared genes or batch dependence and are exploratory; do not use them for confirmatory p-values or superiority claims. Report ceiling effects and ties without claiming equivalence. No winner is chosen by overlap with the atlas target list.

Before a larger confirmatory study, use the pilot to estimate variance and define a practically meaningful improvement and affordable cost. A second experiment with exactly matched realized budgets and a held-out retrieval task set would be needed for stronger architectural claims.

## Secondary LLM assessment

The optional judge uses the 12 native cases, first replicate, all four arms, with anonymized identities and two reversed presentation orders. The rubric is accuracy, integration, contradictions, calibration and next experiment (0–2 each), with quoted justification and separate overclaim flags. The two passes use the same model: they are order-sensitivity checks, **not independent expert judges**. Architecture may remain inferable from response style, especially the vote-based C output. Scores never select the primary winner. No human-expert validation is claimed. Judge calls are a separate explicit opt-in and budget (24 calls, at most $1 each before overshoot).

## Commands

Run from the repository root in its Python environment:

```powershell
python -m pytest tests/test_architecture_benchmark.py -q
python -m benchmarks.architecture.run                 # dry-run only
python -m benchmarks.architecture.judge               # dry-run only
# Only after the user confirms Claude quota has recovered:
python -m benchmarks.architecture.run --execute --claude "C:/path/to/claude.exe"
python -m benchmarks.architecture.run --evaluate-only
python -m benchmarks.architecture.judge --execute --claude "C:/path/to/claude.exe"
```

The runner stores the full prompt, parsed response, usage and timing for each call. Successful outputs are resumed by protocol hash; quota stops are not retried automatically. Resolve a recorded quota-stop file deliberately before resuming. Timeout and parse failures remain in the record and count as failures rather than being silently resampled. Source inputs remain unchanged; the graph was opened read-only. Frozen hashes are checked before execution.

## Remaining work after quota recovery

Verify CLI budget/accounting on development-only smoke inputs; execute the classification pilot; inspect operational completeness; run the secondary judges; finalize cost/time and order-sensitivity summaries; generate an accessible HTML report and figures; publish the results and limitations in the repository and portfolio. Do not publish model score tables before actual executions exist.
