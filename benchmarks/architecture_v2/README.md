# Single agent vs multi-agent: architecture benchmark v2

**Question.** Does splitting the evidence interpretation across specialist agents, with or without an orchestrator,
give better answers than one agent reading everything, and is it worth its cost?

**Answer from this benchmark: no.** On integrated, masked evidence dossiers built from the atlas data, one agent was
right every time (all 7 answers, 144 of 144 dossier-answers across two load levels). The multi-agent arms were never more
accurate, cost 1.3 to 2.9 times as much per case, and lost answers to broken outputs, a failure mode that grows with the
number of calls. The orchestrator was useful only as a corrector of its own specialists.

## Design (decided before running; v1 in PR #10 was replaced because it could not discriminate)

- **24 dossiers** from 12 genes (CA9, CD70, NNMT, SLC16A3, ENPP3, HAVCR1, COL23A1, RAB42, PSAT1, PCK1, CD74, HK2):
  12 native, 12 counterfactual copies with one decisive record changed (`frozen/provenance.json`). Each dossier holds
  13–18 real records of the target over 11 layers (CPTAC RNA/protein, TCGA, confirmatory proteome, purity, single
  cell, surface, normal-tissue staining and proteome, DepMap, treatment, survival) and 12 look-alike records of two
  other masked entities. Gene names are masked.
- **Task:** six yes / no / insufficient sub-questions (protein up, replicated, tumour-intrinsic, surface-accessible,
  vital-organ safety concern, selective dependency) and a verdict (surface target candidate, loss biomarker,
  infiltration signal, needs more data). The definitions (`frozen/protocol.txt`) are given to every arm; the expected
  answers are those definitions applied to the target's records. Traps: transcript share vs share of expressing cells,
  missing vs negative, within-arm vs interaction, other entities' records. Verdict classes 7 / 7 / 6 / 4 (always
  guessing the commonest class: 29 %).
- **Arms** (same model `claude-sonnet-5-5`, effort medium, same 16,000 output-token cap and $2 cap per call, all tools,
  MCP servers, skills and settings disabled):
  - **B, single agent:** one context with the whole dossier answers everything.
  - **C, specialists only:** four contexts (expression; cells; safety; dependency and clinical), each seeing only its
    own record types and answering only its own questions; the verdict is the rule applied to their answers.
  - **D, specialists + orchestrator:** the same four specialist outputs as C, then an orchestrator with the whole
    dossier and the specialists' answers gives the final answers (the atlas's own pattern). C and D share the
    specialist calls; their cost is charged to both.
- **Load:** 4 dossiers per context (main, 3 replicates) and 12 dossiers per context, about 360 records (3 replicates).
  Two variants of one gene never share a context.
- **Failures:** a call whose output is not valid JSON for all its cases is retried once; both views are reported
  (first attempt, and after one retry with the retry's cost added). Missing answers count as wrong.

## Results

Main condition (4 dossiers per context, 3 replicates, 72 dossier-answers per arm):

| Arm | Verdict correct (first attempt) | Verdict correct (after one retry) | Correct among complete answers | Confident errors per dossier | Cost per dossier | Cost per correct verdict | Wall time per dossier | Failed calls (after retry) |
|---|---|---|---|---|---|---|---|---|
| B single agent | **100 %** | **100 %** | 72 / 72 | 0 | **$0.018** | **$0.018** | 5.1 s | 0 / 18 |
| C specialists | 54 % | 82 % | 59 / 60 | 0.03 | $0.026 | $0.032 | 4.7 s | 3 / 72 |
| D specialists + orchestrator | 89 % | 89 % | 64 / 64 | 0 | $0.052 | $0.059 | 11.0 s | 5 / 90 |

Load condition (12 dossiers, about 360 records, per context):

| Arm | Verdict correct (first attempt) | Verdict correct (after one retry) | Correct among complete answers | Confident errors per dossier | Cost per dossier | Cost per correct verdict |
|---|---|---|---|---|---|---|
| B single agent | **100 %** | **100 %** | 72 / 72 | 0 | **$0.015** | **$0.015** |
| C specialists | 61 % | 76 % | 55 / 60 | 0.08 | $0.020 | $0.026 |
| D specialists + orchestrator | 83 % | 83 % | 60 / 60 | 0 | $0.040 | $0.048 |

Paired differences in verdict accuracy (after one retry; bootstrap over the 12 genes, 95 % interval):
C − B −0.18 (−0.26 to −0.11) and D − B −0.11 (−0.17 to −0.06) in the main condition; C − B −0.24 (−0.35 to −0.17) and
D − B −0.17 at 12 dossiers per context. Full tables: `runs*/results.json`.

What the numbers say:

1. **No accuracy gain from splitting.** B made no error. Where C and D answered completely they matched B, never beat
   it, so even the best case for multi-agent is a tie at 1.3–2.9× the cost.
2. **Decomposition alone hurts.** Specialists that see only their slice made real reasoning errors (one in the main
   condition; five under load, with confident wrong yes/no answers), which a fixed combination rule passes straight
   into the verdict.
3. **The orchestrator earns its place only inside a multi-agent design.** D corrected every specialist error it saw
   (64/64 and 60/60 among complete answers), but it is the most expensive arm and a single extra point of failure: a
   broken orchestrator output loses all dossiers of that context (two such batches were still broken after the retry
   in the main condition, and one under load).
4. **More calls, more breakage.** All failures were malformed outputs (a stray `</br>` or `}` inside the JSON, several
   JSON objects, cases left out). B had none in 24 calls; C and D had 11 and 14 first-attempt failures in 96 and 120
   calls. Structured outputs or file writing with validation would remove most of them, but this is the operational
   price of more calls.
5. **Load did not change the picture** up to about 360 records per context: B stayed at 100 %, specialists got worse.

## Recommendation for the atlas

For the interpretation and integration step, a single agent with the whole evidence is enough and is cheaper,
faster and more reliable. Keep the multi-agent pipeline only if a benefit is shown where this benchmark did not look:

- **Retrieval, not interpretation.** Here the evidence was handed over. The production agents find their own
  evidence with dozens of MCP calls, where parallel retrieval and separate contexts could matter. The next experiment
  is a single agent with all MCP tools against the current five agents on the real pipeline (quality by the verifier
  and the sabotage test, tokens, time).
- **Larger or messier evidence** than 360 records per context, or questions without an explicit rubric.

## Limits

One model (Sonnet 5.5) at medium effort; 12 genes from one disease; definitions given explicitly, so this measures
integration under stated rules, not open-ended scientific judgment; B reached the ceiling, so the benchmark cannot show
how much better it might be than the others, only that they were not better; counterfactual dossiers are synthetic
variations, not biology; four dossiers per context are not independent observations (the bootstrap pairs by gene).
Model outputs cost about $9 API-equivalent in total (165 calls including retries). Nothing here validates any target.

## Reproduce

```bash
python -m benchmarks.architecture_v2.build            # only if frozen/ is absent (it refuses to overwrite)
python -m benchmarks.architecture_v2.run              # dry run: settings and call count
python -m benchmarks.architecture_v2.run --execute --claude <path to claude.exe> --reps 1 2 3
python -m benchmarks.architecture_v2.run --execute --claude <path> --reps 1 2 3 --batch-size 12
python -m benchmarks.architecture_v2.evaluate main     # or: first_attempt, load12, load12_first_attempt
```

Running `--execute` again retries only failed calls (successful responses with the same prompt are reused).
