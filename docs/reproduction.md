# Frozen reproduction

`uv run ccrcc run-all` defaults to offline reconstruction. Before deleting any
outputs it verifies the expected bytes and hashes in `reproducibility/inputs.json`.
It skips acquisition commands and rejects HTTP requests through the pipeline's
requests transport. A missing response fails instead of fetching a newer one.
The original source snapshot must be retained locally; it is not uploaded here.

## Run the published analysis

1. Install the locked environment with `uv sync --frozen` (tested: Python 3.12.11,
   Windows; the precise environment is in `reproducibility/acceptance.json`).
2. Restore the frozen local `data/raw`, `data/cache`, literature PubTator documents,
   extracted relations and PMID selection listed in the manifest. Keep the
   repository's saved model answers, reference catalogs and selected narrative.
3. Run `uv run python -m ccrcc.reproducibility` to check the snapshot.
4. Run `uv run ccrcc run-all`. The default fresh run rebuilds processed tables,
   the graph, evaluations, verified findings and `reports/ccrcc_report.html`.

The manifest checks 6,827 data files, 18 project inputs, configuration and uv.lock.
Exact input-tree inventories reject unexpected cached files. The selected analysis
is pinned to `runs/2026-11-01`; a newer local run cannot silently replace it.
Additional replicate findings and closed-book responses are also rejected.

Original AlphaFold coordinates are inputs; `*.plddt50.pdb` files are regenerated
and excluded from the source snapshot. Source URLs, request keys, versions and
retrieval dates are retained where known. Provider releases and historical signed
URLs missing from acquisition records are explicitly marked as unknown. Fresh
upstream downloads are therefore not guaranteed to reproduce the frozen bytes.
`--allow-network` opts into acquisition behavior, not frozen reproduction; changed
inputs require a separate review rather than rewriting expectations during replay.

Saved model analyses and reviewed annotations are inputs. Reproduction does not
regenerate LLM responses. Evaluation regrades frozen QA answers and recomputes
recall, retrospective clinical outcomes, sensitivity and baseline comparisons.
Missing saved answers fail instead of starting a model. QA regex checks and the
statement verifier are not complete scientific validation.

## Verified clean execution

An isolated copy started with code and only manifested data inputs, without the
processed directory or graph. Dependency installation used `uv sync --frozen`.
The final full `uv run ccrcc run-all` completed with exit code 0, including the
recalculation of evaluation and HTML validation. Its artifacts were adopted as
the reviewed report. The original user checkout remained unchanged.

Evidence committed alongside this document:

- `reproducibility/clean-run.log`: full final execution log.
- `reproducibility/clean-comparison.json`: 43 Parquet tables matched exactly,
  including float values and row/column order. Two added structural alias columns
  equal their original legacy columns; no score weights or values changed.
- `reproducibility/acceptance.json`: environment, verification counts and hashes.
- `reproducibility/runtime-hashes.json`: hashes of all 77 runtime source/assets;
  these matched the source actually used in the isolated execution.

The evaluation JSON also matched exactly. All 77 findings were retained. Final
verification checked 130 statements with zero unsupported, semantic or stale
failures. The HTML validator found 462 claims, 1,730 distinct evidence IDs, 1,823
citation links and zero external resource dependencies. These are automatic
checks, not proof of every scientific interpretation.

The first clean attempt exposed a real dependency-order bug: enrichment must
create the gene universe before subtype summaries. The corrected order was
verified by two successful full rebuilds, with the last including final evaluation.
Expected input hashes were not relaxed to obtain these results.

Targeted `.gitattributes` rules preserve frozen text bytes. All 20 hashed project,
configuration and lock files matched after checkout from the staged Git tree
with both `core.autocrlf=false` and `core.autocrlf=true`. Existing CRLF blobs are
preserved byte-for-byte to avoid rewriting historical inputs.

## Re-freeze of the saved agent analysis (2026-11-02)

The agent step was run again after the orchestrator instructions changed (the summary now has to discuss the best
unexplored candidates from the discovery screen, and at least two of them go into the final target list). This is a
reviewed change to saved model analyses, not a refresh of upstream data:

- The selected analysis is now `runs/2026-11-02`. `runs/2026-11-01` stays in the repository as the earlier analysis.
- `reproducibility/inputs.json` was updated only for the entries that describe the saved analysis: the selected run's agent
  files, basis and narrative, and the two replicate findings. Raw data, caches, configuration and the dependency lock keep
  their frozen hashes.
- Final verification of the new run checked 113 statements with no failures detected by the automatic checks (0 unsupported,
  0 semantic, 0 stale). The two repeat runs checked 115 and 118. The final target lists of the three runs have a mean
  Jaccard overlap of 0.905, and COL23A1 and RAB42 appear in all three.
- The files in `reproducibility/clean-run.log`, `clean-comparison.json`, `acceptance.json` and `runtime-hashes.json` describe
  the earlier clean execution. They were not regenerated, because processed data, configuration and code that produce the
  tables did not change. The report templates, help texts and figures did change after that execution.
- A rebuild of the derived data on the developer machine (`run-all` without the frozen check, because three upstream PDB
  reference files changed) produced the same counts as before: 462 claims and 1,730 evidence records in the validator for the
  earlier analysis.
