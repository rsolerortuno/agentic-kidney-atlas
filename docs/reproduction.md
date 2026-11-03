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

## Re-freeze after the HGNC symbol-mapping fix (2026-11-03)

A code review found that three layers were keyed by their source's own gene symbols while the gene universe uses the
current HGNC symbol: TCGA-KIRC (Xena's 2016 symbols), CPTAC cis dosage and the purity table (approved symbols only).
Renamed genes therefore lost those layers when joined. NDUFA4L2 (now COXFA4L2) had TCGA log2FC 6.60 (padj 7.6e-31) in
`layer_tcga.parquet` under its old name, yet scored 0 for replication. 2,158 TCGA rows were affected, 219 of them for genes
strongly changed in CPTAC RNA. DepMap, the treatment cohorts and the subtype tables had the same problem on a smaller
scale. All of them now go through `omics.to_hgnc()`, with the same rules as the CPTAC symbol map (approved symbol first,
then a unique previous or alias symbol), and `tests/test_omics.py` checks the keys.

Effect on the deterministic results:

- Convergence ranking: COXFA4L2 moves from 17th to 2nd and HILPDA enters at 13th. CA9 and COXFA4L2 are the two genes in
  the top ten in at least 80 % of re-weightings (p_top10 0.889 and 0.878).
- Target profile: the top ten are the same genes; COXFA4L2 enters the top 25 (from 221st). Biomarker profile: G6PC1 moves
  from 1,642nd to 12th.
- Discovery: 62 first-in-class opportunities instead of 58 (CALHM6, COLGALT1, MILR1, MTFP1); HILPDA enters the tumour-cell
  track at 7th. COL23A1 and RAB42 still lead it.
- Evaluation: AUROC against clinical targets 0.895 (was 0.896). Ties are now broken at random rather than in convergence
  order, which moved "predicted tractability alone" from 0.847 to 0.857 and removed IL2RB from its top 100.

The same review also changed the checks: the stale-text guard fails closed (a cited id missing from `basis.json` fails;
the ranking ids are now in the basis from the start), the direction check also reads summary and section sentences that
name exactly one cited gene, the T-effector control needs p < 0.05, and genes without Human Protein Atlas data are no longer
treated as safe in the discovery screen (no gene was affected in this run).

The agent step was run again with `scripts/run_informe.ps1` into `runs/2026-11-03`. The headless orchestrator stopped once
at step 4 (its own safety check refused the large `narrative.json` write) and was resumed from step 4 without re-running
the integration or the subagents. Final verification checked 116 statements with no failures detected by the automatic
checks (0 unsupported, 0 semantic, 0 stale). The orchestrator rewrote 48 agent findings that cited numbers without the
structured `numeric_claims` binding, and added one HK2 finding built from its confirmatory-replication record.

- `reproducibility/inputs.json` now selects `runs/2026-11-03`. Its agent files, basis and narrative replace the 2026-11-02
  entries; `runs/2026-11-02` stays in the repository. Fifteen cached Open Targets, STRING and UniProt responses for genes
  that entered the ranking (e.g. G6PC1, ACP3, ADGRF1) were added as inputs.
- The repeat runs in `runs/replicates/` were not re-run; they belong to the 2026-11-02 analysis.
- This machine's data snapshot differs from the manifest in four files that predate this change (`raw/extra/manifest.json`
  and three PDB references in `cache/ternary`), so the rebuild ran the `run-all` steps without the frozen-input check. The
  files in `clean-run.log`, `clean-comparison.json`, `acceptance.json` and `runtime-hashes.json` still describe the earlier
  clean execution.
