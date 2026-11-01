# Building an agentic multi-omics GraphRAG to prioritize targets in ccRCC: a case study

This repository shows how to build a system that combines multi-omics data, literature and external knowledge in an
evidence-backed knowledge graph, and uses LLM agents to pick candidate therapeutic targets and biomarkers for clear cell
renal cell carcinoma (ccRCC). It is a worked example, not a general tool: paths, thresholds and choices belong to this
question and these data.

> **Question:** *What are the candidate therapeutic targets and biomarkers in ccRCC with the strongest convergent evidence?*

Live report: https://rsolerortuno.github.io/ccrcc-multiomics-graphrag/

Outputs: [`reports/ccrcc_report.html`](reports/ccrcc_report.html) (one self-contained file, 13 tabs, a "?" button on every
figure, table and counter) · [`evals/results.md`](evals/results.md) · [`evals/history.md`](evals/history.md) · [`PROGRESS.md`](PROGRESS.md)

## Design in one picture

Core architecture (the treatment, subtype, ADC, immune and second single-cell analyses sit in the independent-layers box). Code computes the statistics, rankings and report assets; Claude Code
extracts literature relations and writes interpretations grounded in the processed results and their evidence.

![ccRCC agentic multi-omics GraphRAG architecture: omics, single-cell and Visium spatial analysis, literature, external evidence, knowledge graph, MCP servers, Claude Code agents and HTML report](assets/ccrcc-architecture.png)

**How to read the diagram**

- **The independent layers combine complementary measurements with shared cohorts.** Somatic recurrence uses CPTAC + TCGA
  (cBioPortal; coding-length correction); CNV cis effects reuse CPTAC; dependency comes from DepMap;
  replication and univariate overall-survival analysis use TCGA-KIRC (Xena); protein-adjusted kinase activity
  uses CPTAC phosphosites + OmniPath; purity effects use CPTAC ESTIMATE. These layers feed the
  [configured convergence and profile components](config.yaml). The target and biomarker profiles exclude
  literature and drug history; convergence v3 still includes a literature stream.
- **Cell, spatial and structure results contextualize candidates; discovery has its own ranking.** cell and spatial context
  and the discovery novelty scores do not enter convergence v3. Discovery combines bulk and scRNA evidence with
  tractability, clinical drugs, off-tumour expression, grade/stage, TCGA replication/prognosis, purity and
  dependency, and expands the spatial/structure shortlist. Spatial analysis uses fixed marker exclusions and
  per-gene seeded permutations. MGD scores prioritize sites; they do not demonstrate degradation.
- **Four analysis agents share four MCP servers.** `agente-omico`, `agente-literatura`, `agente-externo` and
  `agente-celular` run in parallel under `/informe`; `extractor-relaciones` is the separate extraction
  subagent. Agent findings are checked before Opus writes the narrative, and the assembled findings are checked
  again before rendering. `basis.json` fingerprints evidence so changed values invalidate stale text.
- **Rebuild and LLM regeneration are separate operations.** `ccrcc run-all` rebuilds deterministic artifacts
  and revalidates saved extraction/agent outputs; `scripts/run_informe.sh` (or `.ps1`) regenerates the four-agent
  analysis and narrative. relation extraction is a separate saved LLM step. Verification checks citation
  linkage, numeric/direction consistency and evidence freshness; it does not establish biological causality.

Implementation: [`cli.py`](src/ccrcc/cli.py) · [`discovery.py`](src/ccrcc/discovery.py) ·
[`integrate.py`](src/ccrcc/integrate.py) / [`scores.py`](src/ccrcc/scores.py) ·
[MCP configuration](.mcp.json) · [`/informe`](.claude/commands/informe.md) ·
[`verify.py`](src/ccrcc/verify.py) / [`provenance.py`](src/ccrcc/provenance.py) ·
[`render.py`](src/ccrcc/render.py)

## Key design decisions

| Decision | Rationale |
|---|---|
| Raw omics never enter the graph or a prompt | Statistics are computed once, deterministically. Agents query results through tools. |
| Evidence reification | Every claim is stored through an Evidence record (statistic, dataset, date). The SQLite schema enforces `evidence_id NOT NULL` with a foreign key. |
| PubTator3 instead of generic GraphRAG extraction | Entities come normalized. The LLM only extracts relations in a closed schema, and each must quote a sentence verified word for word. |
| Verification is gene-aware and number-aware | A claim about gene X must cite evidence attached to X, and every number and the up or down wording must match that evidence. |
| Convergence means independent evidence streams | Omics layers form one stream so expression alone cannot dominate. Weights are in `config.yaml`, set before looking at results. |
| Two questions, two scores | "Is it a target?" and "is it a biomarker?" need different evidence. Neither profile uses literature or drug history, so both can be tested against clinical outcomes. |
| Robustness is reported | Every rank comes with a 5th to 95th percentile interval from 1,000 random re-weightings. |
| Stale narrative guard | `basis.json` fingerprints the evidence the agents saw. `verify.py final` fails any sentence whose evidence changed afterwards. |
| No paid API | All LLM work runs on Claude Code with a subscription: subagents in `.claude/agents/`, orchestrator `.claude/commands/informe.md`. |

## Data

- **CPTAC ccRCC discovery cohort** (LinkedOmics): 110 tumours, matched normals 75 (RNA) and 84 (proteome, phospho), gene-level copy number. sha256 in `data/raw/manifest.json`.
- **HGNC** for identifier mapping. **PubTator3**: 319 abstracts, 100 used for relation extraction, 198 verbatim-validated relations.
- **Open Targets** (API and 26.09 bulk tables), **ChEMBL** (potent compounds, patent records), **STRING**, **UniProt**, **Human Protein Atlas** (RNA and normal-tissue staining), **AlphaFold DB**, **PDBe-KB**, Europe PMC counts.
- **Independent layers** (sha256 in `data/raw/extra/manifest.json`): somatic mutations (cBioPortal; CPTAC and TCGA-KIRC), DepMap 26Q1, TCGA-KIRC RNA and survival (UCSC Xena), CPTAC phosphosites, ESTIMATE purity, CIBERSORT and xCell (cptac Zenodo record), OmniPath.
- **Treatment response:** IMmotion150 (263 patients, atezolizumab with or without bevacizumab, or sunitinib) and Choueiri 2016 (16 patients, nivolumab), from the cBioPortal datahub.
- **Single cell and space:** Li et al. 2022 (178,817 cells, 12 donors; 38,836 Visium spots, 23 sections; four sections with H&E for the viewer) and Zhang et al. 2021 (20,509 ccRCC cells, 7 donors, plus 5 normal kidneys), both via CZ CELLxGENE.

Raw matrices and abstracts are not redistributed. The scripts download them again.

## What the analysis finds

- **Biology is recovered.** HIF targets (CA9, NDUFA4L2, ANGPTL4, EGLN3, VEGFA) are strongly up and replicate in TCGA-KIRC. Somatic recurrence calls VHL, PBRM1, BAP1, SETD2, KDM5C, MTOR, PTEN and TP53.
- **Candidates:** CA9 is the only gene robustly in the convergence top 10. The target calls include CD70, NNMT, ENPP3, SCARB1, SLC16A3 and HK2. The biomarker calls include PSAT1, PCK1, SLC6A19 and ALDH6A1. EPAS1 is kept on clinical grounds although its RNA barely moves.
- **Main caveat:** no expression-ranked target is a selective DepMap dependency. They are handles for antibodies, conjugates and radioligands, not genes the tumour needs.
- **Verification:** 97 statements in the final findings, 0 unsupported, 0 stale. The report validator checks 429 claims and 1,570 evidence ids, all resolvable.
- **Evaluation** (`evals/results.md`): recall of 23 known genes in the top 100 is 39 % (median rank 171). The target profile ranks RCC clinical targets at AUROC 0.90 and a temporal holdout (first trial 2016 or later) at 0.91. Predicted tractability alone gives 0.85, so the gain is early enrichment. The biomarker profile gets AUROC 0.74 on ClearCode34 (0.75 before the ccA/ccB component, within the confidence interval).
- **What the agents add:** a closed-book LLM gives the textbook list (60 % clinical targets, nothing citable). The pipeline gives data-supported candidates with every statement checked. Three agent runs agree with a mean Jaccard of 0.84 (the closed-book runs: 0.67). The top three candidates differ between runs, and one replicate reported that its agents' tools may have read the main run's files, so treat that figure as an upper bound on independence.
- **Cells and tissue:** CA9, EGFR and MET are made by tumour cells. KDR is endothelial. CD74, CXCR4 and PTPRC are immune. EPAS1 transcripts come mostly from endothelium and fibroblasts. In the held-out check of the spatial deconvolution, CA9 follows the tumour-cell share (rho 0.277), KDR endothelium (0.205), CD74 macrophages (0.304) and PTPRC T cells (0.123). A second single-cell dataset agrees: the dominant source matches in 59 of 63 tumour, 33 of 34 stroma and 183 of 313 immune discovery genes.
- **Discovery:** 410 genes are up in both RNA and protein and tied to one cell type. 204 are novel and 58 are first-in-class opportunities (immune 52, tumour 3, stroma 3). Literature counts use only aliases unique to the gene, and 45 genes have potent preclinical compounds in ChEMBL.
- **Treatment response:** two published signatures reproduce in IMmotion150, so the trial data work. 90, 103 and 13 genes associate with progression in the combination, sunitinib and atezolizumab arms, but none differs between the combination and sunitinib after correction.
- **Subtypes:** PBRM1 and BAP1 reshape the tumour far more than VHL. ccB tumours have a hazard ratio of 2.30 for death in TCGA. PBRM1-mutant tumours raise CA9 and ENPP3, and BAP1-mutant tumours raise HAVCR1 and MET.
- **ADC and radioligands:** the known targets CD70, HAVCR1, CA9, ENPP3 and FOLH1 rank 2, 4, 13, 39 and 488 of 2,017 tumour-up genes.
- **Immune:** 185 CPTAC tumours split into cold, intermediate and hot. PDCD1 ranks 10 in the immune-target profile, but only 1 of 7 checkpoint controls reaches the top 100, so it lists candidates and is not a validated ranking.

## The report

13 tabs: Overview (with what is new), Candidates, Omics, Cells and space, Structure and MGD, Discovery, Tissue viewer, Treatment and subtypes, ADC and immune, Validation and robustness, Literature and drugs, Network, Methods and evidence. The page is light and matte, with no dark mode. Every figure, table and counter has a "?" button that explains what it shows, how to read it and what to keep in mind. The **Tissue viewer** draws four tumour sections on their H&E image: pick a gene, its top-decile spots are ringed, click a spot to see its six neighbours and their cell makeup, and compare the neighbours of all high spots with the section average. The permutation p-values shown next to it come from the Python analysis over all 23 sections.

## Reproduce

Requirements: Windows, macOS or Linux, [uv](https://docs.astral.sh/uv/), and Claude Code logged in with a Claude subscription (agent step only).

```bash
uv sync
uv run ccrcc run-all          # everything deterministic, from downloads to the HTML (about 40 minutes the first time)
uv run pytest                 # 95 tests
uv run ccrcc eval             # recall, outcome gold sets, temporal holdout, robustness, baselines -> evals/results.md
```

Agent step (no API key). It regenerates `runs/<date>/` with fresh analyses and narrative:

```bash
./scripts/run_informe.sh      # Linux / macOS (Windows: ./scripts/run_informe.ps1)
CCRCC_RUN_DIR=replicates/rep1 ./scripts/run_informe.sh   # a replicate that does not touch runs/<date>/
uv run ccrcc baselines        # closed-book LLM runs, agent decisions vs ranking, replicate agreement
```

Single analyses can be run alone: `ccrcc treatment | subtypes | adc | immune | sc-replication | spatial-viewer | discovery`.
Run `/informe` again after any `run-all` that changes evidence, because the stale guard rejects old sentences.

## Repository map

| Path | What it is |
|---|---|
| `src/ccrcc/omics.py`, `literature.py`, `graph.py`, `kg_build.py` | Paired DE, CNV, RNA-protein; PubTator3 and relation validation; the graph |
| `src/ccrcc/layers.py`, `kg_layers.py`, `extra_download.py` | Mutations, cis dosage, DepMap (one `selective_dependency` rule), TCGA, phosphosites, purity; their downloads and graph loading |
| `src/ccrcc/treatment.py`, `subtypes.py`, `adc.py`, `immune.py`, `sc_replication.py` | Treatment response, genetic and ccA/ccB subtypes, ADC suitability, immune profile, second single-cell dataset |
| `src/ccrcc/singlecell.py`, `spatial.py`, `spatial_viewer.py`, `structure.py`, `cell_context.py` | Cell types, tissue neighbourhoods, viewer data, pocket analysis |
| `src/ccrcc/discovery.py`, `safety_clinical.py`, `external.py` | Target discovery, vital-organ expression and grade, Open Targets, ChEMBL, Europe PMC |
| `src/ccrcc/scores.py`, `integrate.py`, `semantic.py`, `provenance.py`, `verify.py`, `assemble.py` | Scores, profiles, sensitivity, verification, staleness guard |
| `src/ccrcc/outcomes.py`, `baselines.py`, `evaluate.py` | Evaluation against clinical outcomes and agent baselines |
| `src/ccrcc/render.py`, `report_figures*.py`, `report_help.py`, `m15_tables.py`, `report_assets/` | The HTML report, its figures, "?" texts and templates |
| `.claude/agents/*.md`, `.claude/commands/informe.md` | Subagents and the orchestrator command |
| `config.yaml` | Thresholds and weights |

## Moving to Neo4j

Only `GraphStore` touches storage. Implement the same methods over the Neo4j driver (nodes as labels, Evidence as nodes,
claims as `(:Entity)-[:SUBJECT_OF]->(:Claim {predicate})-[:OBJECT]->(:Entity)` with `(:Claim)-[:SUPPORTED_BY]->(:Evidence)`).
Agents, MCP servers and the report do not change.

## Decisions made after seeing data

These are not pre-registered. They are logged with their effect in `evals/history.md`.
- A selective dependency needs at least 3 dependent ccRCC lines (24 lines in total).
- Discovery literature counts use only aliases unique to the gene. "Undrugged" is halved for genes with potent ChEMBL compounds or patent records.
- A clinical-precedence term was removed from the target profile because it leaked the answer into the evaluation.
- The numeric check was tightened after it passed a wrong number.

## Limitations

External replication of protein and phosphosites is missing: the confirmatory CPTAC cohort sits in the Proteomic Data Commons, which returned HTTP 503 when this was built. TCGA replicates RNA only. Bulk tissue mixes cells (single-cell data and purity extrapolation help, and the latter assumes a straight-line effect). Only 24 DepMap ccRCC lines, grown in normal oxygen. Single-cell and Visium data come from other cohorts than CPTAC, with few donors and a dissociation bias towards immune cells. Treatment-response tests are weak: one adequately sized trial and no gene that separates the two arms. Pocket-based scores are a heuristic for ordering sites and are only weakly validated. Gold sets are small and clinical targets are mostly druggable by construction. The agents do not yet narrate the treatment, subtype, ADC and immune tabs. LLM steps are not bit-reproducible; their outputs are saved and checked. Not medical advice.
