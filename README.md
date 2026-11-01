# ccRCC Target Atlas: finding drug targets and biomarkers in kidney cancer with data and AI agents

**Live report:** https://rsolerortuno.github.io/ccrcc-multiomics-graphrag/

This project asks one question about clear cell renal cell carcinoma (ccRCC), the most common kidney cancer:

> *Which genes are the best candidates for new drugs or for tests that tell patients apart, judged by how many independent kinds of evidence agree?*

It is a worked example of how to build such a system, not a finished product or medical advice. Everything is open: the code, the data sources, the checks and the places where it fails.

## The idea in plain words

Imagine you are looking for the one weak spot of a very complicated machine, and you have a hundred different sensors. Any single sensor can be wrong. A spot that lights up on one sensor is a hunch. A spot that lights up on many unrelated sensors is worth a closer look.

Kidney cancer works the same way. Almost every ccRCC tumour has lost a gene called *VHL*, which makes the cell act as if it were short of oxygen. Drugs exist, but patients often stop responding. So researchers look for new things to aim at (targets) and for measurements that predict who will do well (biomarkers).

Here we collect many kinds of public data about hundreds of tumours, and for about 19,000 genes we ask the same questions again and again:

- Is the gene higher or lower in tumours than in healthy kidney from the same patient?
- Does a second group of patients show the same thing?
- Which cell makes it: a cancer cell, a blood-vessel cell or an immune cell?
- Does the protein sit on the cell surface where a drug could reach it?
- Is it also high in the heart, liver or lung, where a drug could do harm?
- Do cancer cells die when the gene is switched off?
- Is the gene linked to how well patients do on a drug?
- What do papers and drug databases already say?

A gene that scores well on many of these, in data that were collected separately, is a candidate. The result is a ranked list and an interactive report.

## What you get

| Output | What it is |
|---|---|
| [`reports/ccrcc_report.html`](reports/ccrcc_report.html) | One self-contained web page with 13 tabs. It also runs online at the live link above. |
| [`evals/results.md`](evals/results.md), [`evals/history.md`](evals/history.md) | How well the scores recover known biology, and a log of every change, including changes made after seeing data. |
| [`PROGRESS.md`](PROGRESS.md) | The build log. |
| `runs/` | The saved outputs of each AI-agent run, with their checks. |

The report is written for people who are not specialists. The Overview tab explains the question, the method and each tab. Every figure, table and counter has a round **?** button that explains what it shows, how to read it and what to keep in mind. The Tissue viewer lets you look at real tumour slices with the microscope picture underneath and click a spot to see its neighbours.

## How it works

![Architecture: omics, single-cell and spatial analysis, literature, external evidence, knowledge graph, MCP servers, Claude Code agents and HTML report](assets/ccrcc-architecture.png)

There are four stages. The first two are ordinary code. The third uses language models. The fourth checks the third.

1. **Measure (code only).** Python scripts run the statistics. They compare tumour with matched healthy kidney, look at single cells and tissue slices, test other cohorts and trials, and pull in drug and protein databases. No language model sees raw data or computes a number.
2. **Store with receipts (knowledge graph).** Every result is saved as a *claim* in a small database, together with an *evidence record* that holds the number, the dataset and the date. The database refuses a claim that has no evidence record, so nothing can be asserted without a source.
3. **Ask the agents (Claude Code).** Five AI agents, each with its own tools, read the graph and write short findings:
   - the **omics agent** interprets changes in RNA, protein and phosphorylation;
   - the **literature agent** reads what papers say;
   - the **external agent** looks at drugs and protein databases;
   - the **cell agent** says which cell makes each gene and what surrounds it in the tissue;
   - the **translational agent** says how a candidate would hold up in practice: treatment response, tumour subtypes, antibody-drug conjugates, the immune environment and replication.

   An orchestrator then picks the final candidates and writes the summary.
4. **Check every sentence.** Another script verifies each sentence the agents wrote. The evidence it cites must exist, must belong to the same gene, the numbers must match, an "up" must really be up, and the evidence must not have changed since the agent read it. Sentences that fail are dropped.

The agents only choose what to look at and explain it. The numbers, tables, figures and citations in the report come from code.

## The evidence, layer by layer

| Layer | What it is | Why it matters |
|---|---|---|
| Tumour vs normal | CPTAC: 110 tumours with matched normals, measured as RNA, protein and phosphorylation | The basic question: what changes in the tumour |
| Copy number | Gains and losses of DNA, and whether they move the gene's product | Loss of the 3p chromosome arm is the hallmark of ccRCC |
| Mutations | CPTAC and TCGA tumours, corrected for gene length | Finds drivers (VHL, PBRM1, BAP1, SETD2) |
| Second cohort | TCGA-KIRC RNA and survival; the confirmatory CPTAC cohort for protein and phosphosites | A result that appears in other patients is more believable |
| Cell-line dependency | DepMap CRISPR screens | Shows whether cancer cells need the gene, or whether it is only high |
| Kinase activity | Phosphorylation turned into inferred enzyme activity | Enzymes are drug targets |
| Tumour purity | How much of a sample is cancer cells | Stops infiltrating cells from faking a tumour signal |
| Single cells | 178,817 cells from Li 2022, repeated in Zhang 2021 | Says which cell type makes each gene |
| Tissue slices | 23 Visium sections, four with the microscope image | Says where in the tumour a gene is active and what is next to it |
| Treatment response | IMmotion150 and Choueiri 2016 trials | Links genes to who benefits from which drug |
| Subtypes | VHL, PBRM1, BAP1, SETD2 mutation and ccA/ccB | Finds targets for one genetic subgroup |
| ADC suitability | UniProt surface topology, Human Protein Atlas staining in vital organs, evenness across tumour cells | A target for a conjugate must be on tumour cells and not in vital organs |
| Immune environment | CIBERSORT and xCell estimates of immune cells | Immune-cell genes can be targets too |
| Literature and drugs | PubMed relations, Europe PMC counts, Open Targets, ChEMBL | Tells what is already known and what is new |
| Pockets | AlphaFold structures and cavity detection | A rough guide to where a small molecule might bind. It is only a heuristic |

## What we found

Numbers below come from the current run. The full detail, with caveats, is in the report.

- **The method recovers known biology.** HIF targets such as CA9, NDUFA4L2, ANGPTL4 and VEGFA are strongly up in tumours and replicate in a second cohort. Loss of the 3p arm is the most frequent change. The mutation test calls VHL, PBRM1, BAP1, SETD2, KDM5C, MTOR, PTEN and TP53.
- **Candidates.** CA9 is the one gene that stays in the top ten whatever the weights. Target calls include CD70, NNMT, ENPP3, SCARB1, SLC16A3 and HK2 (see the confirmatory result below). Biomarker calls include PSAT1, PCK1, SLC6A19 and ALDH6A1, mostly kidney genes that tumours lose.
- **A caution that matters.** None of the expression-based targets is a gene the cancer cells need to survive. They are handles for antibodies and radioligands, not weak spots.
- **The protein data replicate.** In the confirmatory cohort, 79 % of the protein changes found in the discovery cohort and 73 % of the phosphosite changes appear again in the same direction (the fold changes correlate at 0.89 and 0.83). The gaps are informative: HK2 rises in discovery and does not in the confirmatory cohort.
- **Treatment response.** Two published gene signatures reproduce in IMmotion150, which shows the trial data are usable. Individual genes are associated with progression in each arm, but no gene separates immunotherapy from sunitinib after correcting for testing 20,000 genes.
- **Subtypes.** PBRM1 and BAP1 mutations reshape the tumour far more than VHL. Patients with a ccB tumour are about 2.3 times as likely to die in TCGA.
- **Antibody and radioligand targets.** The known targets CD70, HAVCR1, CA9 and ENPP3 rank 2, 4, 13 and 39 among 2,017 tumour-up genes, so the score picks them out.
- **Immune environment.** PDCD1 ranks 10th in the immune-target profile, but only one of seven checkpoint controls reaches the top 100, so that profile lists candidates and does not give a validated ranking.
- **The agents add traceability.** A language model asked from memory gives the textbook list with nothing to cite. The pipeline gives candidates backed by data, and every sentence is checked.

## How we know it is not nonsense

- **Controls.** Wherever a score is new, known cases are scored the same way. Known ADC targets rank near the top. Known checkpoints rank unevenly, and the report says so.
- **Held-out clinical targets.** The target profile uses no literature or drug history, and still ranks real RCC drug targets at AUROC 0.90, and 0.887 in the retrospective post-2016 group (6 genes; 1 in the top 100). Trial dates include all phases; genes with incomplete dates are excluded from that group. Features include post-2016 information, so this is not a prospective prediction test.
- **Replication.** Protein and phosphosite results are repeated in a separate patient cohort with a different method. Single-cell attributions are repeated in a second dataset.
- **Stability.** Every rank has an interval from 1,000 random re-weightings.
- **Repeated agent runs (historical).** The original agent step was run three times, each from a clean folder with the tools reading that same folder. The final target lists overlapped with a mean Jaccard of 0.94 and the top four candidates (CA9, NNMT, CD70, SCARB1) were identical in all three. Those runs used the earlier verifier; their zero-failure counts do not establish compliance with the current metric-specific checks.
- **Checks on the text.** Numbers must match a specific cited metric, with sign, unit and source context preserved. Ambiguous multi-layer claims require explicit bindings; unrelated ranking/global values cannot rescue them. Rejected inputs remain visible as failures. A passing result means **no failures detected by the automatic checks**, not complete scientific validation. See [the numeric-claim contract](docs/numeric-claims.md). The current saved agent findings and integrated narrative have been reviewed against this contract; earlier archived runs retain their historical status.

## How the project is organised

```
config.yaml            thresholds and weights, set before looking at results
src/ccrcc/             all the code (details below)
scripts/               the commands the agent step runs
.claude/agents/        the five agents (plus the literature extractor)
.claude/commands/      the orchestrator command that runs them
data/                  downloaded data and results (not in git, rebuilt by the scripts)
runs/                  saved agent outputs and their verification, one folder per run
evals/                 evaluation results and the history of changes
reports/               the final HTML report
tests/                 automated tests
```

| Code in `src/ccrcc/` | What it does |
|---|---|
| `omics.py`, `layers.py`, `confirmatory.py` | Tumour-vs-normal statistics; mutations, copy-number effects, dependency, TCGA, phosphosites, purity; the confirmatory cohort |
| `treatment.py`, `subtypes.py`, `adc.py`, `immune.py` | Treatment response, subtypes, antibody-drug-conjugate suitability, immune environment |
| `singlecell.py`, `sc_replication.py`, `spatial.py`, `spatial_viewer.py` | Cell types, the second single-cell dataset, tissue neighbourhoods, data for the Tissue viewer |
| `structure.py`, `degradation.py`, `protein_viewer.py` | Descriptive pocket features, curated experimental degrons/ternary references and the coordinate-based 3D viewer |
| `discovery.py`, `safety_clinical.py`, `external.py` | The first-in-class screen, vital-organ expression, Open Targets, ChEMBL, Europe PMC |
| `literature.py`, `literature_entities.py` | PubMed ingestion and checked relation extraction |
| `graph.py`, `kg_build.py`, `kg_layers.py`, `kg_cell_structure.py`, `kg_applied.py` | The knowledge graph and the loaders that fill it |
| `scores.py`, `integrate.py` | The convergence score, the target and biomarker profiles, weight sensitivity |
| `verify.py`, `semantic.py`, `provenance.py`, `assemble.py` | The sentence checks and the stale-text guard |
| `mcp_servers/` | The tools the agents use (four MCP servers) |
| `outcomes.py`, `baselines.py`, `evaluate.py` | Tests against clinical targets and against a model that answers from memory |
| `render.py`, `report_*.py`, `report_assets/` | The HTML report: figures, help texts and templates |
| `extra_download.py`, `download.py` | Downloads, each recorded with a checksum |

## Run it yourself

You need [uv](https://docs.astral.sh/uv/) and, for the agent step only, Claude Code logged in with a subscription (no API key).

```bash
uv sync --frozen
uv run python -m ccrcc.reproducibility  # verify the retained input snapshot
uv run ccrcc run-all        # rebuild from frozen local inputs; network access is blocked
uv run pytest               # automated tests
uv run ccrcc eval           # evaluation against clinical targets, robustness, baselines
```

The frozen input inventory is in [reproducibility/inputs.json](reproducibility/inputs.json). Restore the listed data files from the retained local snapshot before running; upstream downloads can change and are not guaranteed to reproduce those bytes. [Reproduction notes](docs/reproduction.md) record coverage, source-version limitations and the successful clean execution. See also the [findings audit](docs/revision-audit.md) and [cleanup record](docs/cleanup.md). The original acquisition commands remain available separately. `run-all --allow-network` opts out of offline reconstruction.

The agent step regenerates the findings and the summary, and takes about half an hour:

```bash
./scripts/run_informe.sh                                   # Windows: ./scripts/run_informe.ps1
CCRCC_RUN_DIR=replicates/rep1 ./scripts/run_informe.sh     # a repeat run that does not touch the main one
uv run ccrcc baselines                                     # compare with a memory-only model and across repeat runs
```

Individual analyses can also run alone, for example `uv run ccrcc treatment`, `subtypes`, `adc`, `immune`, `confirmatory`, `sc-replication`, `spatial-viewer` or `discovery`. After any rebuild that changes evidence, run the agent step again: the stale-text guard rejects sentences written about old numbers.

The online report updates by itself. A workflow in `.github/workflows/pages.yml` publishes `reports/ccrcc_report.html` to GitHub Pages whenever that file changes on `main`.

## Decisions made after seeing data

These were not planned in advance. They are listed so you can judge them, and logged with their effect in `evals/history.md`.

- A selective dependency needs at least three dependent ccRCC cell lines out of 24.
- Literature counts use only gene names that belong to one gene. "Undrugged" is halved for genes with potent compounds or patents in ChEMBL.
- A clinical-precedence term was removed from the target profile because it leaked the answer into the evaluation.
- The number check was tightened after it let a wrong number through.
- Repeat agent runs now start from a clean folder and the agent tools read the same folder, after one repeat run reported it might have read the main run's files.

## Limits

- All data are public and from a few cohorts. Nothing here was tested in a laboratory or in patients.
- Bulk tissue mixes cell types. Single-cell data and purity corrections help, but the single-cell and tissue data come from other patients than the CPTAC cohorts and from few donors, and the dissociation step favours immune cells.
- The confirmatory proteome misses some genes (CA9 among them), because label-free data are sparse in normal kidney. A missing gene is a gap, not a failed replication.
- Cell-line dependency uses only 24 ccRCC lines grown in normal oxygen.
- Treatment-response tests are weak: one adequately sized trial, and no gene separates the two main arms.
- The structural-context score describes monomer pocket features; it is not a degradation predictor. Ligand-contact recovery tests cavity detection only. The 3D viewer displays grey protein surfaces, red selected residues, confidence, transparency and clipping. Experimental degrons and ternary references are separate annotations, not predictions transferred to new candidates.
- Gold sets of clinical targets are small, and clinical targets are mostly druggable by construction.
- Language-model steps cannot be repeated bit for bit. Their outputs are saved and checked.
- Not medical advice.
