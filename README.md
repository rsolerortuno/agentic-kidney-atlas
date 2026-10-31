# Building an agentic multi-omics GraphRAG to prioritize targets in ccRCC — a case study

This repository documents **how to build** an agentic system that combines multi-omics data, literature
and external knowledge in an evidence-reified knowledge graph, and uses LLM agents to identify candidate
therapeutic targets and biomarkers, and why they matter. The disease is clear cell renal cell carcinoma
(ccRCC). It is an analysis and a worked example, not a general-purpose tool: paths, thresholds and choices
are specific to this question and these data.

> **Question:** *What are the candidate therapeutic targets and biomarkers in ccRCC with the strongest
> convergent evidence?*

Outputs: [`reports/ccrcc_report.html`](reports/ccrcc_report.html) (self-contained, open in a browser) ·
[`evals/results.md`](evals/results.md) · [`evals/history.md`](evals/history.md) · [`PROGRESS.md`](PROGRESS.md)

## Design in one picture

```
 CPTAC ccRCC (LinkedOmics)      PubMed via PubTator3          Open Targets · STRING · UniProt
 RNA · proteome · phospho · CNV  normalized annotations        (disk-cached)
          │                          │                                 │
   M2 deterministic stats     M3 closed-schema relation           M5 external client
   paired DE, CNV freq,       extraction (Sonnet subagents)            │
   RNA–protein correlation    + verbatim validation                    │
          └──────────────┬───────────┴─────────────────────────────────┘
                         ▼
        M4 knowledge graph (SQLite + NetworkX behind GraphStore)
        every (subject)-[predicate]->(object) points to an Evidence record
                         │
        M5 MCP servers: omics-mcp · kg-mcp · external-mcp (FastMCP, stdio)
                         │
        M6 Claude Code: orchestrator (Opus 5.5) ──► agente-omico · agente-literatura · agente-externo (Sonnet 5.5)
                         │            scripts/integrate.py (deterministic convergence)
                         │            scripts/verify.py (drops claims without gene-linked evidence)
                         ▼
        M7 Jinja2 HTML report: LLM writes only the interpretive text; numbers, tables, figures, citations are code
```

## Key design decisions (and why)

| Decision | Rationale |
|---|---|
| Raw omics never enter the graph or a prompt | Statistics are computed once, deterministically (paired t-test + BH on matched cases). Agents query results through tools. |
| Evidence reification | `(gene)-[upregulated_in]->(ccRCC)` is stored only through an Evidence record (statistic, dataset, date). The SQLite schema enforces `evidence_id NOT NULL` + FK, so every sentence in the report can be traced. |
| PubTator3 instead of generic GraphRAG extraction | Entities come already normalized (NCBI Gene, MeSH). The LLM only extracts relations in a closed schema, and every relation must quote a sentence that is verified verbatim against the abstract. |
| Verification that is *gene-aware* | An agent claim about gene X is kept only if it cites evidence ids that exist **and** are attached to a claim involving X. A real id about another gene is rejected. |
| Convergence = independent evidence streams | Omics layers are combined into one stream so that differential expression alone (thousands of genes) cannot dominate. Weights live in `config.yaml`. |
| No paid API | All LLM work runs on Claude Code with a Pro subscription: subagents in `.claude/agents/`, orchestrator command `.claude/commands/informe.md`. |

## Data

- **CPTAC ccRCC discovery cohort** (LinkedOmics): 110 tumors; matched normals 75 (RNA), 84 (proteome, phosphoproteome); gene-level copy number. sha256 in `data/raw/manifest.json` (regenerated).
- **HGNC** complete set for identifier mapping (e.g. NDUFA4L2 is HGNC `COXFA4L2` today; the literature uses the old symbol).
- **PubTator3**: 319 abstracts (top 100 for the disease + 5 per top-50 omics gene), 7,270 normalized annotations; 100 abstracts used for relation extraction → 198 verbatim-validated relations.
- **Open Targets Platform** (MONDO_0005005): 1,998 associations with datatype scores; tractability, drugs and Reactome pathways for 371 genes. **STRING** (score ≥ 0.7). **UniProt** Swiss-Prot.

Raw CPTAC matrices and PubMed abstracts are **not redistributed** in this repository. The scripts download them again.

## Results (summary)

See the report for the full, cited analysis. Headlines:
- Expected ccRCC biology is recovered: HIF targets CA9, NDUFA4L2, ANGPTL4, EGLN3 and VEGFA are strongly up in RNA (log2FC 2.9–6.8, padj < 1e-32), and CA9, NDUFA4L2 and ANGPTL4 are also up in protein. 3p is the most frequently lost arm (median 80 % of tumors).
- The knowledge graph holds ~96k evidence-backed claims, with 0 claims lacking resolvable evidence.
- **Prioritized candidates (orchestrator):** CA9 (target and biomarker; the only gene where omics, literature and a late-stage clinical modality all agree), EPAS1/HIF-2α (belzutifan approved; promoted from rank 38), MET and KDR (tractable, approved inhibitors), EGFR (rank 1 by score but kept as a hypothesis: no ccRCC literature or approval), and biomarkers NDUFA4L2/COXFA4L2, ENO2 and PLIN2 (up) and PIGR and ACY1 (down; ACY1 lost with 3p). CD74 and CXCR4 likely reflect immune infiltration.
- **Verification:** 64 narrative and agent statements checked against the graph, **0 unsupported**. The HTML validator found 89 claims and 251 evidence ids, all resolvable.
- **Evaluation:** recall of 23 known ccRCC genes in the top 100 is **22 %**, below the 60 % target, but it is 42× random enrichment (median rank 309/19,161; recall@1000 74 %). [`evals/history.md`](evals/history.md) explains why: known ccRCC genes are mostly *mutational* drivers with little bulk expression change, and mutation data were out of scope. Known-answer questions answered by `claude -p` (Sonnet 5.5) with the MCP servers: **10/10**.

> The agent analysis in `runs/2026-10-31/` was produced headless with `claude -p "/informe"` (Opus 5.5 orchestrator,
> three Sonnet 5.5 subagents, project MCP servers); log in `runs/2026-10-31/informe_headless.log`. An earlier
> session-orchestrated run is kept in `runs/archive/` for comparison.

## Reproduce

Requirements: Windows/macOS/Linux, [uv](https://docs.astral.sh/uv/), Claude Code logged in with a Claude subscription (for the agent step only).

```bash
uv sync
uv run ccrcc run-all          # data → omics → literature tables → enrichment → graph → convergence → verify → HTML
uv run pytest                 # tests (external calls are served from the disk cache after the first run)
```

Agent step (Claude Code, no API key). It regenerates `runs/<date>/` with fresh agent analyses and narrative:

```powershell
./scripts/run_informe.ps1     # = claude -p "/informe" --model claude-opus-5-5 --mcp-config .mcp.json --allowedTools ...
uv run ccrcc eval --qa        # recall + 10 known-answer questions answered by claude -p with the MCP servers
```

Relation extraction (M3) is an LLM step done once with the `extractor-relaciones` subagent. Its validated output
(`data/literature/extracted/*.jsonl`) is versioned and re-validated on every run.

## Repository map

| Path | What it is |
|---|---|
| `src/ccrcc/omics.py` | Paired tumor-vs-normal DE, CNV frequencies, RNA–protein correlation, HGNC mapping |
| `src/ccrcc/literature.py` | PubTator3 ingestion, extraction batches, verbatim validation |
| `src/ccrcc/graph.py` / `kg_build.py` | `GraphStore` (SQLite + NetworkX) and graph construction |
| `src/ccrcc/mcp_servers/` | `omics-mcp`, `kg-mcp`, `external-mcp`, `cell-mcp` (FastMCP) |
| `src/ccrcc/singlecell.py` · `spatial.py` · `structure.py` · `kg_cell_structure.py` · `cell_context.py` | M11 cell-type, spatial and MGD layers |
| `src/ccrcc/integrate.py` · `verify.py` · `assemble.py` · `render.py` | Convergence, verification, findings assembly, HTML |
| `.claude/agents/*.md` · `.claude/commands/informe.md` | Subagents (Sonnet 5.5) and orchestrator command |
| `config.yaml` | Thresholds and convergence weights (v2) |
| `evals/` | Known-answer questions, results, evaluation history |

## Moving to Neo4j

Only `GraphStore` touches storage. Implement the same methods (`add_node`, `add_evidence`, `add_claim`, `neighbors`,
`evidence_for`, `get_evidence`, `check_evidence_ids`, `resolve`, `stats`) over the Neo4j driver: nodes as labels,
Evidence as nodes, claims as `(:Entity)-[:SUBJECT_OF]->(:Claim {predicate})-[:OBJECT]->(:Entity)` with
`(:Claim)-[:SUPPORTED_BY]->(:Evidence)`. Agents, MCP servers and the report do not change.

## Cell-type, spatial and structural context (M11)

Added after the first full run, as annotation layers: they do **not** enter the convergence score (so recall@100 is unchanged
by construction, see `evals/history.md`). Data: Li et al. 2022 (*Cancer Cell*) via CZ CELLxGENE Discover, 178,817 scRNA-seq
cells (12 donors) and 38,836 Visium spots (23 sections); AlphaFold DB models; PDBe-KB ligand sites.

| Question | Method (`src/ccrcc/`) | Output |
|---|---|---|
| Which cell type expresses the target? Is the bulk signal tumor-intrinsic? | `singlecell.py`: donor x tissue x cell-group UMI pseudobulks streamed from the 1.6 GB h5ad; share of the gene's tumour transcripts per group; tumour cells vs normal proximal tubule, donor-paired t-test + BH | `sc_*.parquet` |
| Where is it expressed in situ and what surrounds it, weighted by abundance? | `spatial.py`: NNLS deconvolution with the study's own profiles (targets excluded from markers); NNLS of target CPM on cell shares; ring-1 neighbours of top-decile spots, observed/expected vs section-matched random spots, 1000 permutations, BH | `spatial_*.parquet` |
| Where could a molecular glue degrader bind? | `structure.py`: AlphaFold model trimmed at pLDDT 50, pyKVFinder cavities, MGD site score = confidence, pocket geometry, solvent-accessible lysines within 20 A (acceptor for E2~Ub), exposed non-helical loop (degron shape) | `structure_*` |

What the layers say (details and evidence ids in the report):

- **CA9** is tumour-intrinsic (97.6 % of its tumour transcripts come from tumour cells; log2FC 7.7 vs proximal tubule). **EGFR and MET** too.
- **KDR** is endothelial (95 % of its tumour transcripts); **CD74, CXCR4, PTPRC** are immune. Candidates ranked high by bulk expression that are
  not tumour-cell targets are now identified from data instead of suspected.
- **EPAS1** (HIF-2a) transcripts in this dataset are mostly endothelial and fibroblast (5 % tumour cells): interpret it as a pathway target (belzutifan), not as a tumour-cell expression marker.
- Spatial deconvolution is checked with held-out genes: CA9 vs tumour share rho 0.28, KDR vs endothelium 0.22, CD74 vs macrophage 0.33, positive but moderate.

**MGD scores are a prioritization of sites, not evidence that a degrader exists.** Honest validation against 23 proteins with experimental
ligand complexes: the cavity detector covers a median 44 % of ligand-bound residues; the score places the best-matching cavity first
in 11 of 23, but the top 3 pockets cover a median of only 10 %. A first version had a saturated lysine component (80 % of pockets at 1.0),
fixed from the component's distribution before interpreting any gene. No protein-specific E3 compatibility, ternary complex or
degron-sequence modelling is done; glue-degrader chemistry would need docking / cryo-EM follow-up.

New MCP server `cell-mcp` (celltype_expression, tumor_cell_vs_proximal_tubule, spatial_context, structure_pockets, mgd_ranking) and a
fourth subagent, `agente-celular` (Sonnet 5.5). Not used: Wu et al. 2023 snRNA/snATAC (dbGaP-controlled); GSE175540 / GSE156632 / Mendeley
`g67bkbnhhg` were checked as alternatives, but the CELLxGENE-curated copy of Li 2022 has harmonized Cell Ontology labels.

## Limitations

Single discovery cohort without external validation; bulk tissue (immune genes may reflect infiltration, now quantified by scRNA-seq);
scRNA-seq / Visium come from a different cohort than CPTAC (12 donors; dissociation bias favours immune cells; Visium spots hold several cells);
MGD sites are heuristic and weakly validated (see above);
mutation data out of scope; literature corpus limited and biased toward omics-selected genes; Open Targets
evidence includes non-specific pan-cancer scores; LLM steps are not bit-reproducible (their outputs are saved and
verified). Not medical advice.
