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
- **Prioritized candidates (orchestrator):** CA9 (target and biomarker; strongest convergence across all streams), EPAS1, KDR and MTOR (established targets with genetic and clinical, not expression, support), CD70 and MET (expression plus a clinical modality), and metabolic tissue biomarkers NDUFA4L2/COXFA4L2, SLC16A3, PLIN2, ENO2 and ALDOB. The deterministic ranking's top genes (EGFR, CD74, FAS, CXCR4) are kept as hypotheses because of non-specific genetic evidence and possible immune infiltration.
- **Verification:** 59 narrative and agent statements checked against the graph, **0 unsupported**. The HTML validator found 84 claims and 282 evidence ids, all resolvable.
- **Evaluation:** recall of 23 known ccRCC genes in the top 100 is **22 %**, below the 60 % target, but it is 42× random enrichment (median rank 309/19,161; recall@1000 74 %). [`evals/history.md`](evals/history.md) explains why: known ccRCC genes are mostly *mutational* drivers with little bulk expression change, and mutation data were out of scope. Known-answer questions: 10/10. That result is interim, answered by a Sonnet 5.5 subagent through the same MCP tool functions. The `claude -p` version is pending (see below).

> **Status of the headless run.** The agent analysis in `runs/2026-10-31/` was orchestrated from an interactive
> Claude Code session (Opus 5.5 orchestrator, three Sonnet 5.5 subagents). It used the same agent definitions and
> called the same MCP tool functions through `scripts/mcp_call.py`. The fully headless command
> (`claude -p "/informe"`, `scripts/run_informe.ps1`) is ready but needs the Claude Code CLI to be logged in once.

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
| `src/ccrcc/mcp_servers/` | `omics-mcp`, `kg-mcp`, `external-mcp` (FastMCP) |
| `src/ccrcc/integrate.py` · `verify.py` · `assemble.py` · `render.py` | Convergence, verification, findings assembly, HTML |
| `.claude/agents/*.md` · `.claude/commands/informe.md` | Subagents (Sonnet 5.5) and orchestrator command |
| `config.yaml` | Thresholds and convergence weights (v2) |
| `evals/` | Known-answer questions, results, evaluation history |

## Moving to Neo4j

Only `GraphStore` touches storage. Implement the same methods (`add_node`, `add_evidence`, `add_claim`, `neighbors`,
`evidence_for`, `get_evidence`, `check_evidence_ids`, `resolve`, `stats`) over the Neo4j driver: nodes as labels,
Evidence as nodes, claims as `(:Entity)-[:SUBJECT_OF]->(:Claim {predicate})-[:OBJECT]->(:Entity)` with
`(:Claim)-[:SUPPORTED_BY]->(:Evidence)`. Agents, MCP servers and the report do not change.

## Adding single-cell and spatial data (next step)

Checked as downloadable: Li et al. 2022 (*Cancer Cell*) scRNA-seq + Visium h5ad on Mendeley Data (`g67bkbnhhg`);
GSE175540 (Visium, processed); GSE156632 (scRNA-seq, 7 tumor/normal pairs). Wu et al. 2023 snRNA/snATAC (CPTAC) is
dbGaP-controlled. Suggested integration: cell-type marker and proportion evidence as new Evidence types
(`CellType` nodes with Cell Ontology ids) and a cell-type specificity stream in the convergence score.

## Limitations

Single discovery cohort without external validation; bulk tissue (immune genes may reflect infiltration);
mutation data out of scope; literature corpus limited and biased toward omics-selected genes; Open Targets
evidence includes non-specific pan-cancer scores; LLM steps are not bit-reproducible (their outputs are saved and
verified). Not medical advice.
