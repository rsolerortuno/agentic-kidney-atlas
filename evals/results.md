# Evaluation results

Convergence **v3** (M13). Earlier versions and why each change was made: [history.md](history.md).

## (a) Recall of known ccRCC genes in the top 100 convergence genes

- Known genes: 23 (curated 10 + Open Targets top 20, deduplicated)
- **Recall: 39%** (v2: 22 %; target 60 %)
- Recall with the genetics stream removed: 4%
- recall@50/100/200/500/1000 = 26% / 39% / 52% / 70% / 78%; median rank 171 of 19,161 genes; random expectation 0.52% → 75× enrichment

| Gene | Source | Rank | In top | In top without genetics stream |
|---|---|---|---|---|
| VHL | curated | 16 | ✅ | ❌ |
| PBRM1 | curated | 57 | ✅ | ❌ |
| SETD2 | curated | 49 | ✅ | ❌ |
| BAP1 | curated | 125 | ❌ | ❌ |
| KDM5C | curated | 171 | ❌ | ❌ |
| EPAS1 | curated | 38 | ✅ | ❌ |
| CA9 | curated | 1 | ✅ | ✅ |
| VEGFA | curated | 126 | ❌ | ❌ |
| MTOR | curated | 56 | ✅ | ❌ |
| HIF1A | curated | 2451 | ❌ | ❌ |
| AKT1 | open_targets_top20 | 225 | ❌ | ❌ |
| HNF1A | open_targets_top20 | 497 | ❌ | ❌ |
| TP53 | open_targets_top20 | 40 | ✅ | ❌ |
| HNF1B | open_targets_top20 | 1243 | ❌ | ❌ |
| PTEN | open_targets_top20 | 62 | ✅ | ❌ |
| KDR | open_targets_top20 | 690 | ❌ | ❌ |
| SPDL1 | open_targets_top20 | 364 | ❌ | ❌ |
| FLT4 | open_targets_top20 | 2097 | ❌ | ❌ |
| PIGR | open_targets_top20 | 6 | ✅ | ❌ |
| TERT | open_targets_top20 | 386 | ❌ | ❌ |
| FLCN | open_targets_top20 | 7519 | ❌ | ❌ |
| KMT2C | open_targets_top20 | 5637 | ❌ | ❌ |
| ATM | open_targets_top20 | 663 | ❌ | ❌ |

## (b) Outcomes the scores do not see (M13)

The target and biomarker profiles use measurements only (no literature, no drug history), so they can be scored against clinical outcomes. The overall convergence score includes literature and drug history and is shown only to make that circularity visible. Gold-set construction: `src/ccrcc/outcomes.py`.

### RCC clinical targets (75 genes: phase ≥ 2 in renal cell carcinoma, ≥ 2 trials or approval, drugs with ≤ 5 targets)

| Score | Gold genes | AUROC (95 % CI) | Recall@100 | Recall@500 | Median rank | Enrichment@100 | Gold genes in top 100 |
|---|---|---|---|---|---|---|---|
| target profile | 75 | 0.896 (0.87–0.92) | 17% | 24% | 1,673 | 33.2× | CA9, IDO1, CCND1, ANGPT2, FLT1, ITGAV, EPAS1, EGFR, ITGA5, VEGFA, TOP2A, MET, CDK6 |
| biomarker profile | 75 | 0.705 (0.65–0.76) | 3% | 12% | 4,986 | 5.1× | CA9, TOP2A |
| convergence v3 (incl. literature + drugs) | 75 | 0.890 (0.88–0.91) | 5% | 17% | 2,341 | 10.2× | CA9, EPAS1, MTOR, IDO1 |
| convergence v3 without literature/tractability | 75 | 0.700 (0.64–0.76) | 1% | 12% | 4,288 | 2.6× | CA9 |
| target profile without tractability | 75 | 0.718 (0.65–0.78) | 11% | 20% | 4,025 | 20.4× | CA9, IDO1, CCND1, ANGPT2, FLT1, ITGAV, EPAS1, EGFR |
| predicted tractability alone | 75 | 0.851 (0.83–0.87) | 4% | 8% | 2,656 | 7.7× | NR3C1, TIGIT, IL2RB |

Predicted tractability (strong SM / AB buckets) flags 5,301 genes. It alone gives most of the AUROC — clinical targets are, almost by definition, druggable proteins — so the data layers should be judged by early enrichment: the share of gold genes in the top 100, which within the tractable set would be ~1.9% per tractable gold gene by chance, and by the "without tractability" row.

### Temporal holdout: first RCC trial on or after 2016-01-01 (13 genes: CCND1, CDK4, CDK6, EPAS1, FOLH1, GLS, IDO1, LAG3, PARP1, PARP2, PARP3, TIGIT, TNFSF11)

Earlier clinical targets are removed from the background.

| Score | Gold genes | AUROC (95 % CI) | Recall@100 | Recall@500 | Median rank | Enrichment@100 | Gold genes in top 100 |
|---|---|---|---|---|---|---|---|
| target profile | 13 | 0.908 (0.87–0.94) | 31% | 31% | 1,905 | 58.8× | IDO1, CCND1, EPAS1, CDK6 |
| biomarker profile | 13 | 0.654 (0.52–0.77) | 0% | 15% | 7,011 | 0.0× | — |
| convergence v3 (incl. literature + drugs) | 13 | 0.895 (0.86–0.94) | 15% | 31% | 2,477 | 29.4× | EPAS1, IDO1 |
| convergence v3 without literature/tractability | 13 | 0.662 (0.51–0.82) | 0% | 31% | 6,580 | 0.0× | — |
| target profile without tractability | 13 | 0.713 (0.57–0.83) | 23% | 31% | 5,566 | 44.1× | IDO1, CCND1, EPAS1 |
| predicted tractability alone | 13 | 0.850 (0.81–0.89) | 8% | 8% | 3,206 | 14.7× | TIGIT |

### Prognostic biomarkers: ClearCode34 (33 of 34 genes mapped; Brooks et al. 2014)

ClearCode34 was validated in TCGA, which the biomarker profile's prognosis component also uses; the last row removes that component. A literature-based gold set (verbatim `biomarker_of` relations from the disease-level query) was empty (0 genes): all such relations came from gene-specific queries chosen from the omics ranking, so they were not used.

| Score | Gold genes | AUROC (95 % CI) | Recall@100 | Recall@500 | Median rank | Enrichment@100 | Gold genes in top 100 |
|---|---|---|---|---|---|---|---|
| target profile | 33 | 0.713 (0.62–0.79) | 6% | 18% | 5,315 | 11.6× | EPAS1, TLR3 |
| biomarker profile | 33 | 0.752 (0.68–0.82) | 0% | 3% | 3,632 | 0.0× | — |
| convergence v3 (incl. literature + drugs) | 33 | 0.703 (0.62–0.79) | 3% | 6% | 3,863 | 5.8× | EPAS1 |
| convergence v3 without literature/tractability | 33 | 0.645 (0.55–0.75) | 0% | 6% | 5,031 | 0.0× | — |
| biomarker profile without prognosis | 33 | 0.705 (0.63–0.79) | 0% | 6% | 4,969 | 0.0× | — |

## (c) Robustness to the weights

1,000 random re-weightings (Dirichlet, κ = 20, around the configured weights). Genes in the top 10 in ≥ 80 % of draws:

- **convergence**: CA9
- **target**: SCARB1, CA9, SLC16A3, CD70, NNMT, ENPP3, FABP6
- **biomarker**: PSAT1, NNMT, ALDH6A1

## (d) What the agentic layer adds

| | Closed-book LLM (Sonnet 5.5, no tools) | GraphRAG pipeline |
|---|---|---|
| Runs | 3 | 3 |
| Mean pairwise Jaccard of candidate lists | 0.667 | 0.944 |
| Candidates that are RCC clinical targets | 60%, 60%, 60% | 18%, 18%, 17% |
| Candidates in the data-only target or biomarker top 500 | 80%, 70%, 80% | 100%, 100%, 100% |
| Statements citing checkable evidence | 0 | 91 checked, 0 unsupported |

Closed-book candidates: VHL, EPAS1, VEGFA, KDR, CD274, PBRM1, BAP1, SETD2, CA9, MET; VHL, EPAS1, VEGFA, PDCD1, CD274, PBRM1, BAP1, SETD2, MET, MTOR; EPAS1, VHL, VEGFA, KDR, PDCD1, MTOR, PBRM1, BAP1, SETD2, CA9.

Pipeline final targets: CA9, CD70, NNMT, SCARB1, PSAT1, ENPP3, SLC16A3, ALDH6A1, ALDOB, SLC22A8, EPAS1; CA9, CD70, NNMT, SCARB1, PSAT1, ENPP3, SLC16A3, ALDH6A1, ALDOB, SLC22A8, EPAS1; CA9, CD70, NNMT, ENPP3, SCARB1, EPAS1, SLC16A3, PSAT1, ALDH6A1, ALDOB, ENO2, SLC22A8.

## (e) Known-answer questions (`claude -p` + project MCP servers)

- Model: claude-sonnet-5-5
- **Correct: 10/10**

| ID | Question | Expected | Correct | Answer / error |
|---|---|---|---|---|
| q01 | Which chromosome arm is most frequently lost in the CPTAC ccRCC tumors? | 3p | ✅ | Chromosome arm 3p is the most frequently lost in the CPTAC ccRCC tumors. VHL is lost in 79.1% and PBRM1 in 80.0% of the 110 tumors (evidence IDs EV:omics_cnv:52e5aeb6c197 and EV:omics_cnv:07b6205e67e5).  The copy-number  |
| q02 | Is CA9 protein increased or decreased in ccRCC tumors compared with matched normal kidney in CPTAC? | increased | ✅ | CA9 protein is increased in ccRCC tumors versus matched normal kidney in CPTAC (log2FC +2.27, padj 9.9e-43, 84 pairs; evidence EV:omics_de:5afe18b4024b). |
| q03 | Which antibody targeting CA9 reached phase 3 according to Open Targets? | girentuximab | ✅ | Open Targets lists **girentuximab** (CHEMBL1743027) as the CA9-targeting antibody that reached phase 3 (evidence ID EV:open_targets_drug:646a0aec5b64). |
| q04 | Which gene has the highest Open Targets association score with clear cell renal carcinoma? | VHL | ✅ | VHL has the highest Open Targets association score with clear cell renal carcinoma, at 0.6466 (evidence ID EV:open_targets_assoc:730901c424b4). |
| q05 | Is EGLN3 measured in the CPTAC ccRCC proteome used in this project? | no | ✅ | No. The project's CPTAC ccRCC results for EGLN3 include only an RNA entry (log2FC 4.41, padj 1.4e-42). There is no proteome entry for it. |
| q06 | In what fraction of CPTAC ccRCC tumors is VHL affected by copy-number loss (to two decimals)? | 0.79 | ✅ | VHL has copy-number loss in 0.79 of CPTAC ccRCC tumors (loss frequency 0.791 across 110 tumors, on chromosome arm 3p; evidence ID EV:omics_cnv:52e5aeb6c197). |
| q07 | Which approved drug targets EPAS1 (HIF-2alpha) according to Open Targets? | belzutifan | ✅ | According to Open Targets, the approved drug targeting EPAS1 (HIF-2alpha) is belzutifan (CHEMBL4585668, max stage APPROVAL, evidence ID EV:open_targets_drug:199c475802ce). |
| q08 | Which gene ranks first in this project's deterministic convergence ranking? | CA9 | ✅ | **CA9** ranks first in the project's deterministic convergence ranking, with a score of 4.52 across 5 evidence streams. |
| q09 | Is UMOD protein increased or decreased in ccRCC tumors versus matched normal kidney in CPTAC? | decreased | ✅ | UMOD protein is decreased in ccRCC tumors versus matched normal kidney in CPTAC (log2FC −2.235, padj 3.18e-33, 84 pairs; evidence ID EV:omics_de:e6e9e4ef8e40). |
| q10 | Name a high-confidence STRING interaction partner of CA9 that is a HIF alpha transcription factor. | HIF1A or EPAS1 | ✅ | HIF1A is a high-confidence STRING partner of CA9 (score 0.972, evidence ID EV:string:9e9468e265e0). EPAS1, which encodes HIF-2α, also scores high at 0.923. |
