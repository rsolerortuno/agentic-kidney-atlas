# Evaluation results

## (a) Recall of known ccRCC genes in the top 100 convergence genes

- Known genes: 23 (curated 10 + Open Targets top 20, deduplicated)
- **Recall: 22%**
- Recall with the Open Targets genetics stream removed (circularity check): 4%
- Supplementary: recall@50/100/200/500/1000 = 17% / 22% / 35% / 52% / 74%; median rank 309 of 19,161 genes; random expectation 0.52% → 42× enrichment

| Gene | Source | Rank | In top | In top without genetics stream |
|---|---|---|---|---|
| VHL | curated | 67 | ✅ | ❌ |
| PBRM1 | curated | 188 | ❌ | ❌ |
| SETD2 | curated | 205 | ❌ | ❌ |
| BAP1 | curated | 203 | ❌ | ❌ |
| KDM5C | curated | 546 | ❌ | ❌ |
| EPAS1 | curated | 38 | ✅ | ❌ |
| CA9 | curated | 10 | ✅ | ✅ |
| VEGFA | curated | 1355 | ❌ | ❌ |
| MTOR | curated | 301 | ❌ | ❌ |
| HIF1A | curated | 1406 | ❌ | ❌ |
| AKT1 | open_targets_top20 | 309 | ❌ | ❌ |
| HNF1A | open_targets_top20 | 814 | ❌ | ❌ |
| TP53 | open_targets_top20 | 199 | ❌ | ❌ |
| HNF1B | open_targets_top20 | 7806 | ❌ | ❌ |
| PTEN | open_targets_top20 | 581 | ❌ | ❌ |
| KDR | open_targets_top20 | 21 | ✅ | ❌ |
| SPDL1 | open_targets_top20 | 107 | ❌ | ❌ |
| FLT4 | open_targets_top20 | 702 | ❌ | ❌ |
| PIGR | open_targets_top20 | 5 | ✅ | ❌ |
| TERT | open_targets_top20 | 516 | ❌ | ❌ |
| FLCN | open_targets_top20 | 1475 | ❌ | ❌ |
| KMT2C | open_targets_top20 | 1425 | ❌ | ❌ |
| ATM | open_targets_top20 | 1519 | ❌ | ❌ |

## (b) Known-answer questions (`claude -p` + project MCP servers)

- Model: interim: Claude Sonnet 5.5 subagent via MCP tool functions (claude -p pending login)
- **Correct: 10/10**

| ID | Question | Expected | Correct | Answer / error |
|---|---|---|---|---|
| q01 | Which chromosome arm is most frequently lost in the CPTAC ccRCC tumors? | 3p | ✅ | Chromosome arm 3p is the most frequently lost (about 79-81% of tumors for VHL, PBRM1, SETD2, BAP1). |
| q02 | Is CA9 protein increased or decreased in ccRCC tumors compared with matched normal kidney in CPTAC? | increased | ✅ | CA9 protein is increased in tumors versus matched normal (log2FC +2.27, padj 9.9e-43). |
| q03 | Which antibody targeting CA9 reached phase 3 according to Open Targets? | girentuximab | ✅ | Girentuximab (CHEMBL1743027) reached phase 3. |
| q04 | Which gene has the highest Open Targets association score with clear cell renal carcinoma? | VHL | ✅ | VHL has the highest Open Targets association score (0.6466). |
| q05 | Is EGLN3 measured in the CPTAC ccRCC proteome used in this project? | no | ✅ | No, EGLN3 is not measured in the proteome (only RNA data exist for it). |
| q06 | In what fraction of CPTAC ccRCC tumors is VHL affected by copy-number loss (to two decimals)? | 0.79 | ✅ | 0.79 (VHL loss frequency 0.791, 3p, n=110 tumors). |
| q07 | Which approved drug targets EPAS1 (HIF-2alpha) according to Open Targets? | belzutifan | ✅ | Belzutifan (CHEMBL4585668) is the approved drug targeting EPAS1. |
| q08 | Which gene ranks first in this project's deterministic convergence ranking? | EGFR | ✅ | EGFR ranks first (score 3.87). |
| q09 | Is UMOD protein increased or decreased in ccRCC tumors versus matched normal kidney in CPTAC? | decreased | ✅ | UMOD protein is decreased in tumors versus matched normal (log2FC -2.24, padj 3.2e-33). |
| q10 | Name a high-confidence STRING interaction partner of CA9 that is a HIF alpha transcription factor. | HIF1A or EPAS1 | ✅ | EPAS1 (HIF-2alpha; STRING score 0.923); HIF1A (0.972) is also a partner. |

## Analysis

Recall@100 is below the 60 % target. The cause is structural, not a bug: see [history.md](history.md) for the evaluated versions, the single methodological change, and why mutational drivers with unchanged bulk expression (VHL, PBRM1, SETD2, BAP1, KDM5C…) cannot reach the top 100 of an expression-convergence ranking without a somatic-mutation layer.
