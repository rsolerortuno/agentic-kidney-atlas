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

- Model: claude-sonnet-5-5
- **Correct: 10/10**

| ID | Question | Expected | Correct | Answer / error |
|---|---|---|---|---|
| q01 | Which chromosome arm is most frequently lost in the CPTAC ccRCC tumors? | 3p | ✅ | Chromosome arm 3p is lost most often: the VHL locus is lost in 79.1% of the 110 CPTAC tumors (EV:omics_cnv:52e5aeb6c197). The tool only reports one gene at a time, so I compared three genes rather than every arm. The oth |
| q02 | Is CA9 protein increased or decreased in ccRCC tumors compared with matched normal kidney in CPTAC? | increased | ✅ | CA9 protein is increased in ccRCC tumors versus matched normal kidney in CPTAC (log2FC +2.27, padj 9.9e-43, 84 pairs; evidence ID EV:omics_de:5afe18b4024b). |
| q03 | Which antibody targeting CA9 reached phase 3 according to Open Targets? | girentuximab | ✅ | Girentuximab, the anti-CA9 antibody, reached Phase 3 in Open Targets (evidence ID EV:open_targets_drug:646a0aec5b64). |
| q04 | Which gene has the highest Open Targets association score with clear cell renal carcinoma? | VHL | ✅ | VHL has the highest Open Targets association score with clear cell renal carcinoma, at 0.6466 (evidence ID EV:open_targets_assoc:730901c424b4). |
| q05 | Is EGLN3 measured in the CPTAC ccRCC proteome used in this project? | no | ✅ | No. The omics tool returns only an RNA result for EGLN3 (log2FC 4.41, padj 1.4e-42, evidence EV:omics_de:1a89557965e3) and no proteome result, so EGLN3 does not appear to be measured in the CPTAC ccRCC proteome used here |
| q06 | In what fraction of CPTAC ccRCC tumors is VHL affected by copy-number loss (to two decimals)? | 0.79 | ✅ | VHL is affected by copy-number loss in 0.79 of CPTAC ccRCC tumors (loss frequency 0.791 across 110 tumors, evidence ID EV:omics_cnv:52e5aeb6c197). |
| q07 | Which approved drug targets EPAS1 (HIF-2alpha) according to Open Targets? | belzutifan | ✅ | Open Targets lists belzutifan (CHEMBL4585668) as the approved drug targeting EPAS1 (HIF-2α), with max stage APPROVAL (evidence ID EV:open_targets_drug:199c475802ce). |
| q08 | Which gene ranks first in this project's deterministic convergence ranking? | EGFR | ✅ | EGFR ranks first in the deterministic convergence ranking, with a score of 3.87 across 5 evidence streams. |
| q09 | Is UMOD protein increased or decreased in ccRCC tumors versus matched normal kidney in CPTAC? | decreased | ✅ | UMOD protein is decreased in ccRCC tumors versus matched normal kidney in CPTAC (log2FC −2.24, padj 3.2e-33, 84 pairs; evidence EV:omics_de:e6e9e4ef8e40). |
| q10 | Name a high-confidence STRING interaction partner of CA9 that is a HIF alpha transcription factor. | HIF1A or EPAS1 | ✅ | HIF1A is a high-confidence STRING partner of CA9, with a combined score of 0.972 (evidence ID EV:string:9e9468e265e0). EPAS1, which encodes HIF-2α, also scores high at 0.923. |

## Analysis

Recall@100 is below the 60 % target. The cause is structural, not a bug: see [history.md](history.md) for the evaluated versions, the single methodological change, and why mutational drivers with unchanged bulk expression (VHL, PBRM1, SETD2, BAP1, KDM5C…) cannot reach the top 100 of an expression-convergence ranking without a somatic-mutation layer.
