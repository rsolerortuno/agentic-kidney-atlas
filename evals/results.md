# Evaluation results

Convergence **v3** (M13). Earlier versions and why each change was made: [history.md](history.md).

## (a) Recall of known ccRCC genes in the top 100 convergence genes

- Known genes: 23 (curated 10 + Open Targets top 20, deduplicated)
- **Recall: 39%** (v2: 22 %; target 60 %)
- Recall with the genetics stream removed: 4%
- recall@50/100/200/500/1000 = 26% / 39% / 52% / 65% / 78%; median rank 176 of 19,161 genes; random expectation 0.52% → 75× enrichment

| Gene | Source | Rank | In top | In top without genetics stream |
|---|---|---|---|---|
| VHL | curated | 18 | ✅ | ❌ |
| PBRM1 | curated | 58 | ✅ | ❌ |
| SETD2 | curated | 50 | ✅ | ❌ |
| BAP1 | curated | 126 | ❌ | ❌ |
| KDM5C | curated | 176 | ❌ | ❌ |
| EPAS1 | curated | 39 | ✅ | ❌ |
| CA9 | curated | 1 | ✅ | ✅ |
| VEGFA | curated | 127 | ❌ | ❌ |
| MTOR | curated | 57 | ✅ | ❌ |
| HIF1A | curated | 2560 | ❌ | ❌ |
| AKT1 | open_targets_top20 | 231 | ❌ | ❌ |
| HNF1A | open_targets_top20 | 509 | ❌ | ❌ |
| TP53 | open_targets_top20 | 41 | ✅ | ❌ |
| HNF1B | open_targets_top20 | 1287 | ❌ | ❌ |
| PTEN | open_targets_top20 | 63 | ✅ | ❌ |
| KDR | open_targets_top20 | 709 | ❌ | ❌ |
| SPDL1 | open_targets_top20 | 306 | ❌ | ❌ |
| FLT4 | open_targets_top20 | 2190 | ❌ | ❌ |
| PIGR | open_targets_top20 | 7 | ✅ | ❌ |
| TERT | open_targets_top20 | 395 | ❌ | ❌ |
| FLCN | open_targets_top20 | 7606 | ❌ | ❌ |
| KMT2C | open_targets_top20 | 5745 | ❌ | ❌ |
| ATM | open_targets_top20 | 680 | ❌ | ❌ |

## (b) Outcomes the scores do not see (M13)

The target and biomarker profiles use measurements only (no literature, no drug history), so they can be scored against clinical outcomes. The overall convergence score includes literature and drug history and is shown only to make that circularity visible. Gold-set construction: `src/ccrcc/outcomes.py`.

### RCC clinical targets (75 genes: phase ≥ 2 in renal cell carcinoma, ≥ 2 trials or approval, drugs with ≤ 5 targets)

| Score | Gold genes | AUROC (95 % CI) | Recall@100 | Recall@500 | Median rank | Enrichment@100 | Gold genes in top 100 |
|---|---|---|---|---|---|---|---|
| target profile | 75 | 0.895 (0.87–0.92) | 17% | 24% | 1,688 | 33.2× | CA9, IDO1, CCND1, ANGPT2, FLT1, ITGAV, EPAS1, EGFR, ITGA5, VEGFA, TOP2A, MET, CDK6 |
| biomarker profile | 75 | 0.685 (0.63–0.73) | 4% | 11% | 5,296 | 7.7× | CA9, TOP2A, FLT1 |
| convergence v3 (incl. literature + drugs) | 75 | 0.886 (0.87–0.90) | 5% | 17% | 2,442 | 10.2× | CA9, EPAS1, MTOR, IDO1 |
| convergence v3 without literature/tractability | 75 | 0.696 (0.64–0.75) | 1% | 12% | 4,416 | 2.6× | CA9 |
| target profile without tractability | 75 | 0.717 (0.65–0.78) | 11% | 20% | 4,089 | 20.4× | CA9, IDO1, CCND1, ANGPT2, FLT1, ITGAV, EPAS1, EGFR |
| predicted tractability alone | 75 | 0.857 (0.83–0.88) | 0% | 13% | 2,660 | 0.0× | — |
| target profile, tractability without ligand history | 75 | 0.797 (0.74–0.85) | 15% | 28% | 1,273 | 28.1× | CA9, IDO1, ANGPT2, FLT1, ITGAV, EGFR, ITGA5, VEGFA, TOP2A, MET, CDK6 |
| predicted tractability without ligand history, alone | 75 | 0.774 (0.72–0.82) | 0% | 12% | 2,533 | 0.0× | — |

Predicted tractability (strong SM / AB buckets) flags 5,301 genes. It alone gives most of the AUROC — clinical targets are, almost by definition, druggable proteins — so the data layers should be judged by early enrichment: the share of gold genes in the top 100, which within the tractable set would be ~1.9% per tractable gold gene by chance, and by the "without tractability" row.

The strong small-molecule bucket includes Open Targets' "High-Quality Ligand", which comes from ChEMBL medicinal chemistry: drug-discovery history, not a property of the tumour. Without it (pocket and location buckets only, 3,968 genes) the target profile's AUROC is 0.797 instead of 0.895. Much of the gap between the data-only rows and the full profile is therefore known chemistry, and the headline AUROC should be read with that in mind.

### Retrospective post-2016 evaluation: first RCC trial on or after 2016-01-01 (6 genes: CCND1, PARP1, PARP2, PARP3, TIGIT, TNFSF11)

Other eligible clinical targets, including those with incomplete dates, are removed from the background. First observed RCC trial dates use all phases. Features are not frozen before 2016: this is retrospective, not prospective validation.

PARP1/PARP2/PARP3 share one gene family and drug class, so the set holds about 4 independent targets, not 6; the bootstrap interval treats every gene as independent and is too narrow.

| Score | Gold genes | AUROC (95 % CI) | Recall@100 | Recall@500 | Median rank | Enrichment@100 | Gold genes in top 100 |
|---|---|---|---|---|---|---|---|
| target profile | 6 | 0.886 (0.84–0.94) | 17% | 17% | 2,198 | 31.8× | CCND1 |
| biomarker profile | 6 | 0.668 (0.55–0.82) | 0% | 17% | 6,977 | 0.0× | — |
| convergence v3 (incl. literature + drugs) | 6 | 0.861 (0.82–0.92) | 0% | 17% | 3,027 | 0.0× | — |
| convergence v3 without literature/tractability | 6 | 0.628 (0.45–0.80) | 0% | 17% | 7,657 | 0.0× | — |
| target profile without tractability | 6 | 0.643 (0.48–0.82) | 17% | 17% | 6,602 | 31.8× | CCND1 |
| predicted tractability alone | 6 | 0.847 (0.79–0.90) | 0% | 0% | 2,669 | 0.0× | — |
| target profile, tractability without ligand history | 6 | 0.666 (0.46–0.88) | 0% | 17% | 5,657 | 0.0× | — |
| predicted tractability without ligand history, alone | 6 | 0.555 (0.33–0.75) | 0% | 0% | 8,811 | 0.0× | — |

### Prognostic biomarkers: ClearCode34 (33 of 34 genes mapped; Brooks et al. 2014)

ClearCode34 was validated in TCGA, which the biomarker profile's prognosis component also uses; the last row removes that component. A literature-based gold set (verbatim `biomarker_of` relations from the disease-level query) was empty (0 genes): all such relations came from gene-specific queries chosen from the omics ranking, so they were not used.

| Score | Gold genes | AUROC (95 % CI) | Recall@100 | Recall@500 | Median rank | Enrichment@100 | Gold genes in top 100 |
|---|---|---|---|---|---|---|---|
| target profile | 33 | 0.713 (0.62–0.79) | 6% | 18% | 4,598 | 11.6× | EPAS1, TLR3 |
| biomarker profile | 33 | 0.739 (0.68–0.81) | 0% | 3% | 4,256 | 0.0× | — |
| convergence v3 (incl. literature + drugs) | 33 | 0.701 (0.62–0.79) | 3% | 6% | 3,999 | 5.8× | EPAS1 |
| convergence v3 without literature/tractability | 33 | 0.644 (0.55–0.75) | 0% | 6% | 5,140 | 0.0× | — |
| biomarker profile without prognosis | 33 | 0.686 (0.61–0.77) | 0% | 6% | 5,507 | 0.0× | — |

## (c) Robustness to the weights

1,000 random re-weightings (Dirichlet, κ = 20, around the configured weights). Genes in the top 10 in ≥ 80 % of draws:

- **convergence**: CA9, COXFA4L2
- **target**: SCARB1, CA9, SLC16A3, CD70, NNMT, ENPP3, FABP6
- **biomarker**: PSAT1, PCK1, SLC6A19

Re-weighting never removes a stream. Leaving one stream out at a time (its weight set to 0), these top-25 genes drop out of the top 25 (the stream they depend on in brackets):

- **convergence**: SCARB1 (genetics), PIGR (genetics), NNMT (literature), UMOD (literature), FABP1 (omics), NPHS2 (omics), MT1G (omics/literature), HILPDA (omics/literature/replication), ABAT (literature), ENO2 (omics/literature), FCER1G (omics/literature), SLC16A3 (omics/literature/replication), VHL (cnv/literature/genetics), SLC12A1 (omics/tractability/replication), SLC22A6 (omics/tractability/replication), CASR (omics/tractability/replication), PLIN2 (omics/literature/replication), TYROBP (omics/literature/replication), SLC6A19 (omics/literature/tractability/replication), SLC47A2 (omics/cnv/tractability/replication)
- **target**: SLC16A3 (omics_up), HK2 (omics_up), SLC2A1 (omics_up), IGFBP3 (tumor_intrinsic/tractability_predicted), ALDOC (tumor_intrinsic/tractability_predicted), LOX (omics_up/tumor_intrinsic/tractability_predicted), PYGL (omics_up/tractability_predicted), VWF (omics_up/tumor_intrinsic/tractability_predicted), ADM (tumor_intrinsic/tractability_predicted), FCGR3A (omics_up/tractability_predicted), IDO1 (omics_up/tumor_intrinsic/tractability_predicted), LDHA (omics_up/tumor_intrinsic/tractability_predicted), COXFA4L2 (omics_up/tumor_intrinsic), CCND1 (omics_up/tumor_intrinsic/genetics/tractability_predicted)
- **biomarker**: SLPI (prognosis), DUSP9 (omics), VWF (prognosis), L1CAM (omics/replication), CA4 (prognosis/grade), G6PC1 (omics/prognosis), HSD11B2 (omics/replication/prognosis/grade), CYP2J2 (omics/replication/ccab), SPTBN2 (omics/replication/ccab), UCHL1 (replication/prognosis/grade/ccab), FCGR3A (omics/replication/grade/measurable), CA9 (omics/replication/ccab/measurable), FCER1G (omics/prognosis/grade/measurable), PROM2 (omics/replication/prognosis/grade/ccab), ALDOB (omics/replication/prognosis), PTH1R (omics/replication/grade/ccab), PLCL1 (omics/prognosis/grade/measurable), STRA6 (omics/replication/prognosis/ccab), COL23A1 (omics/replication/ccab/measurable/cnv)

## (d) What the agentic layer adds

| | Closed-book LLM (Sonnet 5.5, no tools) | GraphRAG pipeline |
|---|---|---|
| Runs | 3 | 3 |
| Mean pairwise Jaccard of candidate lists | 0.667 | 0.862 |
| Candidates that are RCC clinical targets | 60%, 60%, 60% | 15%, 15%, 14% |
| Candidates in the data-only target or biomarker top 500 | 80%, 70%, 80% | 100%, 100%, 100% (by construction: the agents pick candidates from these profiles) |
| Statements citing checkable evidence | 0 | 118 checked, 0 unsupported |

Closed-book candidates: VHL, EPAS1, VEGFA, KDR, CD274, PBRM1, BAP1, SETD2, CA9, MET; VHL, EPAS1, VEGFA, PDCD1, CD274, PBRM1, BAP1, SETD2, MET, MTOR; EPAS1, VHL, VEGFA, KDR, PDCD1, MTOR, PBRM1, BAP1, SETD2, CA9.

Pipeline final targets: CA9, CD70, NNMT, SLC16A3, ENPP3, SCARB1, EPAS1, PSAT1, PCK1, COXFA4L2, ALDOB, COL23A1, RAB42; CA9, CD70, NNMT, SCARB1, EPAS1, ENPP3, SLC16A3, PSAT1, PCK1, ALDH6A1, COXFA4L2, COL23A1, RAB42; CA9, CD70, NNMT, EPAS1, ENPP3, SLC16A3, SCARB1, PSAT1, PCK1, ALDOB, COXFA4L2, SLC22A8, COL23A1, RAB42.

The pipeline's Jaccard mostly reflects the deterministic ranking the agents start from, not agreement reached independently; it shows the agent step is stable, not that its choices are right.

### No-agent baseline

Rule: alternate target profile, biomarker profile, discovery tumour track; same length as the agent list. Overlap with the agents' final list (Jaccard): 0.529. RCC clinical targets: 8% of the baseline list versus 15% of the agents' list.

Baseline list: SCARB1, PSAT1, COL23A1, CA9, PCK1, RAB42, SLC16A3, ALDH6A1, TNFAIP6, CD70, SLC6A19, FABP6, NNMT.

Added by the agents: ENPP3, EPAS1, COXFA4L2, ALDOB. Left out by the agents: ALDH6A1, TNFAIP6, SLC6A19, FABP6.

This is what the agent step changes. Whether the changes are improvements needs the reasons in the report and, ultimately, experiments; the overlap alone does not say.

## (e) How often the verifier catches a wrong statement

Each of the 118 statements of `runs/2026-11-03` that pass verification was copied and corrupted one way at a time; a corruption counts as caught when the statement is rejected or loses evidence ids (the rule `verify.py final` applies).

| Corruption | Applied | Caught | Rate |
|---|---|---|---|
| Number: last digit of the first printed number changed | 72 | 72 | 100% |
| Direction word flipped (higher/lower, up/down, ...) | 54 | 15 | 28% |
| ... of which the statement cites significant omics_de evidence | 19 | 14 | 74% |
| Statement attributed to another gene | 105 | 105 | 100% |
| Evidence of another gene cited | 105 | 105 | 100% |

The direction check only reads tumour-versus-normal wording against cited omics_de records and printed log2FC signs, so a flipped direction in a literature or treatment sentence is not caught. A changed last digit is accepted when the new value still rounds to the stored one. Missed examples are listed in `evals/verifier_sabotage.json`.

## (f) Known-answer questions (`claude -p` + project MCP servers)

- Model: saved claude-sonnet-5-5 answers (regex regraded)
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
