# Evaluation history

The primary metric (spec M8) is recall of 23 known ccRCC genes (10 curated + Open Targets top 20,
deduplicated) in the top 100 of a 19,161-gene convergence ranking; target ≥ 60 %.

| Version | What changed | Recall@100 | Without genetics/OT stream |
|---|---|---|---|
| v1 | Pre-registered design: omics stream (RNA + protein + phospho + concordance) 2.0, CNV 0.5, local literature 1.0, Open Targets overall score 2.0, tractability 0.5 | 13 % | 4 % |
| v1.1 | **Bug fix**: druggability was only fetched for genes chosen with stale weights, so approved-drug targets (EPAS1/belzutifan, MTOR) scored tractability 0. Enrichment now uses the same scoring as `integrate.py` and always includes the Open Targets top 100 | 22 % | 4 % |
| v2 | **One methodological change, decided before measuring**: the Open Targets *overall* score also contains literature and expression evidence that other streams already measure (double counting). Replaced by a genetics stream = max(OT germline genetic association, OT somatic mutation); literature = max(local verbatim-validated relations / 3, OT text-mining literature) to correct the selection bias of the local corpus (PubMed searched only for the top omics genes) | 22 % | 4 % |

No further weight tuning was done against this metric, to avoid fitting the ranking to the evaluation set.

## Why recall@100 stays at 22 %

1. **Expression convergence and driver status are nearly orthogonal in ccRCC.** Most known genes (VHL, PBRM1,
   SETD2, BAP1, KDM5C, TP53, PTEN, ATM, KMT2C…) are drivers through *mutation* or 3p loss. Their bulk RNA or
   protein levels barely change (e.g. VHL protein log2FC −0.17, PBRM1 −0.27), so the omics stream (the
   largest) gives them ~0. Thousands of genes have strong, concordant tumor-vs-normal changes (549 genes have
   a perfect omics component), and they fill the top 100.
2. **The spec excludes somatic mutation data** (only the 8 LinkedOmics files: clinical, CNV, RNA, proteome,
   phosphoproteome). The only genetic signal left is Open Targets plus 3p loss frequency, which is shared by
   ~490 genes on 3p and therefore not specific.
3. **A uniform Open Targets somatic score (0.456)** is attached to many cancer-census genes (EGFR, CD74, FAS,
   CXCR4, MET, PTPRC…), inflating the genetics stream for genes with no ccRCC-specific mutation evidence.

The ranking still enriches the known genes strongly: 42× over random at top 100, median rank 309/19,161
(top 1.6 %), recall@500 52 %, recall@1000 74 %.

## What would fix it (future work, outside this spec)

- Add the CPTAC/TCGA somatic mutation layer as its own stream (mutation frequency, driver calls).
- Score genes against two separate questions, "driver/mechanism" and "expression-level target or biomarker",
  instead of one convergence score. The evaluation list mixes both.
- Down-weight non-specific Open Targets evidence (uniform pan-cancer somatic scores).
