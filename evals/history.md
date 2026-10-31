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

## M11: single-cell, spatial and structure layers (2026-10-31)

These layers **annotate** candidates; they do not change the convergence score, so recall@100 (22 %) is unchanged by
construction. Decisions taken before looking at results:

- Datasets: Li et al. 2022 (*Cancer Cell*) as curated by CZ CELLxGENE Discover: scRNA-seq (270,855 cells, 12 donors;
  178,817 used after excluding blood, adrenal, fat, oncocytoma and the authors' "MHCII" RCC sub-cluster, which
  CELLxGENE annotated as monocytes) and 23 Visium sections (38,836 spots) from tumour core and tumour-normal interface.
- Cell-type attribution uses raw UMI pseudobulks (donor × tissue × group), not single cells, so p-values are over donors.
- Spatial deconvolution excludes the candidate genes from the marker sets; correlations of CA9/KDR/CD74/PTPRC with the
  matching deconvolved cell share are therefore a held-out check (all positive, tested in `tests/test_cell_structure.py`).
- Neighbour enrichment is observed/expected against section-matched random spots, so abundant cell types are not favoured.

**Design flaw found and fixed in the MGD score.** First run: the lysine component saturated at 4 accessible lysines and
80 % of pockets reached 1.0, so it carried no information (the ranking was effectively confidence × geometry). Fixed by
saturating at 8 (`structure.lysine_saturation`), a change made from the component's distribution, not from any gene's
rank. Remaining caveat: the score is a heuristic prioritization of sites with no ground truth for glue degraders; the only
validation available is whether detected cavities recover experimentally ligand-bound residues (PDBe-KB).
