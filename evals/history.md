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


## M12: target discovery (2026-10-31)

Annotation only; recall@100 is unchanged. Decisions and findings:

- Pool: concordant tumor-up RNA (log2FC >= 1) and protein (>= 0.5), dominant scRNA cell group holding >= 40 % of tumour transcripts and
  expressing in >= 10 % of its cells. Tracks by dominant group: tumor cell / stroma / immune.
- Novel = no clinical-stage drug in Open Targets and <= 15 Europe PMC records for gene AND ccRCC terms. First-in-class adds SM/AB tractability and
  <= 50 nTPM in every vital organ (HPA). Weights are in `config.yaml`; none was tuned against an outcome (there is none: this is hypothesis generation).
- Found and fixed: (1) Open Targets modality codes (SM, AB, PR, OC) were compared to names, so no gene was druggable; PR is called for ~90 % of genes, so
  it now only adds 0.2. (2) FABP6 passed as first-in-class with 4,264 nTPM in small intestine; vital-organ expression now blocks that status.
  (3) The omics component could go negative; clipped and covered by a unit test.
- MGD score positive controls (known glue substrates, scored like the candidates): GSPT1 3rd and CSNK1A1 10th of 54 proteins, but SALL4 36th,
  ZFP91 47th, IKZF3 48th and IKZF1 last: beta-hairpin degrons are captured, CRBN zinc-finger degrons are not (low pLDDT, no cavity). The score is
  therefore only informative for proteins with a ligandable cavity.
- Process error: rebuilding the graph while `ccrcc.structure` was reading it locked SQLite and crashed that run; stages that read the graph must not
  overlap with `kg_build`.

## M13: v3 convergence, independent layers, outcome evaluation (2026-10-31)

Requested after a review of v2. All weights and saturation constants below were written into `config.yaml` **before**
measuring anything; nothing was tuned afterwards against recall or the outcome gold sets. Two changes were made after
looking at results, both to remove a leak or a bug, and are listed under "Corrections".

| Version | What changed | Recall@100 | Without genetics stream | Median rank |
|---|---|---|---|---|
| v2 | (above) | 22 % | 4 % | 309 |
| v3 | continuous omics; CNV × cis effect; measured somatic recurrence replaces the OT somatic datatype; tractability for all genes; + DepMap dependency (1.0) and TCGA replication (1.0) | **39 %** | 4 % | **171** |

recall@50/200/500/1000 = 26 / 52 / 70 / 78 % (v2: 17 / 35 / 52 / 74 %).

**Why v3 recovers more known genes.** The somatic-recurrence layer (CPTAC + TCGA-KIRC, 511 tumours, binomial against a
coding-length background, hypermutators excluded) calls exactly the canonical drivers: VHL (247/511), PBRM1, BAP1,
SETD2, KDM5C, MTOR, PTEN, TP53, and leaves long genes such as TTN (80 mutated tumours, obs/exp 0.93) and MUC16 out.
VHL moves from rank 67 to 16, SETD2 to 49, PBRM1 to 57, MTOR to 56, PTEN to 62, TP53 to 40. The genes still missed are
not expression- or mutation-recurrent in these cohorts (HIF1A, a tumour suppressor here; FLCN and KMT2C) or are known
only from germline / pathway knowledge (HNF1A/B, FLT4, ATM).

**Why recall@100 is still not the right target.** The 60 % goal mixes drivers, targets and biomarkers in one list, and
the top 20 Open Targets genes in it are partly what the genetics stream measures (removing that stream leaves 4 %). M13
therefore adds two question-specific profiles and evaluates them against outcomes they do not see (`results.md` (b)):

- **Target profile vs RCC clinical targets** (75 genes; phase ≥ 2 in RCC, ≥ 2 trials or approval, drugs with ≤ 5
  targets): AUROC 0.90, 13 in the top 100 (33× enrichment). **Temporal holdout** (13 targets whose first RCC trial
  began in 2016 or later: CCND1, CDK4/6, EPAS1, FOLH1, GLS, IDO1, LAG3, PARP1-3, TIGIT, TNFSF11): AUROC 0.91, 4 of 13
  in the top 100 (IDO1, CCND1, EPAS1, CDK6).
- **Honest reading of that AUROC.** Predicted tractability alone reaches 0.85: clinical targets are druggable proteins
  by construction. What the data layers add is early enrichment: with random order inside the tractable set, 4 % of
  clinical targets fall in the top 100; the target profile puts 17 % there (holdout: 8 % vs 31 %). Without the
  tractability component the target profile still places 11 % (holdout 23 %) in the top 100, at AUROC 0.72.
- **Biomarker profile vs ClearCode34** (Brooks et al. 2014; 33 of 34 genes mapped): AUROC 0.75 (0.70 without the TCGA
  prognosis component, which shares a cohort with ClearCode34's validation). Weak early enrichment: the profile ranks
  strong tumour-vs-normal changes, ClearCode34 separates ccA from ccB tumours, which is a within-tumour question.
- The overall convergence score includes literature and drug history; its clinical-target AUROC (0.89) is partly
  circular and is reported only to show that.

**Robustness.** 1,000 Dirichlet re-weightings (κ = 20) per ranking. CA9 is the only gene in the convergence top 10 in
≥ 80 % of draws; the target profile has seven such genes (SCARB1, CA9, SLC16A3, CD70, NNMT, ENPP3, FABP6), the
biomarker profile three (PSAT1, NNMT, ALDH6A1). Several high convergence ranks (ALDOB, SLC22A8, PIGR, VHL) have 90 %
rank intervals spanning 1 to > 150: their position depends on the chosen weights.

**Semantic verification.** Statements must now quote numbers that exist in the evidence they cite (or the gene's row of
`convergence.json`, or a count from the project summaries) and must not contradict the direction of the cited omics
evidence. Applied to the archived run that produced the published M12 report (`runs/archive/headless_m12_2026-10-31/`,
written during M11), 22 of 96 statements fail; the v2 verifier had passed all of them. Two causes: numbers taken from a
record that is not cited (e.g. CA9's phospho log2FC with only the RNA and protein ids; KDR's whole-model pLDDT cited
against a pocket record) and numbers that went stale when M12 rebuilt the spatial layer with more genes (different
target-high spots and permutation draws) while the narrative was not regenerated. In the new official run 3 of 94
statements fail (two citation errors, one stale narrative copy of them) and are dropped; 91 remain, 0 unsupported.

**Corrections made after looking at results** (none changes a weight):

1. The target profile first used "predicted tractability OR clinical precedence". Clinical precedence is drug history
   and leaks into the clinical-target evaluation, so the component now uses predicted SM / AB tractability only.
2. The first version of the numeric check accepted any number found in the project summary files; with M13 that pool
   contained thousands of decimals and matched almost anything (a deliberately wrong "log2FC 4.2" passed the unit
   test). It now accepts only integer counts from summaries. A second bug exempted every number after a word ending in
   "p", "q", "m" or "v", so p- and q-values were never checked; exemptions are now whole words.

**What the agents add (`results.md` (d)).** A closed-book run of the same question (Sonnet 5.5, no tools, 3 runs) lists
VHL, EPAS1, VEGFA, PBRM1, BAP1, SETD2, CD274/PDCD1, KDR, MTOR, CA9, MET: the textbook answer (60 % of its genes are RCC
clinical targets; mean pairwise Jaccard 0.67 across runs) with no checkable evidence. Only 18 % of the pipeline's genes
are clinical targets, but all of them are in the data-only target or biomarker top 500; the closed-book union and
the pipeline list share two genes, CA9 and EPAS1 (Jaccard 0.10). The closed-book model recalls what has been tried; the pipeline proposes what these data support,
and the orchestrator promoted EPAS1 (rank 38, the approved belzutifan target) back into the list on clinical grounds. The pipeline's final lists are driven by the data (CA9, CD70, NNMT, SCARB1, PSAT1, ENPP3,
SLC16A3, ALDH6A1, ALDOB, SLC22A8, EPAS1) and every statement cites verified evidence. Three independent `/informe`
runs on the same graph (official + `runs/replicates/rep1`, `rep2`) wrote different text and different agent findings
(summary text similarity 0.61 between two of them) but converged on the same 11 genes; rep2 added ENO2 and reordered
the middle of the list (mean pairwise Jaccard 0.94; CA9, CD70, NNMT first in all three). 91, 93 and 87 statements
were checked, 0 unsupported in each; the semantic check dropped 3, 0 and 2 statements respectively.
