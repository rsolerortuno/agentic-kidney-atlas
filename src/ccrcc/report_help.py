"""Plain-language help for the "?" button on every figure, table and counter of the report.

Each entry has three short parts: what the panel shows, how to read it, and what to keep in mind. Written for someone who
knows some biology but not this pipeline. Nothing here is computed; numbers always come from the tables and figures.
"""

HELP: dict[str, dict[str, str]] = {
    # ------------------------------------------------------------------ overview counters
    "kpi_universe": {"what": "How many genes got a score. It is every gene that has a valid approved symbol and at least one measurement.",
                     "read": "A bigger number does not mean a better analysis. It only says how wide the net was.",
                     "watch": "Genes missing from the proteomics or phospho data still get scores from the layers that do cover them."},
    "kpi_claims": {"what": "Statements stored in the knowledge graph, such as 'gene X is up in ccRCC'. Each one points to a record with the number, the dataset and the date behind it.",
                   "read": "The report only cites claims that have such a record. If you click a bracketed number in a paragraph, you land on that record.",
                   "watch": "Most claims come from the genome-wide tests, so the count is large without meaning that each claim is interesting."},
    "kpi_cells": {"what": "Single cells from Li et al. 2022 used to ask which cell type makes each transcript.",
                  "read": "More cells give steadier cell-type shares, but the number of donors (12) matters more than the number of cells.",
                  "watch": "Tissue is dissociated before sequencing, which favours immune cells over tumour cells."},
    "kpi_spots": {"what": "Visium spots, each a 55 micrometre circle of tissue that holds a handful of cells, from 23 tissue sections.",
                  "read": "A spot is a local average, not a single cell.",
                  "watch": "Spots near each other share cells, so neighbouring spots are not independent."},
    "kpi_proteins": {"what": "Proteins with an AlphaFold structural model that went through the pocket analysis.",
                     "read": "Only a count of how many structures were analysed.", "watch": "Models are predictions, not measured structures."},
    "kpi_fic": {"what": "Genes that pass every filter for a first-in-class opportunity: no clinical drug, little ccRCC literature, a druggable protein, low expression in vital organs, replication in TCGA, and (outside the immune track) an up-regulation that survives correction for tumour purity.",
                "read": "Treat each one as a hypothesis to test, not a finding.", "watch": "Most of them are immune-cell genes, because the filters leave many of those standing."},
    # ------------------------------------------------------------------ candidates
    "top25": {"what": "The 25 genes with the highest convergence score. Each column is one stream of evidence on a 0 to 1 scale.",
              "read": "A high score with only one or two filled columns means a single type of evidence is carrying the gene. The rank interval shows how far the rank moves when the weights are changed at random.",
              "watch": "Literature and drug columns reward what is already well studied, so they push known genes up."},
    "stream_stack": {"what": "For the top genes, how much each evidence stream adds to the final score.",
                     "read": "Each bar is one gene. The coloured segments are the streams. Wide colour variety means several kinds of data agree.",
                     "watch": "Weights are set in config.yaml before looking at results, but they are still a choice."},
    "component_heatmap": {"what": "The same streams as a grid, one row per gene.",
                          "read": "Dark cells are strong evidence. Empty cells are streams with nothing for that gene.", "watch": "Missing is not the same as negative. Some genes simply have no data in a layer."},
    # ------------------------------------------------------------------ omics
    "volcano_rna": {"what": "Every gene's change in tumour compared with matched normal kidney (RNA), against how sure we are of that change.",
                    "read": "Right is higher in tumour, left is lower. Higher up means a smaller adjusted p-value. Blue dots are the top convergence genes.",
                    "watch": "With 75 matched pairs, tiny differences become highly significant. Look at the horizontal distance too."},
    "volcano_protein": {"what": "The same comparison for protein abundance measured by mass spectrometry.",
                        "read": "Same reading as the RNA plot. Fewer genes appear because proteomics sees fewer genes.", "watch": "Protein is closer to function than RNA, but it covers fewer genes."},
    "rna_vs_protein": {"what": "Whether a gene that goes up in RNA also goes up in protein.",
                       "read": "Points along the diagonal agree. Points off the diagonal change in one layer only.", "watch": "Disagreement can be biology (protein stability) or measurement noise."},
    "arm_loss": {"what": "How often each chromosome arm is lost in the tumours.",
                 "read": "Taller bars are more frequent losses. Arm 3p, home of VHL and PBRM1, should be the highest.", "watch": "Arm-level calls come from copy-number segments, so they say little about single genes."},
    "grade": {"what": "How strongly each gene's expression follows tumour grade.",
              "read": "Positive means higher in more aggressive tumours, negative the opposite.", "watch": "There is no survival data in this cohort, so grade is the only clinical read-out."},
    # ------------------------------------------------------------------ cells and space
    "dotplot": {"what": "For each candidate, which cell type expresses it in the tumour.",
                "read": "Dot size is the share of cells of that type with the gene detected. Dot colour is the average level.", "watch": "Dissociation under-captures tumour cells, so the plot can understate tumour expression."},
    "share_stack": {"what": "Where the tumour's transcripts of each gene come from, split by cell type.",
                    "read": "A bar that is mostly one colour means one cell type makes the gene. A gene made by tumour cells is a tumour target. One made by endothelium or immune cells is not.",
                    "watch": "Shares depend on how many cells of each type were captured."},
    "neighbors": {"what": "Which cell types sit next to spots where a gene is highly expressed.",
                  "read": "Red means the cell type is found next to the gene's high spots more often than chance. Blue means less often. Chance is estimated from random spots of the same sections.",
                  "watch": "The effect is about 100 micrometres, a ring of six spots, so it is a neighbourhood and not a cell-to-cell contact."},
    "map_core": {"what": "A tumour-core section drawn from the spot data, coloured by a cell share or by a gene.",
                 "read": "Pick a variable in the menu above the map. Brighter means more.", "watch": "Use the Tissue viewer tab if you want the actual H&E picture under the spots."},
    "map_interface": {"what": "The same map for a section at the tumour-normal border.",
                      "read": "Compare with the core map to see how cell mixes change at the border.", "watch": "One section is one patient."},
    "sc_replication": {"what": "A check of the cell-type call in a second, independent single-cell dataset (Zhang et al. 2021).",
                       "read": "Each dot is a gene from the discovery screen. On the diagonal, both datasets give the same tumour-cell share.", "watch": "Zhang has fewer cells and fewer donors, and the two studies group cell types differently, so seven shared groups were used."},
    # ------------------------------------------------------------------ structure
    "mgd_bubble": {"what": "Best pocket score per protein, from the AlphaFold models.", "read": "Higher means larger weighted structural components, not a greater probability of degradation.", "watch": "The score is a heuristic for ordering sites. It does not show that any compound exists."},
    "mgd_components": {"what": "How the pocket score is built from its parts.", "read": "Each bar is one component and its weight.", "watch": "Weights were fixed before ranking."},
    "ligand_validation": {"what": "A check against proteins whose bound ligands are known from experiments.", "read": "Higher bars mean the detected pockets overlap the real ligand site more often.", "watch": "The agreement is modest."},
    "control_ranking": {"what": "Where known cases rank under the same scoring.", "read": "Ranks describe monomer features, not compound-dependent degradation.", "watch": "Known substrates alone, without matched negatives, do not validate a predictor."},
    # ------------------------------------------------------------------ discovery
    "discovery_scatter": {"what": "The discovery screen: how much ccRCC literature each gene has against how much it changes in tumour.",
                          "read": "Candidates in the lower right are strongly up in tumour and little studied.", "watch": "Literature counts depend on gene names and can be off for ambiguous symbols. Version 3 uses only aliases that belong to one gene."},
    "discovery_bars": {"what": "The weights of each component of the novelty score.", "read": "Longer bars count more.", "watch": "Weights were set before looking at the ranking."},
    "safety": {"what": "Expression of each candidate in vital organs (Human Protein Atlas RNA).", "read": "Genes to the left are quiet in vital organs, so a drug is less likely to hit healthy tissue.", "watch": "This is RNA, not protein. The ADC tab adds protein staining."},
    "discovery_table": {"what": "The full screen, filterable by track and status.", "read": "Use the chips and the search box. A star next to a TCGA value means it did not replicate. A down arrow means the change disappears when tumour purity is set to one.",
                        "watch": "Rank is computed within each track."},
    "chemical_matter": {"what": "Whether potent compounds are already known for the protein, even if no drug is in the clinic.",
                        "read": "'Preclinical compounds' means at least three molecules with activity of 1 micromolar or better in ChEMBL, or patent-derived activity records. Those genes get half the undrugged credit.",
                        "watch": "ChEMBL covers published and patented medicinal chemistry, not private company programs."},
    # ------------------------------------------------------------------ validation
    "tcga": {"what": "Whether the CPTAC change in each gene shows up in a separate cohort (TCGA-KIRC).", "read": "Points on the diagonal replicate. The correlation between cohorts is the summary.", "watch": "TCGA is RNA only, so protein is not replicated yet."},
    "purity": {"what": "How the fold change would look in a pure tumour, once the effect of normal cells mixed into the sample is removed.", "read": "A gene that keeps its change at purity one is made by tumour cells. One that loses it was carried by infiltrating cells.", "watch": "The correction assumes a straight-line effect of purity."},
    "mutations": {"what": "Genes mutated more often than their length would predict, across CPTAC and TCGA-KIRC tumours.", "read": "Far to the right means more mutations than expected. VHL, PBRM1, BAP1 and SETD2 should stand out.", "watch": "Hypermutated tumours were excluded from the background rate."},
    "dependency": {"what": "Whether cancer cell lines need the gene to survive (DepMap CRISPR screens), ccRCC lines against all others.", "read": "Below the diagonal means ccRCC lines depend on the gene more than other lines do. A selective dependency now needs at least three dependent ccRCC lines.", "watch": "Only 24 ccRCC lines, grown in normal oxygen."},
    "kinases": {"what": "Kinases whose known targets are more phosphorylated in tumours.", "read": "Longer bars mean higher inferred activity.", "watch": "Inferred from substrate lists, so well studied kinases are favoured."},
    "rank_target": {"what": "How stable each gene's rank is in the target profile when the weights are changed at random.", "read": "A short bar near the top is a robust call. A long bar means the rank depends on the weights.", "watch": "The random weights are drawn around the chosen ones, not over every possible choice."},
    "rank_biomarker": {"what": "The same stability check for the biomarker profile.", "read": "Same reading as the target plot.", "watch": "The biomarker profile now also includes the ccA/ccB axis."},
    "rank_conv": {"what": "The same stability check for the overall convergence score.", "read": "Same reading as the target plot.", "watch": "Few genes stay in the top ten in most draws."},
    "outcomes": {"what": "How well each score ranks genes that are real clinical targets, using data the score never saw.", "read": "AUROC of 0.5 is chance and 1.0 is perfect. The retrospective group uses first observed RCC trials across all phases, from 2016 onward. Features are not frozen before 2016; this is not a prospective prediction test.", "watch": "Clinical targets are mostly druggable proteins, so a score that simply flags druggability looks good too. The tractability-only row shows that."},
    # ------------------------------------------------------------------ literature
    "lit_relations": {"what": "The kinds of relations extracted from the abstracts, such as 'up-regulated in' or 'biomarker of'.", "read": "Taller bars mean more extracted relations of that type.", "watch": "Every relation quotes a sentence that was checked word for word against the abstract."},
    "lit_years": {"what": "When the abstracts were published.", "read": "Shows whether the corpus is recent or dated.", "watch": "The corpus favours genes that omics already ranked high."},
    "drug_stage": {"what": "The most advanced clinical stage of drugs aimed at each top gene (Open Targets).", "read": "Longer bars are later stages. Approval is the end.", "watch": "A drug against a gene in another disease also counts."},
    # ------------------------------------------------------------------ network
    "network": {"what": "Genes, drugs and ccRCC as a graph, with the evidence behind each edge.", "read": "Drag nodes to move them. Click an edge to see the evidence id. Circles are genes, diamonds are drugs, the rectangle is the disease.", "watch": "STRING links show known interactions, not a causal order."},
    # ------------------------------------------------------------------ M15
    "treatment_arms": {"what": "For every gene, how much its expression predicts time to progression in two trial arms of IMmotion150.",
                       "read": "Each dot is a gene. Negative numbers mean higher expression goes with longer benefit. A gene below the diagonal does better on the immunotherapy combination than on sunitinib. A gene above it does better on sunitinib.",
                       "watch": "About 90 patients per arm. The test of a real difference between arms is weaker than the test inside one arm, and no gene passes it after correction."},
    "treatment_controls": {"what": "A sanity check on the trial data using three published gene signatures.",
                           "read": "A hazard ratio below 1 means a higher signature goes with longer progression-free survival. The angiogenesis signature should help on sunitinib and the T-effector signature should help on the combination.",
                           "watch": "Both known patterns reproduce here, so the gene-level results are at least built on working data."},
    "treatment_table": {"what": "Top single genes by arm and, in the last table, genes whose effect differs between arms.", "read": "Log HR is per standard deviation of expression. Adjusted p-values are corrected over about 20,000 genes.",
                        "watch": "Treat all single-gene hits as leads. The Choueiri 2016 trial has 16 patients, so it can only check direction."},
    "subtype_counts": {"what": "How many genes differ between tumours with and without a mutation in each driver gene.", "read": "Taller bars mean the driver reshapes the tumour more. BAP1 and PBRM1 change far more than VHL, probably because many tumours labelled VHL wild type also lack VHL function (loss or methylation), which blurs that comparison.",
                       "watch": "Few tumours carry a BAP1 or SETD2 mutation (15 and 14 in CPTAC), so TCGA is used to confirm."},
    "subgroup_table": {"what": "Genes that are higher in mutant tumours in CPTAC RNA, CPTAC protein and TCGA, and also up in tumour against normal.",
                       "read": "These are possible targets for one genetic subgroup. Differences are in log2 units, mutant minus wild type.", "watch": "A gene can be higher in BAP1-mutant tumours just because BAP1-mutant tumours contain more immune cells."},
    "ccab": {"what": "Each tumour's position on the ccA to ccB scale, from the 34 genes of the ClearCode34 classifier.", "read": "Tumours above zero are called ccA (better outcome), below zero ccB. Hazard ratio of ccB against ccA in TCGA is the sanity check.", "watch": "This is a score-based approximation of the published classifier, not the original centroid method."},
    "adc": {"what": "Surface proteins as candidates for antibody-drug conjugates or radioligands.", "read": "Right means more increased in tumour. Higher means more tumour cells carry it, even in the weakest quarter of donors. Redder dots show staining in vital organs.", "watch": "Dots with no antibody data in the Human Protein Atlas are drawn as neutral."},
    "adc_table": {"what": "The same candidates ranked by suitability: surface location, tumour increase, absence from vital organs by staining, and how evenly tumour cells express them.",
                  "read": "Known ccRCC targets (CA9, CD70, ENPP3, HAVCR1, FOLH1) are marked as controls.", "watch": "If a control ranks low, look at the surface, staining and homogeneity columns to see which term pulled it down."},
    "immune_classes": {"what": "Immune cell scores in tumours split into cold, intermediate and hot by xCell.", "read": "The scores should rise from cold to hot.", "watch": "CIBERSORT and xCell are estimates from bulk RNA."},
    "immune_table": {"what": "Genes that look like immune targets: made by immune cells, tied to immune infiltration, linked to response on atezolizumab + bevacizumab, up in tumour, and tractable.",
                     "read": "Known checkpoints are shown as controls. Their ranks tell how much to trust the profile.", "watch": "Drug history is not used, so checkpoints are not guaranteed to rank high."},
    "tissue_viewer": {"what": "Four real tissue sections with their H&E picture, and the Visium spots on top.",
                      "read": "Choose a gene and a section. Spots in the top tenth of that gene are circled. Click a spot to see its six neighbours and what they are made of. The bars compare the neighbours of the high spots with the whole section.",
                      "watch": "The p-values shown below the viewer come from the permutation test over all 23 sections. The live bars are descriptive."},
}


def get(key: str) -> dict[str, str] | None:
    return HELP.get(key)


from .report_help_more import EXTRA_HELP  # noqa: E402

HELP.update(EXTRA_HELP)
