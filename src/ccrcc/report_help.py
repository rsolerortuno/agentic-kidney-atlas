"""Plain-language help for the "?" button on every figure, table and counter of the report.

Each entry has three parts. "what": what the panel is and why it is there. "read": the axes, colours and symbols, and
what to look for. "watch": what we see in this analysis and what could mislead. Written for someone who knows a little
biology but not this pipeline. Nothing here is computed; numbers always come from the tables and figures.
"""

HELP: dict[str, dict[str, str]] = {
    # ------------------------------------------------------------------ overview counters
    "kpi_universe": {
        "what": "The number of human genes that received a score.",
        "read": "A gene is counted if it has an approved gene name and at least one measurement in the data.",
        "watch": "A bigger number only means a wider search. Genes missing from the proteomics still get scores from the other layers."},
    "kpi_claims": {
        "what": "Statements stored in the knowledge graph, such as 'gene X is higher in tumours'.",
        "read": "Each statement is linked to an evidence record with the value, the dataset and the date. The bracketed numbers in the text open those records.",
        "watch": "Most come from tests run on every gene, so the count is large without every claim being interesting."},
    "kpi_cells": {
        "what": "Single cells from Li et al. 2022 that tell us which cell type makes each gene.",
        "read": "Each cell is one dot in the data. More cells make the cell-type shares steadier.",
        "watch": "The number of patients (12) matters more than the number of cells. Tissue is broken up before sequencing, which favours immune cells."},
    "kpi_spots": {
        "what": "Visium spots: 55 micrometre circles of tissue whose RNA is read together, from 16 tissue sections.",
        "read": "A spot holds a handful of cells, so it is a local average and not a single cell.",
        "watch": "Neighbouring spots share cells, so they are not independent measurements."},
    "kpi_proteins": {
        "what": "Proteins with a predicted 3D structure (AlphaFold) that went through the pocket analysis.",
        "read": "It is only a count of structures analysed.",
        "watch": "Models are predictions, not measured structures."},
    "kpi_fic": {
        "what": "Genes that pass every filter of the discovery screen for a first-in-class opportunity.",
        "read": "To pass, a gene must have no clinical drug, little ccRCC literature, a druggable protein, low expression in vital organs, the same change in TCGA and, outside the immune track, a change that survives correction for tumour purity.",
        "watch": "Most of them are immune-cell genes, because the filters leave many of those standing. Each one is a hypothesis to test, not a finding."},
    # ------------------------------------------------------------------ candidates
    "top25": {
        "what": "The 25 genes with the highest convergence score, with the score split into its parts.",
        "read": "Each row is a gene. The columns Omics to Repl. are the seven evidence streams on a 0 to 1 scale, and the score adds them up after weighting. The text below the table explains how every column is calculated. The rank interval is the range the gene's rank takes when the weights are changed at random.",
        "watch": "A high score with only one or two filled columns means one kind of evidence carries the gene. The literature and drug columns favour genes that are already well studied."},
    "stream_stack": {
        "what": "How much each evidence stream adds to the score of each of the top 25 genes.",
        "read": "Vertical axis: genes, best at the top. Horizontal axis: contribution to the score (weight times component). Each colour is one stream, so a bar made of many colours means several kinds of evidence agree.",
        "watch": "CA9 is long and multicoloured. Weights were set before the results were seen, but they are still a choice."},
    "component_heatmap": {
        "what": "The same seven streams as a grid, one row per gene.",
        "read": "Horizontal axis: the seven streams. Vertical axis: genes. Colour: strength of the component from 0 (pale) to 1 (dark).",
        "watch": "An empty cell means no evidence from that stream. It does not mean evidence against. Some genes simply have no data in a layer."},
    # ------------------------------------------------------------------ omics
    "volcano_rna": {
        "what": "Every gene's RNA change in tumour compared with matched healthy kidney, against how sure we are of that change.",
        "read": "Horizontal axis: log2 fold change (0 = no change, +1 = twice as high in tumour, -1 = half). Vertical axis: how small the adjusted p-value is (higher = more certain). Blue diamonds are the top convergence genes. Genes far to the right and high up are strongly and reliably higher in tumours.",
        "watch": "With 75 matched pairs, even tiny differences are significant, so look at the horizontal distance too. The two arms, left and right, are lost kidney genes and gained tumour genes."},
    "volcano_protein": {
        "what": "The same comparison for protein abundance measured by mass spectrometry (84 pairs).",
        "read": "Same axes as the RNA plot. Fewer dots appear because proteomics sees fewer genes.",
        "watch": "Protein is closer to what a drug hits than RNA, but it covers fewer genes and compresses large changes."},
    "rna_vs_protein": {
        "what": "Whether a gene that changes in RNA changes in protein too.",
        "read": "Horizontal axis: RNA log2 fold change. Vertical axis: protein log2 fold change. Dots near the diagonal agree. Dots in the top right are up in both and dots in the bottom left are down in both.",
        "watch": "Disagreement can be real biology (the protein is stable or unstable) or measurement noise, so genes that agree are the safer candidates."},
    "arm_loss": {
        "what": "How often each chromosome arm is lost in the 110 tumours.",
        "read": "Horizontal axis: chromosome arms. Vertical axis: the median share of genes on that arm that are lost across tumours (0 to 1). The highlighted bar is arm 3p.",
        "watch": "Arm 3p, which carries VHL and PBRM1, should be the tallest, and it is. Arm-level calls are coarse and say little about a single gene."},
    "grade": {
        "what": "Whether a gene goes up or down as tumours become more aggressive, measured by histologic grade G1 (mildest) to G4 (most aggressive).",
        "read": "Each pair of bars is one gene: blue is RNA, orange is protein. The horizontal axis is the Spearman correlation between the gene and grade. To the right of the vertical line, the gene is higher in more aggressive tumours. To the left, it is lower. A bar near 0 means no link.",
        "watch": "The caption under the chart names the strongest cases. Correlations are modest (the largest are around 0.4), and the cohort has no survival data, so grade is the only clinical read-out here. The survival analysis uses TCGA."},
    # ------------------------------------------------------------------ cells and space
    "dotplot": {
        "what": "For each candidate, which cell types express it in the tumour.",
        "read": "Vertical axis: genes. Horizontal axis: cell types. Dot size: the share of cells of that type in which the gene is detected. Dot colour: the average expression level.",
        "watch": "A big dark dot under 'ccRCC tumor cell' marks a tumour-cell gene. Under endothelial or immune columns, the signal comes from the microenvironment. Dissociation under-captures tumour cells."},
    "share_stack": {
        "what": "Where the tumour's RNA of each gene comes from, split by cell type.",
        "read": "Each bar is one gene and adds up to 100 %. Each colour is a cell type. Horizontal axis: share of the gene's transcripts in tumour tissue.",
        "watch": "A mostly-blue bar means tumour cells make the gene, which suits a tumour target. A mostly other-colour bar means a blood-vessel or immune gene. Shares depend on how many cells of each type were captured."},
    "neighbors": {
        "what": "Which cell types sit next to the places where a gene is strongly expressed.",
        "read": "Rows: genes. Columns: cell types. Colour: log2 of observed over expected share of that cell type in the ring of neighbouring spots. Red means more than chance, blue means less. Cells with a label passed the statistical test (adjusted p < 0.05).",
        "watch": "The ring is about 100 micrometres, so this describes a neighbourhood and not cell-to-cell contact. Chance is estimated from random spots of the same sections."},
    "map_core": {
        "what": "A real tumour-core section drawn spot by spot.",
        "read": "Each hexagon is a spot, placed where it was in the tissue. Use the menu above the map to colour by a cell type's estimated share or by a gene's expression. Darker means more.",
        "watch": "In the tumour core the tumour-cell share is high almost everywhere. The Tissue viewer tab shows the same kind of data on top of the microscope picture."},
    "map_interface": {
        "what": "A section at the border between tumour and normal kidney.",
        "read": "Same reading as the core map. Compare the two: the tumour-cell share is high in one corner and low across the rest.",
        "watch": "One section is one patient, so do not generalise from a single map."},
    "sc_replication": {
        "what": "A check of the cell-type call in a second, independent single-cell dataset (Zhang et al. 2021).",
        "read": "Each dot is a discovery-screen gene. Horizontal axis: the share of its tumour RNA that comes from tumour cells in Li 2022. Vertical axis: the same share in Zhang 2021. Dots on the diagonal agree. Colour is the gene's track.",
        "watch": "Tumour and stroma genes agree well. Immune genes agree less, because the two studies capture very different amounts of immune cells. Both studies group cells differently, so seven shared groups were used."},
    # ------------------------------------------------------------------ structure
    "mgd_bubble": {
        "what": "The best predicted pocket of each protein, described by size and structural-context score.",
        "read": "Horizontal axis: pocket volume on a log scale. Vertical axis: the structural-context score from 0 to 1. Bubble size: number of exposed lysines near the pocket. Colour: gene group.",
        "watch": "The score only describes the pocket's surroundings. It does not show that a drug exists or that the protein can be degraded."},
    "mgd_components": {
        "what": "How the structural-context score is built from its parts.",
        "read": "Each bar is a protein. Segments are the weighted components: model confidence, pocket geometry, nearby lysines and an exposed loop. Horizontal axis: how much each adds to the score.",
        "watch": "Weights were fixed before ranking. A high score from one component alone is weaker than a balanced one."},
    "ligand_validation": {
        "what": "A test of the pocket detector against proteins whose bound ligands are known from experiments.",
        "read": "Each bar is a protein. Horizontal axis: the fraction of experimentally ligand-bound residues that fall inside a detected pocket (all pockets, and the top three). Longer bars mean better coverage.",
        "watch": "Coverage is modest. This tests the pocket detector only, not whether a protein can be degraded."},
    "control_ranking": {
        "what": "Known reference proteins scored the same way as the candidates.",
        "read": "Each bar is a reference protein, with the score on the vertical axis. If references score no better than other proteins, the score does not single them out.",
        "watch": "Some references score low, which limits how far the score can be trusted."},
    # ------------------------------------------------------------------ discovery
    "discovery_scatter": {
        "what": "The discovery screen on one page: how well studied each gene is against how strongly it rises in tumours.",
        "read": "Horizontal axis: ccRCC papers plus 1, on a log scale (left = little studied). Vertical axis: protein log2 fold change in tumour (up = higher in tumour). Marker size: how specific the gene is to one cell type. Colour: track. The best candidates are towards the top left.",
        "watch": "Literature counts depend on gene names and can be off for ambiguous symbols. Only aliases that belong to one gene are counted."},
    "discovery_bars": {
        "what": "The pieces of the novelty score.",
        "read": "Each bar is one gene. Segments are the components (change in tumour, cell specificity, little literature, no drug, druggability, safety and others). Horizontal axis: contribution to the score.",
        "watch": "Weights were set before ranking. Scores are compared inside a track, not across tracks."},
    "safety": {
        "what": "Whether the novel candidates are also expressed in vital organs.",
        "read": "Horizontal axis: the highest RNA level in any vital organ (Human Protein Atlas, nTPM, log scale). Vertical axis: novelty score. Dots left of the dashed line are quiet in vital organs.",
        "watch": "This is RNA and not protein. The ADC tab adds protein staining. Most novel tumour-cell genes are expressed in at least one vital organ."},
    "discovery_table": {
        "what": "The full screen, filterable by track and status.",
        "read": "Use the chips to pick a track or status and the box to search. Each row is a gene. A star next to a TCGA value means it did not replicate, and a down arrow means the change disappears when tumour purity is set to one.",
        "watch": "Rank is computed inside each track, so a rank of 1 in the immune track is not comparable to a rank of 1 in the tumour track."},
    "chemical_matter": {
        "what": "Whether potent compounds are already known for the protein, even when no drug is in the clinic.",
        "read": "'Preclinical compounds' means at least three molecules with activity of 1 micromolar or better in ChEMBL, or patent records. Those genes get half the undrugged credit.",
        "watch": "ChEMBL covers published and patented chemistry, not private company programmes."},
    # ------------------------------------------------------------------ validation
    "tcga": {
        "what": "Whether the change seen in CPTAC shows up in a different set of patients (TCGA-KIRC).",
        "read": "Horizontal axis: RNA log2 fold change in CPTAC. Vertical axis: the same in TCGA-KIRC. Dots on the diagonal replicate. The title gives the correlation across all genes.",
        "watch": "The correlation is high (about 0.8), so the tumour-versus-normal signal is not a quirk of one cohort. TCGA has RNA only, so protein is checked in the confirmatory cohort."},
    "purity": {
        "what": "Whether a change belongs to tumour cells or to the normal and immune cells mixed into the sample.",
        "read": "Each gene has two points. The first is the change measured in the whole sample. The second is the change extrapolated to a pure tumour (100 % cancer cells). If the point stays far from zero, tumour cells make the change. If it moves towards zero, infiltrating cells were carrying it.",
        "watch": "The extrapolation assumes a straight-line effect of purity, so treat it as a model and not a measurement."},
    "mutations": {
        "what": "Genes mutated more often than their length alone would predict.",
        "read": "Horizontal axis: observed over expected number of mutated tumours (log scale). Vertical axis: how small the adjusted p-value is. Genes to the right and high up are recurrently mutated drivers. VHL, PBRM1, BAP1 and SETD2 should stand out.",
        "watch": "Long genes mutate by chance, which is why the expected count uses gene length. Hypermutated tumours were left out of the background rate."},
    "dependency": {
        "what": "Whether cancer cell lines need a gene to survive, comparing ccRCC lines with all others.",
        "read": "Horizontal axis: the average dependency score in other lines. Vertical axis: the same in ccRCC lines (lower means the cells need the gene more). A dot below the diagonal means ccRCC lines depend on the gene more than other lines do. Orange dots pass the selective-dependency rule (at least three dependent ccRCC lines).",
        "watch": "Only 24 ccRCC lines, grown in normal oxygen. No top expression target is a selective dependency, so they are handles for antibodies and not weak spots."},
    "kinases": {
        "what": "Enzymes (kinases) whose known targets are more phosphorylated in tumours.",
        "read": "Each bar is a kinase. Horizontal axis: the activity score (z), higher means more active in tumours. Coloured bars pass the statistical cut.",
        "watch": "Inferred from lists of known substrates, so well-studied kinases are favoured."},
    "rank_target": {
        "what": "How stable each gene's position in the target profile is when the weights are changed at random.",
        "read": "Each row is a gene. The dot is the median rank and the bar covers the 5th to 95th percentile of 1,000 re-weightings. Horizontal axis: rank on a log scale (left = better). A short bar near the left is a robust call.",
        "watch": "A long bar means the position depends on the chosen weights. The weights are varied around the chosen ones, not over every possible choice."},
    "rank_biomarker": {
        "what": "The same stability check for the biomarker profile.",
        "read": "Same axes as the target plot.",
        "watch": "The biomarker profile includes the ccA/ccB axis, which reorders some calls."},
    "rank_conv": {
        "what": "The same stability check for the overall convergence score.",
        "read": "Same axes as the target plot.",
        "watch": "Only a few genes stay near the top whatever the weights, with CA9 the clearest."},
    "outcomes": {
        "what": "How well each score ranks genes the clinic later used as targets, using information the score never saw.",
        "read": "Each group of bars is a set of true clinical targets. Vertical axis: AUROC (0.5 = random, 1.0 = perfect), with the 95 % interval as a whisker. Bars are different scores.",
        "watch": "Clinical targets are mostly druggable proteins, so a score that only flags druggability also looks good. The tractability-only bar shows that. The post-2016 set is small and retrospective."},
    # ------------------------------------------------------------------ literature
    "lit_relations": {
        "what": "The kinds of relations extracted from the abstracts.",
        "read": "Vertical axis: relation type (for example 'up-regulated in'). Horizontal axis: number of relations. Every relation quotes a sentence that was checked word for word against the abstract.",
        "watch": "The corpus is small and favours genes that omics already ranked high."},
    "lit_years": {
        "what": "When the abstracts were published.",
        "read": "Horizontal axis: year. Vertical axis: number of publications in the local corpus.",
        "watch": "It shows how recent the evidence is, not how important."},
    "drug_stage": {
        "what": "The most advanced clinical stage of drugs aimed at each top gene (Open Targets).",
        "read": "Each row is a gene. Longer bars mean later stages. Approval is the end of the line.",
        "watch": "A drug for the gene in another disease also counts."},
    # ------------------------------------------------------------------ network
    "network": {
        "what": "The top genes, the drugs that act on them and ccRCC, drawn as a graph.",
        "read": "Circles are genes, diamonds are drugs and the rectangle is the disease. Lines are drug to target, protein interactions and up or down changes in ccRCC. Drag nodes to move them and click a line to see its evidence id.",
        "watch": "It shows neighbours and clusters, not cause and effect."},
    # ------------------------------------------------------------------ treatment, subtypes, ADC, immune
    "treatment_arms": {
        "what": "For every gene, how much its expression relates to time to progression in two arms of the IMmotion150 trial.",
        "read": "Each dot is a gene. Horizontal axis: log hazard ratio on sunitinib. Vertical axis: log hazard ratio on atezolizumab plus bevacizumab. Below zero means longer benefit with higher expression. A dot below the diagonal does relatively better on the immunotherapy combination, and one above does relatively better on sunitinib.",
        "watch": "About 90 patients per arm. The cloud is centred and no gene separates the two arms after correcting for 20,000 tests, so any single gene is a lead."},
    "treatment_controls": {
        "what": "A sanity check of the trial data with three published gene signatures.",
        "read": "Horizontal axis: signature. Vertical axis: hazard ratio of progression per standard deviation of the signature, on a log scale (below 1 means longer benefit). Colours are the arms, and the printed p-values come from Cox regression.",
        "watch": "The angiogenesis signature helps on sunitinib and the T-effector signature helps on the combination, both as published, so the gene-level results rest on working data."},
    "treatment_table": {
        "what": "The top single genes in each arm and the genes whose effect differs between arms.",
        "read": "Log HR is per standard deviation of expression: negative means longer benefit with higher expression. q is the adjusted p-value over about 20,000 genes. The last column says the same in words.",
        "watch": "Treat every single-gene hit as a lead. None of the genes in the third table has q below 0.1. The nivolumab trial has 16 patients, so it can only check direction."},
    "subtype_counts": {
        "what": "How many genes differ between tumours with and without a mutation in each driver gene.",
        "read": "Horizontal axis: the driver gene. Vertical axis: number of genes that differ (q < 0.1). Colours are the data type: CPTAC RNA, CPTAC protein and TCGA RNA. A taller bar means the driver reshapes the tumour more.",
        "watch": "PBRM1 and BAP1 change far more than VHL, probably because many tumours labelled VHL wild type also lack VHL function. BAP1 and SETD2 mutants are few in CPTAC (about 15), so TCGA confirms."},
    "subgroup_table": {
        "what": "Genes that are higher in mutant tumours in all three datasets and also higher in tumour than in normal kidney.",
        "read": "These are candidates for one genetic subgroup. The three numeric columns are mutant minus wild type in log2 units. 'Tractable' and 'Clinical drug' come from Open Targets.",
        "watch": "A gene can be higher in BAP1-mutant tumours just because those tumours contain more immune cells."},
    "ccab": {
        "what": "Each tumour's position on the ccA to ccB scale, from the 34 ClearCode34 genes.",
        "read": "Vertical axis: ccA score minus ccB score. The violins show how the tumours of each cohort are spread. Above 0 is ccA (better outlook), below 0 is ccB (worse outlook).",
        "watch": "This is a score-based approximation of the published classifier. ccB tumours have worse survival in TCGA, which is the sanity check."},
    "adc": {
        "what": "Surface proteins as candidates for antibody-drug conjugates or radioligands.",
        "read": "Each dot is a surface protein. Horizontal axis: change in tumour (protein log2 fold change, or RNA divided by 2). Vertical axis: share of tumour cells expressing it in the weakest quarter of donors. Colour: staining in vital organs from none (pale) to high (red). Good candidates are far right, high up and pale.",
        "watch": "Dots without antibody data in the Human Protein Atlas are drawn neutral. Known ccRCC targets are labelled as controls."},
    "adc_table": {
        "what": "The same candidates ranked by suitability.",
        "read": "Suitability adds four parts: surface topology, increase in tumour, absence from vital organs and evenness across tumour cells. The vital-organ column gives the strongest staining and the tissue. Kidney staining is shown but not scored.",
        "watch": "If a control ranks low, check which column pulled it down."},
    "immune_classes": {
        "what": "Immune cell scores in tumours split into cold, intermediate and hot.",
        "read": "Horizontal axis: immune class (tertiles of the xCell immune score). Vertical axis: CIBERSORT absolute score of each cell type. Boxes show the middle half of tumours.",
        "watch": "Scores should rise from cold to hot, and they do. They are estimates from bulk RNA, not cell counts."},
    "immune_table": {
        "what": "Genes that look like immune targets.",
        "read": "The score adds the immune share of the gene's tumour RNA, its link to CD8 T cells, Tregs, M2 macrophages and the immune score, response on the combination, increase in tumour and tractability. Known checkpoints serve as controls.",
        "watch": "Drug history is not used, so checkpoints are not guaranteed to rank high, and only one of seven reaches the top 100."},
    "tissue_viewer": {
        "what": "Four real tissue sections with their microscope (H&E) picture and the Visium spots on top.",
        "read": "Choose a gene and a section. Spots in the top tenth of the gene's expression are ringed. Click a spot to outline its six neighbours. Orange bars: mean cell makeup of the neighbours of all ringed spots. Blue bars: the same over the whole section.",
        "watch": "The p-values under the viewer come from the permutation test over all 16 sections. The live bars only describe the section on screen."},
    # ------------------------------------------------------------------ confirmatory
    "confirmatory_protein": {
        "what": "A check of every protein change against a second, independent set of patients measured with a different method.",
        "read": "Each dot is a protein. Horizontal axis: change in the discovery cohort. Vertical axis: change in the confirmatory cohort. Dots along the diagonal agree and dots in the top right or bottom left replicate in direction. Blue dots replicate, grey dots do not.",
        "watch": "The two methods compress fold changes differently, so judge direction and significance and not the size. Proteins missing in normal kidney are dropped when fewer than 10 matched pairs remain."},
    "confirmatory_phospho": {
        "what": "The same check for phosphosites, the places on proteins where a phosphate tag sits.",
        "read": "Same axes as the protein plot. Each dot is one site on one protein.",
        "watch": "Fewer pairs are available than for proteins, and multi-site peptides were skipped."},
    "confirmatory_table": {
        "what": "The candidates one by one, in the discovery and in the confirmatory cohort.",
        "read": "'Replicates' means the same direction and an adjusted p-value below 0.05 in the confirmatory cohort.",
        "watch": "Genes absent or too sparse in the confirmatory proteome are not listed. That is a gap and not a failure."},
    "applied_findings": {
        "what": "Short findings written by an AI agent from the tools that serve this tab, after automatic checking.",
        "read": "Each sentence cites evidence records. Open a bracketed number to see the exact values.",
        "watch": "A sentence is kept only if its numbers and its direction match the cited record."},
    # ------------------------------------------------------------------ other tables
    "cell_table": {
        "what": "For each candidate, the cell type with the highest expression and how much of the gene's tumour RNA that cell type makes.",
        "read": "Columns: top cell type, share of its cells expressing the gene, share of the gene's tumour RNA, tumour-cell share, a specificity index (tau, 0 = everywhere, 1 = one cell type) and the tumour versus proximal tubule comparison.",
        "watch": "Shares depend on how many cells of each type the experiment captured."},
    "spatial_table": {
        "what": "In the tissue slices, which cell type explains each gene's signal and which neighbours are enriched.",
        "read": "Positive log2 values mean more of that cell type next to the gene's high spots than chance, among results with adjusted p below 0.05.",
        "watch": "Spots mix cells and neighbouring spots overlap."},
    "pocket_table": {
        "what": "The best predicted pocket of each protein with its size, depth, confidence and score.",
        "read": "Check the confidence (pLDDT) first. Then volume, depth, nearby lysines and the score.",
        "watch": "Predicted structures, a heuristic score and a modest validation. Not evidence that a drug exists."},
    "controls_table": {
        "what": "Known reference proteins scored the same way as the candidates.",
        "read": "Their score and rank among all proteins.",
        "watch": "References do not all rank well, which limits the score."},
    "profile_table": {
        "what": "The top 25 genes of the target or biomarker profile with each component.",
        "read": "Components run from 0 to 1 and are explained in the box above the table. The rank interval and the chance of staying in the top ten show how stable the call is.",
        "watch": "Neither profile uses literature or drug history, so a well-known gene can rank low."},
    "baseline_table": {
        "what": "A comparison of a language model that answers from memory with the full pipeline.",
        "read": "Rows are measures such as overlap between runs and the share of candidates that are known clinical targets.",
        "watch": "Only a few runs of each, so small differences mean little."},
    "agents_table": {
        "what": "Where the orchestrating agent agreed with or departed from the numeric ranking.",
        "read": "'Promoted' means the gene is in the final list but outside the convergence top ten.",
        "watch": "The reasons are in the candidate cards on the Overview tab."},
}
