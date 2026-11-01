"""Help entries for panels added later (confirmatory cohort, new tables, agent findings). Merged into HELP by report_help."""

EXTRA_HELP: dict[str, dict[str, str]] = {
    "confirmatory_protein": {
        "what": "A check of every protein change against a second, independent set of patients: the confirmatory CPTAC cohort, measured with a different mass-spectrometry method (label-free DIA).",
        "read": "Each dot is a protein. The horizontal axis is the change in the discovery cohort and the vertical axis the change in the confirmatory cohort. Dots along the diagonal agree. Dots in the top right and bottom left replicate in direction.",
        "watch": "The two methods compress fold changes differently, so judge direction and significance, not whether the numbers are equal. Proteins that are missing in normal kidney are dropped when fewer than 10 matched pairs remain."},
    "confirmatory_phospho": {
        "what": "The same check for phosphosites, the places on proteins where a phosphate tag sits.",
        "read": "Same reading as the protein plot. Each dot is one site on one protein.",
        "watch": "Fewer pairs are available than for proteins, and multi-site peptides were skipped."},
    "confirmatory_table": {
        "what": "The candidates one by one: change in discovery, change in the confirmatory cohort, and whether the call replicates.",
        "read": "Replicates means the same direction and an adjusted p-value below 0.05 in the confirmatory cohort.",
        "watch": "Genes that are absent or too sparse in the confirmatory proteome are not listed. That is a gap in the data and does not count as a failure."},
    "applied_findings": {
        "what": "Short findings written by an AI agent from the tools that serve this tab, after automatic checking.",
        "read": "Each sentence cites evidence records. Open a bracketed number to see the exact values.",
        "watch": "A sentence is kept only if its numbers and its direction match the cited record."},
    "cell_table": {
        "what": "For each candidate, the cell type with the highest expression and how much of the gene's tumour RNA that cell type makes.",
        "read": "A high tumour-cell share means a tumour-cell target. The last column compares the tumour cells with the normal cell they came from.",
        "watch": "Shares depend on how many cells of each type the experiment caught."},
    "spatial_table": {
        "what": "In the tissue slices, which cell type explains each gene's signal, and which neighbours are significantly enriched next to its high spots.",
        "read": "Positive log2 values mean more of that cell type next to the gene than chance.",
        "watch": "Spots mix cells and neighbouring spots overlap."},
    "pocket_table": {
        "what": "The best predicted pocket of each protein with its size, depth, confidence and score.",
        "read": "The score combines coordinate confidence, cavity geometry, accessible lysines and exposed loops. Inspect each component and the local pLDDT.",
        "watch": "Ligand-contact recovery assesses cavity detection only. The score does not predict binding or degradation."},
    "controls_table": {
        "what": "Monomer structural summaries for known molecular-glue substrates.",
        "read": "Compare the descriptive structural features; rank does not measure degradation susceptibility.",
        "watch": "No matched experimental negative set is included. This is not predictive validation."},
    "profile_table": {
        "what": "The top 25 genes of the target or biomarker profile with each component.",
        "read": "Components run from 0 to 1. The rank interval and the probability of staying in the top ten show how stable the call is.",
        "watch": "Neither profile uses literature or drug history."},
    "baseline_table": {
        "what": "A comparison of a language model that answers from memory with the full pipeline.",
        "read": "The pipeline gives candidates that are supported by data, with every statement checked. The memory-only model gives the textbook list.",
        "watch": "Only three runs of each, so differences between small numbers mean little."},
    "agents_table": {
        "what": "Where the orchestrator agreed with or departed from the numeric ranking.",
        "read": "Promoted means the gene is in the final list but outside the convergence top ten.",
        "watch": "The reasons are in the candidate cards on the Overview tab."},
    "whatsnew": {
        "what": "A short list of what this version adds, grouped by purpose.",
        "read": "Each item names the tab where you can see it.",
        "watch": "One criterion was changed after seeing data. It is listed so you can judge it."},
}
