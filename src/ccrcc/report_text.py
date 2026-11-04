"""Plain-language text blocks of the report: tab intros and the guided tour on the Overview page.

Written for someone who knows a little biology but not this pipeline. Nothing here is computed.
"""

PLAIN: dict[str, str] = {
    "candidates": "This tab ranks every gene by how many separate lines of evidence point to it. A gene scores high when different kinds of data agree: it changes in tumours, a copy-number change really moves its product, papers support it, genetics point to it, the protein looks druggable, cancer cells depend on it, and a second patient cohort sees the same change. The table lists the top 25 and the plots show where each score comes from.",
    "omics": "This tab is the starting point of the whole analysis. It compares tumour tissue with healthy kidney tissue taken from the same patients, for RNA, protein and phosphorylation (the on and off switches on proteins). A gene that is much higher in tumours is a possible drug target or marker. A gene that is much lower is usually a kidney gene that the tumour has lost.",
    "cells": "A tumour sample is a mix of cancer cells, blood vessels, connective tissue and immune cells. A gene that looks high in the mix may be made by any of them. This tab uses single-cell data to say which cell makes each gene, and tissue maps to say where in the tumour it sits and what is next to it. That matters because a drug aimed at a gene made by blood-vessel cells is a different drug from one aimed at cancer cells.",
    "structure": "Explore candidate protein shapes predicted by AlphaFold and their detected cavities, alongside experimental molecular-glue complexes. Grey shows the protein and red marks the selected residues. Pocket scores summarize structural features; they do not predict binding or degradation. Experimental degron annotations are shown separately.",
    "discovery": "Most top-ranked genes are already well known. This tab looks for the opposite: genes that rise in tumours, are made by one cell type, have little cancer literature and no clinical drug. Each gene is checked for safety (is it also high in vital organs?), druggability and independent confirmation. Genes that pass every filter are called first-in-class opportunities. They are hypotheses to test, not proven targets.",
    "viewer": "This is the one place where you can look at real tissue. Pick a gene and see where it is switched on across a tumour slice, then click a spot to see which cells are around it.",
    "response": "Two questions that matter to patients and doctors. Does a gene help tell who will do better on which drug? And does a gene behave differently in tumours with different genetic changes or different aggressiveness? The answers come from two clinical trials and from the two large tumour cohorts.",
    "adcimmune": "Antibody-drug conjugates and radioligands are drugs that hunt for a protein on the cancer cell surface and deliver a toxin or radiation. They only work when the protein is on tumour cells, on most of them, and missing from vital organs. The first half of this tab scores genes on that. The second half asks which genes belong to the immune cells in the tumour, since those can be treatment targets too.",
    "validation": "A ranking is only worth trusting if it survives hard questions. This tab checks the candidates against other cohorts and other kinds of data, shows how much the order changes when the weights change, and tests whether the scores can recover real clinical targets they were never shown.",
    "literature": "This tab covers what is already known. The literature part shows what PubMed abstracts say about the top genes. The drug part shows which drugs already exist, how far they got in clinical trials, and how druggable each protein looks.",
    "network": "The network draws the top genes, the drugs that act on them and the genes known to interact. It helps you see clusters and neighbours. It is not proof of mechanism.",
    "agents": "The atlas is written by a team of AI agents rather than by one. This tab tests whether that is worth it. It explains, step by step, two experiments that compared one agent with several, how the answers were marked, what each option cost and what went wrong, and how the team was then improved so that its reports are as complete as one agent's while staying cheaper and faster. It is about how the atlas is made, not about any gene.",
    "methods": "This tab holds the receipts: every dataset with its checksum, the models used, the limits of the analysis, and a searchable list of the evidence records the report cites.",
}

TOUR: list[tuple[str, str, str]] = [
    ("candidates", "Candidates", "Which genes have the most agreeing evidence? The ranking, how each score is built and how stable the order is."),
    ("omics", "Omics", "What changes in tumours compared with healthy kidney? RNA, protein, copy number and grade."),
    ("cells", "Cells & space", "Which cell makes each gene, and what sits next to it in the tissue? Repeated in a second dataset."),
    ("structure", "Structure & MGD", "Do the proteins have pockets a small molecule could use? A heuristic, shown with its validation."),
    ("discovery", "Discovery", "Which genes are new, undrugged and still safe? The first-in-class screen with all its filters."),
    ("viewer", "Tissue viewer", "Real tissue with its H&E picture. Pick a gene, click a spot, see its neighbours."),
    ("response", "Treatment & subtypes", "Do genes track response to immunotherapy or sunitinib, or differ with VHL, PBRM1, BAP1, SETD2 and ccA/ccB?"),
    ("adcimmune", "ADC & immune", "Which proteins suit antibody-drug conjugates or radioligands, and which genes belong to the immune microenvironment?"),
    ("validation", "Validation & robustness", "Do the calls hold in other cohorts and under different weights? Includes the confirmatory cohort."),
    ("literature", "Literature & drugs", "What do papers say, and which drugs already exist?"),
    ("network", "Network", "Genes, drugs and interactions drawn as a graph."),
    ("agents", "Agent benchmark", "Is a team of AI agents better than one? Two experiments, the costs, the errors and the fixes."),
    ("methods", "Methods & evidence", "Datasets, limits and every cited evidence record."),
]
