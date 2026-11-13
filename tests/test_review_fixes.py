"""Regressions found in the 2026-11-13 review: immune classes, single-cell claim direction, evidence descriptions."""
import json
import sqlite3

import pandas as pd
import pytest

from ccrcc.paths import GRAPH_DB, PROCESSED


@pytest.mark.skipif(not (PROCESSED / "immune_classes.parquet").exists(), reason="run ccrcc.immune")
def test_immune_classes_use_tumours_only():
    from ccrcc.omics import LAYERS, read_matrix
    cls = pd.read_parquet(PROCESSED / "immune_classes.parquet")
    tumours = set(read_matrix(LAYERS["rna"][0]).columns)
    assert not cls["case"].str.endswith("-A").any()          # adjacent normal samples of the WashU files
    assert set(cls["case"]) <= tumours and len(cls) == len(tumours)
    n = cls["immune_class"].value_counts()
    assert set(n.index) == {"cold", "intermediate", "hot"} and n.max() - n.min() <= 1    # thirds of the tumours
    s = json.loads((PROCESSED / "immune_summary.json").read_text())
    assert s["tumours"] == len(cls)


@pytest.mark.skipif(not GRAPH_DB.exists(), reason="build the graph")
def test_single_cell_direction_is_claimed_only_when_significant():
    c = sqlite3.connect(GRAPH_DB)
    rows = c.execute("SELECT c.predicate, e.props FROM claims c JOIN evidence e ON e.id = c.evidence_id "
                     "WHERE e.type = 'sc_tumor_vs_pt'").fetchall()
    if not rows:
        pytest.skip("single-cell layer not loaded")
    for pred, props in rows:
        p = json.loads(props)
        if p["padj"] < 0.05:
            assert pred == ("upregulated_in" if p["log2fc"] > 0 else "downregulated_in")
        else:
            assert pred == "not_changed_in"


def test_every_m16_record_type_has_a_description():
    from ccrcc.render import evidence_view
    np_rec = {"evidence_id": "EV:normal_proteome:x", "type": "normal_proteome", "source": "s", "gene": "CA9",
              "vital_tissues_detected": 3, "vital_protein_max_rel": 1.2, "vital_protein_max_tissue": "Stomach",
              "vital_protein_ts_max": 4.5, "vital_protein_ts_tissue": "Stomach", "vital_protein_enriched": True}
    pr_rec = {"evidence_id": "EV:prism_response:x", "type": "prism_response", "source": "s", "gene": "KDR",
              "compound": "axitinib", "moa": "VEGFR inhibitor", "compound_targets": 3, "n_compounds": 20,
              "selectivity": 0.4, "mean_lfc_ccrcc": -1.1, "mean_lfc_other": -0.7, "n_ccrcc_lines": 12, "padj": 0.2}
    for rec, words in ((np_rec, ("CA9", "Stomach", "4.5")), (pr_rec, ("KDR", "axitinib", "0.4"))):
        detail = evidence_view(rec)["detail"]
        assert detail and all(w in detail for w in words)


@pytest.mark.skipif(not (PROCESSED / "confirmatory_protein.parquet").exists(), reason="run ccrcc.confirmatory")
def test_confirmatory_tables_use_current_symbols():
    from ccrcc.omics import load_hgnc
    current = set(load_hgnc()["symbol"])
    p = pd.read_parquet(PROCESSED / "confirmatory_protein.parquet")
    assert p["gene"].is_unique
    assert {"AARS1", "ADSS1", "MMUT"} <= set(p["gene"]) and not {"AARS", "ADSSL1", "MUT"} & set(p["gene"])
    assert (p["dataset_symbol"] != p["gene"]).sum() >= 150         # renamed proteins are kept, not dropped
    s = pd.read_parquet(PROCESSED / "confirmatory_phospho.parquet")
    assert s["site_id"].is_unique and s["gene"].isin(current).all()


@pytest.mark.skipif(not (PROCESSED / "spatial_summary.json").exists(), reason="run ccrcc.spatial")
def test_section_count_in_the_texts_matches_the_data():
    from ccrcc.report_help import HELP
    n = json.loads((PROCESSED / "spatial_summary.json").read_text())["sections"]
    texts = " ".join(v for h in HELP.values() for v in h.values() if isinstance(v, str))
    import re
    assert set(re.findall(r"(\d+) (?:tissue )?sections", texts)) <= {str(n), "4"}

