"""M12 checks: target discovery scoring, safety and clinical layers."""
import json

import numpy as np
import pandas as pd
import pytest

from ccrcc.config import load
from ccrcc.discovery import _gene_terms, score
from ccrcc.graph import GraphStore
from ccrcc.paths import PROCESSED

needs_data = pytest.mark.skipif(not (PROCESSED / "discovery_candidates.parquet").exists(),
                                reason="run `ccrcc run-all` first")


def _frame(**over):
    base = dict(gene=["A", "B"], rna_log2fc=[3.0, 3.0], protein_log2fc=[1.5, 1.5], source_share=[0.9, 0.9],
                sc_tvp_log2fc=[4.0, 4.0], sc_tvp_padj=[0.01, 0.01], ccrcc_papers=[2, 500], total_papers=[50, 5000],
                clinical_drugs=[[], []], any_drugs=[[], []], modalities=[["AB"], ["AB"]], surface_or_secreted=[True, True],
                vital_max_ntpm=[2.0, 2.0], clin_rho=[0.4, 0.4])
    base.update(over)
    return pd.DataFrame(base)


def test_understudied_beats_well_studied():
    s = score(_frame()).set_index("gene")
    assert s.loc["A", "novelty_score"] > s.loc["B", "novelty_score"]
    assert bool(s.loc["A", "novel"]) and not bool(s.loc["B", "novel"])
    assert s.loc["A", "status"] == "first-in-class opportunity"


def test_clinical_drug_removes_novelty():
    s = score(_frame(clinical_drugs=[["X"], []])).set_index("gene")
    assert not bool(s.loc["A", "novel"]) and s.loc["A", "status"] == "has clinical drug (repurposing)"


def test_vital_organ_expression_blocks_first_in_class():
    s = score(_frame(vital_max_ntpm=[400.0, 2.0])).set_index("gene")
    assert s.loc["A", "status"].startswith("novel, off-tumor risk")
    assert s.loc["A", "c_safety"] < s.loc["B", "c_safety"]


def test_untractable_is_harder_to_drug():
    s = score(_frame(modalities=[[], ["AB"]], surface_or_secreted=[False, True])).set_index("gene")
    assert s.loc["A", "status"] == "novel, harder to drug" and s.loc["A", "c_druggable"] == 0


def test_components_are_bounded():
    s = score(_frame(rna_log2fc=[20.0, -3.0], protein_log2fc=[9.0, 0.0], clin_rho=[0.99, -0.8]))
    for k in ("omics", "cell_specificity", "understudied", "undrugged", "druggable", "safety", "aggressiveness"):
        assert s[f"c_{k}"].between(0, 1).all(), k


def test_gene_terms_use_specific_aliases_only():
    hg = pd.DataFrame({"prev_symbol": ["NDUFA4L2|ab"], "alias_symbol": ["xx|HIGD2B1"]}, index=["COXFA4L2"])
    t = _gene_terms("COXFA4L2", hg)
    assert 'TITLE_ABS:"COXFA4L2"' in t and 'TITLE_ABS:"NDUFA4L2"' in t and '"ab"' not in t and '"xx"' not in t


@needs_data
def test_pool_satisfies_its_own_thresholds_and_novel_definition():
    cfg = load()["discovery"]
    d = pd.read_parquet(PROCESSED / "discovery_candidates.parquet")
    assert (d.rna_log2fc >= cfg["min_rna_log2fc"]).all() and (d.protein_log2fc >= cfg["min_protein_log2fc"]).all()
    assert (d.source_share >= cfg["min_source_share"]).all()
    novel = d[d.novel]
    assert (novel.ccrcc_papers <= cfg["novel_max_ccrcc_papers"]).all() and (novel.clinical_drugs == "").all()
    assert set(d.track) <= {"tumor", "stroma", "immune"} and d.novelty_score.between(0, 1).all()
    assert d.gene.is_unique


@needs_data
def test_known_established_targets_are_not_called_novel():
    d = pd.read_parquet(PROCESSED / "discovery_candidates.parquet").set_index("gene")
    for g in ("CA9", "MET", "ANGPTL4", "NNMT"):
        if g in d.index:
            assert not bool(d.loc[g, "novel"]), g       # heavily studied in ccRCC (and CA9/MET have clinical drugs)


@needs_data
def test_shortlist_is_novel_tractable_and_covered_by_structure():
    short = json.loads((PROCESSED / "discovery_shortlist.json").read_text())
    d = pd.read_parquet(PROCESSED / "discovery_candidates.parquet").set_index("gene")
    assert short and d.loc[short, "novel"].all() and d.loc[short, "druggable_core"].all()
    st = {r["gene"] for r in json.loads((PROCESSED / "structure_summary.json").read_text())["results"] if r["modeled"]}
    assert len(set(short) - st) <= 2          # at most a couple without a UniProt/AlphaFold entry


@needs_data
def test_hpa_and_clinical_layers_behave():
    n = pd.read_parquet(PROCESSED / "normal_tissue.parquet").set_index("ensembl")
    assert n.loc["ENSG00000163631", "vital_max_ntpm"] > 1000          # ALB: liver
    c = pd.read_parquet(PROCESSED / "clinical_associations.parquet")
    assert {"rna", "protein"} <= set(c.layer) and {"grade", "stage"} == set(c.variable)
    assert c.rho.abs().max() <= 1 and c.padj.dropna().between(0, 1).all()


@needs_data
def test_discovery_evidence_is_in_the_graph_and_resolvable():
    store = GraphStore()
    d = pd.read_parquet(PROCESSED / "discovery_candidates.parquet")
    row = store.conn.execute("SELECT count(*) FROM evidence WHERE type = 'discovery_profile'").fetchone()[0]
    assert row >= 0.9 * len(d)
    assert store.stats()["claims_without_resolvable_evidence"] == 0
