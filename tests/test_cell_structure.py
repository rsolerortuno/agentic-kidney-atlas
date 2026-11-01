"""M11 checks: single-cell, spatial and structure layers behave correctly and recover known biology."""
import json

import numpy as np
import pytest

from ccrcc.graph import GraphStore
from ccrcc.paths import PROCESSED
from ccrcc.spatial import neighborhood_enrichment
from ccrcc.structure import _fit, secondary_structure
from ccrcc.singlecell import cell_group, tau

pytestmark = pytest.mark.skipif(not (PROCESSED / "structure_summary.json").exists(), reason="run `ccrcc run-all` first")


def sc():
    return {t["gene"]: t for t in json.loads((PROCESSED / "sc_summary.json").read_text())["targets"]}


def test_cell_group_mapping():
    assert cell_group("RCC", "EMT", "unknown") == "RCC"
    assert cell_group("RCC", "MHCII", "MHC-II-positive classical monocyte") is None   # possible doublets excluded
    assert cell_group("T-cell", "CD4+Treg", "x") == "Treg" and cell_group("T-cell", "CD8+T_EM", "x") == "CD8_T"
    assert cell_group("Myeloid", "MHC-II TAM", "macrophage") == "Macrophage"


def test_tau_bounds():
    assert tau(np.array([1.0, 0, 0, 0])) == 1.0 and tau(np.array([1.0, 1, 1, 1])) == 0.0


def test_known_biology_in_single_cell():
    s = sc()
    assert s["CA9"]["top_groups"][0]["cell_group"] == "RCC" and s["CA9"]["tumor_cell_share"] > 0.8
    assert s["CA9"]["tumor_vs_pt"]["log2fc"] > 2 and s["CA9"]["tumor_vs_pt"]["padj"] < 0.05
    assert s["KDR"]["top_groups"][0]["cell_group"] == "Endothelial" and s["KDR"]["tumor_cell_share"] < 0.05
    assert s["CD74"]["tumor_cell_share"] < 0.1 and s["PTPRC"]["tumor_cell_share"] < 0.05   # immune genes


def test_spatial_deconvolution_validated_by_held_out_genes():
    v = json.loads((PROCESSED / "spatial_summary.json").read_text())["validation"]
    for k in ("spearman_CA9_vs_RCC", "spearman_KDR_vs_Endothelial", "spearman_CD74_vs_Macrophage"):
        assert v[k] > 0.1, k
    assert v["mean_share_core"]["RCC"] > v["mean_share_interface"]["RCC"]


def test_neighborhood_enrichment_detects_planted_signal():
    rng = np.random.default_rng(0)
    n = 600
    sections = np.repeat(["s1", "s2"], n // 2)
    NP = rng.dirichlet(np.ones(3) * 5, size=n)
    high = np.zeros(n, bool)
    high[rng.choice(n, 80, replace=False)] = True
    NP[high, 0] += 0.3
    NP /= NP.sum(axis=1, keepdims=True)
    e = neighborhood_enrichment(high, sections, NP, np.ones(n, bool), rng, 300)
    assert e["log2_enrichment"][0] > 0.3 and e["p"][0] < 0.01
    assert abs(e["log2_enrichment"][1:]).max() < e["log2_enrichment"][0]


def test_secondary_structure_helix_and_fit():
    t = np.arange(30)
    helix = np.c_[2.3 * np.cos(np.radians(100) * t), 2.3 * np.sin(np.radians(100) * t), 1.5 * t]
    assert (secondary_structure(helix) == "H").mean() > 0.8
    assert _fit(300, 150, 1500) == 1.0 and _fit(5000, 150, 1500) < 0.2 and _fit(0, 150, 1500) == 0.0


def test_structure_positive_controls():
    s = {r["gene"]: r for r in json.loads((PROCESSED / "structure_summary.json").read_text())["results"]}
    assert sum(r["modeled"] for r in s.values()) == len(s) >= 25
    # enzymes with experimental ligand complexes: the cavity detector (all cavities) must cover the ligand residues
    for g in ("ENO2", "ACY1", "FHIT"):
        assert s[g]["ligand_site_validation"]["recall"] >= 0.5, g
    assert s["ENO2"]["mean_plddt"] > 90 and s["IKZF1"]["mean_plddt"] < 60   # IKZF1 is mostly disordered
    for r in s.values():
        for c in r.get("pockets", []):
            assert 0 <= c["mgd_score"] <= 1 and set(c["components"]) == {"confidence", "geometry", "lysines", "degron"}


def test_graph_layers_have_resolvable_evidence():
    st = GraphStore().stats()
    assert st["claims_without_resolvable_evidence"] == 0
    for k in ("sc_celltype_expression", "sc_tumor_vs_pt", "spatial_expresser", "spatial_neighborhood",
              "alphafold_model", "structural_context"):
        assert st["evidence_by_type"][k] > 0
    assert st["nodes_by_type"]["CellType"] >= 17 and st["nodes_by_type"]["Pocket"] > 0


def test_cell_mcp_tools_return_citable_evidence():
    from ccrcc.mcp_servers import cell
    ids = [r["evidence_id"] for r in cell.celltype_expression("CA9") if r.get("evidence_id")]
    assert ids and not GraphStore().check_evidence_ids(ids)["invalid"]
    assert cell.tumor_cell_vs_proximal_tubule("CA9")["evidence_id"]
    assert cell.spatial_context("CA9")["enriched_neighbors"][0]["cell_group"] == "RCC"
    sp = cell.structure_pockets("KDR")
    assert sp["pockets"][0]["evidence_id"] and sp["model_evidence_id"]
    assert cell.structural_context_ranking(5)[0]["best_structural_context_score"] >= cell.structural_context_ranking(5)[-1]["best_structural_context_score"]


def test_cell_mcp_discovery_tools():
    from ccrcc.mcp_servers import cell
    rows = cell.novel_targets(n=5)
    assert rows and all(r["status"] == "first-in-class opportunity" and r["evidence_id"] for r in rows)
    assert not GraphStore().check_evidence_ids([r["evidence_id"] for r in rows])["invalid"]
    assert cell.novel_targets(track="tumor", status="any", n=3)[0]["track"] == "tumor"
    prof = cell.discovery_profile(rows[0]["gene"])
    assert prof["gene"] == rows[0]["gene"] and prof["evidence_id"]
    assert "not in the discovery pool" in cell.discovery_profile("VHL").get("message", "") or cell.discovery_profile("VHL")["gene"] == "VHL"
