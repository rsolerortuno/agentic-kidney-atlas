"""M13 checks: extra layers, convergence v3, profiles, sensitivity, semantic verification, outcome evaluation.

Unit tests (statistics, regexes) run without data; the biology checks need `ccrcc run-all` first.
"""
import json

import numpy as np
import pandas as pd
import pytest
from scipy import stats

from ccrcc.layers import cox_univariate, paired_de
from ccrcc.outcomes import auroc
from ccrcc.paths import PROCESSED
from ccrcc import semantic

needs_data = pytest.mark.skipif(not (PROCESSED / "layers_summary.json").exists(), reason="run `ccrcc run-all` first")


# ------------------------------------------------------------------ statistics (no data needed)
def test_paired_de_matches_scipy():
    rng = np.random.default_rng(0)
    t = pd.DataFrame(rng.normal(1, 1, (5, 30)), columns=[f"c{i}" for i in range(30)])
    n = pd.DataFrame(rng.normal(0, 1, (5, 30)), columns=t.columns)
    n.iloc[0, :5] = np.nan                                   # missing values are dropped pairwise
    de = paired_de(t, n).set_index("feature")
    for i in range(5):
        ok = n.iloc[i].notna()
        ref = stats.ttest_rel(t.iloc[i][ok], n.iloc[i][ok])
        assert de.loc[i, "pvalue"] == pytest.approx(ref.pvalue, rel=1e-6)
        assert de.loc[i, "log2fc"] == pytest.approx((t.iloc[i][ok] - n.iloc[i][ok]).mean())


def test_cox_matches_statsmodels():
    from statsmodels.duration.hazard_regression import PHReg
    rng = np.random.default_rng(1)
    n = 300
    z = rng.normal(size=(n, 3))
    T = np.round(rng.exponential(np.exp(-0.7 * z[:, 0])) * 50)   # ties on purpose
    E = (rng.random(n) < 0.7).astype(int)
    beta, p = cox_univariate((z - z.mean(0)) / z.std(0), T.astype(float), E)
    for j in range(3):
        zz = (z[:, j] - z[:, j].mean()) / z[:, j].std()
        ref = PHReg(T, zz[:, None], status=E, ties="breslow").fit()
        assert beta[j] == pytest.approx(ref.params[0], abs=1e-4)
        assert p[j] == pytest.approx(ref.pvalues[0], rel=1e-3)
    assert beta[0] > 0.3 and p[0] < 1e-6                     # the simulated effect is recovered


def test_auroc():
    s = np.array([5, 4, 3, 2, 1.0])
    assert auroc(s, np.array([1, 1, 0, 0, 0], bool)) == 1.0
    assert auroc(s, np.array([0, 0, 0, 1, 1], bool)) == 0.0
    assert auroc(np.ones(4), np.array([1, 0, 1, 0], bool)) == 0.5


# ------------------------------------------------------------------ semantic verification (no data needed)
EV = [{"type": "omics_de", "gene": "CA9", "layer": "rna", "log2fc": 6.3123, "padj": 1.72e-49, "significant": True,
       "n_pairs": 75}]


def test_numbers_match_with_rounding_percent_and_scientific():
    ok = semantic.check("CA9 is up in tumours (log2FC 6.31, padj 1.7e-49, 75 pairs).", "CA9", EV)
    assert ok["unmatched_numbers"] == [] and ok["direction_contradiction"] is None
    sc = [{"type": "sc_celltype_expression", "gene": "CA9", "share_of_tumor_transcripts": 0.976}]
    assert semantic.check("97.6 % of its tumour transcripts", "CA9", sc)["unmatched_numbers"] == []


def test_wrong_number_and_wrong_direction_are_caught():
    bad = semantic.check("CA9 is up with log2FC 4.2.", "CA9", EV)
    assert bad["unmatched_numbers"] == ["4.2"]
    contra = semantic.check("CA9 is strongly down-regulated in tumours.", "CA9", EV)
    assert contra["direction_contradiction"]


def test_significance_threshold_is_not_a_claim():
    assert semantic.numbers_in("q 0.061, not significant at 0.05") == [
        {"text": "0.061", "value": 0.061, "decimals": 3, "sci": False, "percent": False}]


def test_identifiers_and_exempt_numbers_are_not_checked():
    txt = "CA9 (HIF-2α target, 3p, PMID:36419844, Phase 3, XL-820, IODINE I 124 GIRENTUXIMAB) in 2019."
    assert semantic.numbers_in(txt) == []


# ------------------------------------------------------------------ biology recovered by the new layers
@needs_data
def test_mutation_layer_recovers_drivers_and_corrects_length():
    m = pd.read_parquet(PROCESSED / "layer_mutations.parquet").set_index("gene")
    for g in ("VHL", "PBRM1", "SETD2", "BAP1", "KDM5C"):
        assert m.loc[g, "padj"] < 1e-10, g
    assert m.loc["VHL", "freq"] > 0.4
    assert m.loc["TTN", "padj"] > 0.5        # long gene, mutated often, not recurrent after length correction


@needs_data
def test_dependency_layer_lineage_and_hif():
    d = pd.read_parquet(PROCESSED / "layer_dependency.parquet").set_index("gene")
    assert d.loc["PAX8", "selectivity"] > 0.5 and d.loc["PAX8", "padj"] < 0.01      # renal lineage factor
    assert d.loc["HNF1B", "selectivity"] > 0.5
    assert d.loc["HIF1A", "mean_effect_ccrcc"] > 0                                  # tumour suppressor in ccRCC
    lines = json.loads((PROCESSED / "layer_dependency_lines.json").read_text())["ccrcc_lines"]
    assert {"786-O", "769-P"} <= {x["name"] for x in lines}


@needs_data
def test_tcga_replicates_hif_targets():
    t = pd.read_parquet(PROCESSED / "layer_tcga.parquet").set_index("gene")
    for g in ("CA9", "NDUFA4L2", "ANGPTL4", "EGLN3", "VEGFA"):
        assert t.loc[g, "log2fc"] > 2 and t.loc[g, "padj"] < 1e-10, g


@needs_data
def test_purity_separates_infiltration_from_tumour_expression():
    p = pd.read_parquet(PROCESSED / "layer_purity.parquet").set_index("dataset_symbol")
    de = pd.read_parquet(PROCESSED / "de_results.parquet")
    rna = de[de.layer == "rna"].set_index("gene")["log2fc"]
    assert p.loc["PTPRC", "rna_purity_rho"] < -0.5 and p.loc["PTPRC", "rna_lfc_at_purity1"] < rna["PTPRC"] / 2
    assert p.loc["CD74", "rna_lfc_at_purity1"] < rna["CD74"] / 2
    assert p.loc["CA9", "rna_lfc_at_purity1"] >= rna["CA9"]


@needs_data
def test_kinase_activity_protein_adjusted():
    k = pd.read_parquet(PROCESSED / "layer_kinase_activity.parquet")
    assert not k["gene"].isin(["EGF", "AVP", "GH1"]).any()     # ligands are not kinases
    assert k["n_substrate_sites"].min() >= 5


# ------------------------------------------------------------------ convergence v3, profiles, sensitivity
@needs_data
def test_v3_components_bounded_and_continuous():
    from ccrcc.integrate import compute
    t = compute()
    comps = [c for c in t.columns if c.startswith("c_")]
    assert (t[comps].min().min() >= 0) and (t[comps].max().max() <= 1.0000001)
    # continuous omics: far fewer ties at the maximum than the binary v2 component (549 genes at 1.0)
    assert (t["c_omics"] >= 0.999).sum() < 100
    # 3p passengers no longer share one CNV bonus
    assert t["c_cnv"].round(3).nunique() > 50
    rank = {g: i + 1 for i, g in enumerate(t["gene"])}
    assert rank["CA9"] <= 5 and rank["VHL"] <= 50


@needs_data
def test_profiles_answer_different_questions():
    from ccrcc.integrate import compute
    t = compute()
    tgt = t.sort_values("score_target", ascending=False)["gene"].head(100).tolist()
    bio = t.sort_values("score_biomarker", ascending=False)["gene"].head(100).tolist()
    assert "CA9" in tgt[:10]
    down = t.set_index("gene")["rna_log2fc"]
    assert all(down.get(g, 0) > 0 for g in tgt[:50])            # targets are tumour-up genes
    assert any(down.get(g, 0) < 0 for g in bio[:50])            # biomarkers include lost genes


@needs_data
def test_sensitivity_intervals():
    s = pd.read_parquet(PROCESSED / "sensitivity_convergence.parquet")
    assert (s["rank_p05"] <= s["rank_median"]).all() and (s["rank_median"] <= s["rank_p95"]).all()
    assert s["p_top10"].between(0, 1).all() and s.iloc[0]["p_top10"] > 0.5


@needs_data
def test_new_layers_in_graph_have_evidence():
    from ccrcc.graph import GraphStore
    st = GraphStore().stats()
    assert st["claims_without_resolvable_evidence"] == 0
    for ev in ("somatic_recurrence", "depmap_dependency", "tcga_de", "tcga_survival", "kinase_activity", "purity"):
        assert st["evidence_by_type"].get(ev, 0) > 0, ev
