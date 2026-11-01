"""Tests for the treatment, subtype, ADC, immune, replication, literature and report-help additions."""
import json
import re

import numpy as np
import pandas as pd
import pytest

from ccrcc.paths import PROCESSED

needs = lambda name: pytest.mark.skipif(not (PROCESSED / name).exists(), reason=f"{name} not built")  # noqa: E731


def test_selective_dependency_needs_three_lines():
    from ccrcc.layers import selective_dependency
    cfg = {"min_lines": 10, "padj": 0.1, "min_dependent_lines": 3}
    dep = pd.DataFrame({"n_ccrcc_lines": [24, 24, 24, 24, 5], "padj": [0.01, 0.01, 0.01, 0.5, 0.01],
                        "selectivity": [0.4, 0.4, 0.4, 0.4, 0.4], "n_ccrcc_dependent": [5, 2, 3, 6, 6],
                        "common_essential": [False, False, False, False, False]})
    assert selective_dependency(dep, cfg).tolist() == [True, False, True, False, False]
    dep.loc[0, "common_essential"] = True
    assert not selective_dependency(dep, cfg).iloc[0]


def test_selective_dependency_falls_back_to_fraction():
    from ccrcc.layers import selective_dependency
    cfg = {"min_lines": 10, "padj": 0.1, "min_dependent_lines": 3}
    dep = pd.DataFrame({"n_ccrcc_lines": [24, 24], "padj": [0.01, 0.01], "selectivity": [0.4, 0.4],
                        "frac_ccrcc_dependent": [3 / 24, 1 / 24], "common_essential": [False, False]})
    assert selective_dependency(dep, cfg).tolist() == [True, False]


def test_mann_whitney_auc_matches_scipy():
    from scipy import stats

    from ccrcc.treatment import _auc
    rng = np.random.default_rng(1)
    x = rng.normal(size=(5, 40))
    y = np.array([1] * 15 + [0] * 25)
    x[:, y == 1] += 0.8
    auc, p = _auc(x, y)
    for i in range(5):
        u = stats.mannwhitneyu(x[i, y == 1], x[i, y == 0], alternative="two-sided", use_continuity=False)
        assert auc[i] == pytest.approx(u.statistic / (15 * 25))
        assert p[i] == pytest.approx(u.pvalue, rel=0.02)


@needs("treatment_summary.json")
def test_published_signatures_reproduce():
    s = json.loads((PROCESSED / "treatment_summary.json").read_text())
    assert s["controls_reproduced"]["angiogenesis_longer_pfs_on_sunitinib"]
    assert s["controls_reproduced"]["t_effector_favours_atezo_bev_over_sunitinib"]


@needs("treatment_summary.json")
def test_no_diverged_cox_fits_are_reported():
    d = pd.read_parquet(PROCESSED / "layer_treatment.parquet")
    assert d["cox_log_hr"].abs().max() < 10


def test_surface_score_cases():
    from ccrcc.adc import surface_score
    tm = {"n_transmembrane": 1, "extracellular_residues": 300, "gpi_anchor": False, "locations": ["Cell membrane"]}
    assert surface_score(tm) == 1.0
    assert surface_score({**tm, "extracellular_residues": 5}) == 0.7
    assert surface_score({"n_transmembrane": 0, "extracellular_residues": 0, "gpi_anchor": True,
                          "locations": ["Cell membrane"]}) == 0.7
    assert surface_score({"n_transmembrane": 0, "extracellular_residues": 0, "gpi_anchor": False,
                          "locations": ["Cytoplasm"]}) == 0.0
    assert surface_score(None) == 0.0


@needs("adc_suitability.parquet")
def test_adc_controls_rank_in_top_decile():
    s = json.loads((PROCESSED / "adc_summary.json").read_text())
    assert s["control_median_percentile"] < 0.1
    d = pd.read_parquet(PROCESSED / "adc_suitability.parquet")
    assert d["suitability"].between(0, 1).all()


def test_alias_terms_are_disambiguated():
    from ccrcc.discovery import alias_index, gene_terms_v2
    hg = pd.DataFrame({"prev_symbol": ["OLDA|SHAR", None, "SHAR"], "alias_symbol": ["ALIASONE", None, None]},
                      index=["GENEA", "OLDA", "GENEB"])
    idx = alias_index(hg)
    t = gene_terms_v2("GENEA", hg, idx)
    assert '"GENEA"' in t and '"ALIASONE"' in t
    assert "SHAR" not in t          # shared by two genes
    assert "OLDA" not in t          # approved symbol of another gene


def test_short_word_symbols_need_context():
    from ccrcc.discovery import gene_terms_v2
    hg = pd.DataFrame({"prev_symbol": [None], "alias_symbol": [None]}, index=["MET"])
    assert "gene OR protein" in gene_terms_v2("MET", hg, {})
    hg2 = pd.DataFrame({"prev_symbol": [None], "alias_symbol": [None]}, index=["CA9"])
    assert "gene OR protein" not in gene_terms_v2("CA9", hg2, {})


@needs("immune_classes.parquet")
def test_immune_classes_are_tertiles():
    d = pd.read_parquet(PROCESSED / "immune_classes.parquet")
    c = d["immune_class"].value_counts()
    assert set(c.index) == {"cold", "intermediate", "hot"}
    assert c.max() - c.min() <= 3
    m = d.groupby("immune_class")["immune_score"].mean()
    assert m["cold"] < m["intermediate"] < m["hot"]


@needs("subtype_ccab_calls.parquet")
def test_ccab_control_worse_survival_for_ccb():
    s = json.loads((PROCESSED / "subtypes_summary.json").read_text())["ccab"]
    assert s["control_tcga_os_ccB_vs_ccA"]["hr"] > 1.5
    assert s["control_tcga_os_ccB_vs_ccA"]["pvalue"] < 0.001


@needs("subtype_ccab_de.parquet")
def test_ccab_component_excludes_clearcode_genes():
    from ccrcc.outcomes import clearcode34
    from ccrcc.scores import components
    u = pd.read_parquet(PROCESSED / "gene_universe.parquet")
    u["ot_genetic"] = 0.0
    t = components(u)
    assert (t.loc[t["gene"].isin(clearcode34()), "c_ccab"] == 0).all()
    assert t["c_ccab"].between(0, 1).all() and (t["c_ccab"] > 0).any()


@needs("sc_replication.parquet")
def test_second_dataset_confirms_ca9_is_tumour_cell():
    d = pd.read_parquet(PROCESSED / "sc_replication.parquet").set_index("gene")
    assert d.loc["CA9", "li_top"] == "RCC" and d.loc["CA9", "zhang_top"] == "RCC"
    assert d.loc["KDR", "zhang_top"] == "Endothelial"


@needs("spatial_viewer.json")
def test_viewer_payload_is_consistent():
    v = json.loads((PROCESSED / "spatial_viewer.json").read_text())
    assert len(v["sections"]) >= 2 and "CA9" in v["genes"]
    for s in v["sections"]:
        n = s["n_spots"]
        assert len(s["x"]) == len(s["y"]) == len(s["neighbors"]) == n
        assert all(len(a) == n for a in s["shares"].values())
        assert all(len(a) == n for a in s["expr"].values())
        assert s["image"].startswith("data:image/jpeg;base64,")
        assert max(max(nb, default=0) for nb in s["neighbors"]) < n
        # shares are percentages of a composition
        tot = np.sum([s["shares"][g] for g in v["groups"]], axis=0)
        ok = (np.abs(tot - 100) <= 7) | (tot == 0)     # spots with no marker signal keep an all-zero composition
        assert ok.all() and (tot > 0).mean() > 0.9


def test_every_template_plot_and_help_key_has_text():
    import pathlib

    from ccrcc.report_help import HELP
    base = pathlib.Path(__file__).parents[1] / "src" / "ccrcc" / "report_assets"
    text = "".join(p.read_text(encoding="utf-8") for p in base.glob("*.j2"))
    keys = set(re.findall(r"plot\('(\w+)'", text)) | set(re.findall(r"q\(['\"](\w+)['\"]\)", text))
    missing = sorted(k for k in keys if k not in HELP)
    assert not missing, missing
    for k, h in HELP.items():
        assert set(h) == {"what", "read", "watch"} and all(h.values()), k
        assert "—" not in "".join(h.values()), k          # no em dashes in help text


def test_report_template_is_light_only():
    import pathlib
    t = (pathlib.Path(__file__).parents[1] / "src" / "ccrcc" / "report_assets" / "report.html.j2").read_text(encoding="utf-8")
    assert "prefers-color-scheme:dark" not in t.replace(" ", "")
    assert 'data-theme="dark"' not in t
