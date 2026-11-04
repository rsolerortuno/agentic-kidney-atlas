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


@needs("confirmatory_protein.parquet")
def test_confirmatory_cohort_is_independent_and_replicates():
    from ccrcc.omics import LAYERS, read_matrix
    s = json.loads((PROCESSED / "confirmatory_summary.json").read_text())
    assert s["protein"]["replication_rate"] > 0.6 and s["protein"]["spearman_log2fc"] > 0.7
    assert s["phosphosite"]["spearman_log2fc"] > 0.6
    d = pd.read_parquet(PROCESSED / "confirmatory_protein.parquet")
    assert d["conf_n_pairs"].min() >= 10
    # discovery tumours are excluded: the confirmatory tumours are a different set of patients
    import re
    raw = pd.read_csv(PROCESSED.parents[1] / "data" / "raw" / "confirmatory" / "CPTAC_ccRCC_global_proteome_DIA.tsv", sep="\t", nrows=2)
    disc = set(read_matrix(LAYERS["protein"][0]).columns)
    used = {c[:-2] for c in raw.columns if re.match(r"^C3[LN]-\d+_T$", c)} - disc
    assert len(used) >= 100 and d["n_tumour_conf"].iloc[0] == len(used)


def test_confirmatory_phosphosite_names_are_parsed():
    import re
    name = "NP_001347480_PABPN1_S19"
    parts = name.split("_")
    assert parts[-2] == "PABPN1" and re.match(r"^[STY]\d+$", parts[-1])
    assert ";" in "NP_001271265_EIF3J_S11;S13".split("_")[-1]       # multi-site peptides are skipped


@needs("adc_suitability.parquet")
def test_applied_evidence_is_attached_to_the_gene():
    from ccrcc.graph import GraphStore, evidence_id
    store = GraphStore()
    for ev_type in ("treatment_response", "adc_profile", "immune_profile"):
        eid = evidence_id(ev_type, "CD70")
        rec = store.get_evidence([eid])
        if not rec:
            pytest.skip("applied evidence not loaded: run `ccrcc applied-graph`")
        n = store.conn.execute("SELECT count(*) FROM claims WHERE evidence_id = ?", (eid,)).fetchone()[0]
        assert n == 1 and rec[0]["gene"] == "CD70"


def test_new_agent_is_registered_everywhere():
    import pathlib

    from ccrcc.verify import AGENT_FILES
    root = pathlib.Path(__file__).parents[1]
    assert "agent_applied.json" in AGENT_FILES
    agent = (root / ".claude" / "agents" / "agente-aplicado.md").read_text(encoding="utf-8")
    cmd = (root / ".claude" / "commands" / "informe.md").read_text(encoding="utf-8")
    assert "agente-aplicado" in cmd and "agent_applied.json" in cmd
    for tool in ("treatment_response", "subtype_profile", "adc_profile", "immune_profile", "sc_replication_check",
                 "confirmatory_replication"):
        assert tool in agent
    from ccrcc.mcp_servers import cell, omics
    for tool in ("treatment_response", "subtype_profile", "adc_profile", "immune_profile", "sc_replication_check"):
        assert hasattr(cell, tool)
    assert hasattr(omics, "confirmatory_replication")


def test_mcp_config_passes_the_replicate_folder_to_the_servers():
    import pathlib
    cfg = json.loads((pathlib.Path(__file__).parents[1] / ".mcp.json").read_text())
    for name, srv in cfg["mcpServers"].items():
        assert "CCRCC_RUN_DIR" in srv.get("env", {}), name


def test_overview_text_blocks_are_complete():
    from ccrcc.report_text import PLAIN, TOUR
    tabs = {"candidates", "omics", "cells", "structure", "discovery", "viewer", "response", "adcimmune", "validation",
            "literature", "network", "agents", "methods"}
    assert set(PLAIN) == tabs and {t[0] for t in TOUR} == tabs
    assert not any("—" in v for v in PLAIN.values())

def test_section_menu_lists_every_tab_once():
    import re
    from pathlib import Path
    from ccrcc.report_text import NAV, NAV_TEXT, TOUR
    keys = [k for _, subs in NAV for _, ks in subs for k in ks]
    assert len(keys) == len(set(keys))
    assert set(keys) == {t[0] for t in TOUR} | {"overview"} == set(NAV_TEXT)
    tpl = Path(__file__).resolve().parents[1] / "src/ccrcc/report_assets/report.html.j2"
    panels = set(re.findall(r'id="tab-([a-z]+)"', tpl.read_text(encoding="utf-8")))
    panels |= {"agents", "viewer", "adcimmune", "response"}  # panels defined in included templates
    assert set(keys) <= panels
