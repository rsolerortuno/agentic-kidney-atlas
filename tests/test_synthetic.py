"""End-to-end checks on small synthetic inputs, so CI exercises the pipeline without the (git-ignored) source data.

Covered: verification of a findings file against a toy graph (linkage, numbers, direction, stale basis), the verifier
sabotage test, leave-one-stream-out, the no-agent baseline list, the pooled-scale treatment interaction and HGNC
re-keying with a toy HGNC table.
"""
import json

import numpy as np
import pandas as pd
import pytest

from ccrcc import provenance, verify
from ccrcc.graph import GraphStore


@pytest.fixture
def toy(tmp_path, monkeypatch):
    store = GraphStore(tmp_path / "g.sqlite")
    with store.batch():
        store.add_node("MONDO:0005005", "Disease", "ccRCC")
        for i, (g, lfc) in enumerate((("GENEA", 2.5), ("GENEB", -1.75), ("GENEC", 3.25)), start=1):
            store.add_node(f"HGNC:{i}", "Gene", g)
            store.add_evidence(f"EV:de:{g}", "omics_de", "DATASET:TOY", None, gene=g, layer="protein", log2fc=lfc,
                               padj=1e-6, significant=True)
            store.add_claim(f"HGNC:{i}", "upregulated_in" if lfc > 0 else "downregulated_in", "MONDO:0005005", f"EV:de:{g}")
    monkeypatch.setattr(verify, "GraphStore", lambda: store)
    monkeypatch.setattr(verify, "resolve_gene", lambda g: g.upper() if isinstance(g, str) else None)
    run = tmp_path / "run"
    run.mkdir()
    findings = {"agent_findings": {"omics": {"findings": [
        {"gene": "GENEA", "statement": "GENEA protein is up in tumour (log2FC 2.5).", "evidence_ids": ["EV:de:GENEA"]},
        {"gene": "GENEB", "statement": "GENEB protein is lower in tumour (log2FC -1.75).", "evidence_ids": ["EV:de:GENEB"]},
        {"gene": "GENEC", "statement": "GENEC protein is higher in tumour (log2FC 3.25).", "evidence_ids": ["EV:de:GENEC"]},
    ]}}, "narrative": {"summary": [{"text": "GENEA is elevated in tumours while other genes vary.",
                                    "evidence_ids": ["EV:de:GENEA"]}], "sections": {}},
        "targets": []}
    (run / "findings.json").write_text(json.dumps(findings), encoding="utf-8")
    provenance.start_basis(run, store, ["EV:de:GENEA", "EV:de:GENEB", "EV:de:GENEC"])
    return store, run, findings


def test_verify_final_on_toy_run(toy):
    store, run, findings = toy
    rep = verify.verify_final(run / "findings.json")
    assert rep["statements"] == 4 and rep["unsupported"] == 0, rep["detail"]


def test_verify_final_catches_wrong_direction_and_unknown_basis_id(toy):
    store, run, findings = toy
    findings["agent_findings"]["omics"]["findings"][0]["statement"] = "GENEA protein is down in tumour."
    with store.batch():
        store.add_evidence("EV:de:new", "omics_de", "DATASET:TOY", None, gene="GENEC", layer="rna", log2fc=1.0,
                           padj=0.01, significant=True)
        store.add_claim("HGNC:3", "upregulated_in", "MONDO:0005005", "EV:de:new")
    findings["agent_findings"]["omics"]["findings"][2]["evidence_ids"] = ["EV:de:new"]     # never in basis.json
    (run / "findings.json").write_text(json.dumps(findings), encoding="utf-8")
    rep = verify.verify_final(run / "findings.json")
    assert rep["unsupported"] == 2 and rep["stale"] == 1 and rep["semantic_failures"] == 2


def test_sabotage_reports_detection_rates(toy, tmp_path):
    from ccrcc import sabotage
    store, run, _ = toy
    res = sabotage.run(run, out_path=tmp_path / "sab.json")
    c = res["corruptions"]
    assert res["statements_passing"] == 4
    assert c["number"]["applied"] == 3 and c["number"]["rate"] == 1.0
    assert c["gene"]["rate"] == 1.0 and c["evidence"]["rate"] == 1.0
    assert c["direction"]["applied_with_de_evidence"] >= 3 and c["direction"]["rate_with_de_evidence"] == 1.0


def test_leave_one_stream_out_flags_single_stream_genes():
    from ccrcc.scores import leave_one_out
    rng = np.random.default_rng(0)
    t = pd.DataFrame({"gene": [f"G{i}" for i in range(100)], "c_a": rng.random(100) * 0.2, "c_b": rng.random(100) * 0.2})
    t.loc[0, "c_a"] = 1.0                       # G0 leads only because of stream a
    t.loc[1, ["c_a", "c_b"]] = 0.8              # G1 is supported by both
    out = leave_one_out(t, {"a": 1.0, "b": 1.0}, top=5).set_index("gene")
    assert out.loc["G0", "leaves_top_without"] == ["a"]
    assert not out.loc["G1", "single_stream_dependent"]


def test_no_agent_list_alternates_sources(monkeypatch, tmp_path):
    from ccrcc import baselines
    monkeypatch.setattr(baselines, "PROCESSED", tmp_path)
    pd.DataFrame({"gene": ["D1", "D2"], "track": ["tumor", "tumor"], "rank_in_track": [1, 2]}).to_parquet(
        tmp_path / "discovery_candidates.parquet")
    s = pd.DataFrame({"gene": ["T1", "T2", "B1", "B2"], "score_target": [4, 3, 0, 0], "score_biomarker": [0, 0, 4, 3]})
    assert baselines.no_agent_list(5, s) == ["T1", "B1", "D1", "T2", "B2"]


def test_pooled_interaction_detects_an_arm_specific_gene():
    from ccrcc.treatment import _interaction
    rng = np.random.default_rng(1)
    n = 200
    arm = np.where(np.arange(n) < n // 2, "io", "other")
    g1 = rng.normal(size=n)                     # protective on io only
    g2 = rng.normal(size=n)                     # no effect
    hazard = np.exp(np.where(arm == "io", -1.0 * g1, 0.0))
    t = rng.exponential(1 / hazard)
    clin = pd.DataFrame({"arm": arm, "pfs_time": t, "pfs_event": 1, "responder": np.nan},
                        index=[f"P{i}" for i in range(n)])
    x = pd.DataFrame([g1 + 5, g2 + 5], index=["G1", "G2"], columns=clin.index)
    out = _interaction(x, clin, "io", "other", "toy", "io_vs_other").set_index("gene")
    assert out.loc["G1", "cox_log_hr"] < -0.5 and out.loc["G1", "cox_pvalue"] < 1e-3
    assert out.loc["G2", "cox_pvalue"] > 0.01
