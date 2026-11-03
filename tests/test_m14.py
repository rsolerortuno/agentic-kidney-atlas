"""M14 checks: reproducible spatial permutations, stale-narrative guard, discovery v2 gates."""
import json

import numpy as np
import pandas as pd
import pytest

from ccrcc import provenance
from ccrcc.graph import GraphStore
from ccrcc.paths import PROCESSED
from ccrcc.spatial import gene_rng, neighborhood_enrichment
from ccrcc.verify import Verifier


def test_gene_rng_independent_of_other_genes():
    a = gene_rng(7, "CA9").integers(0, 10**9, 5)
    _ = gene_rng(7, "KDR").integers(0, 10**9, 50)     # analysing another gene first must not matter
    b = gene_rng(7, "CA9").integers(0, 10**9, 5)
    assert (a == b).all()
    assert not (gene_rng(7, "CA9").integers(0, 10**9, 5) == gene_rng(7, "KDR").integers(0, 10**9, 5)).all()


def test_neighbourhood_enrichment_reproducible():
    rng = np.random.default_rng(0)
    n = 400
    NP = rng.random((n, 3))
    sections = np.repeat([0, 1], n // 2)
    high = rng.random(n) < 0.1
    has_nb = np.ones(n, bool)
    e1 = neighborhood_enrichment(high, sections, NP, has_nb, gene_rng(7, "X"), 200)
    e2 = neighborhood_enrichment(high, sections, NP, has_nb, gene_rng(7, "X"), 200)
    assert np.allclose(e1["p"], e2["p"]) and np.allclose(e1["observed"], e2["observed"])


def test_stale_evidence_is_detected(tmp_path):
    store = GraphStore(tmp_path / "g.sqlite")
    with store.batch():
        store.add_node("HGNC:1", "Gene", "GENEA")
        store.add_node("MONDO:0005005", "Disease", "ccRCC")
        store.add_evidence("EV:t:1", "omics_de", "test", "2026-01-01", gene="GENEA", log2fc=1.5)
        store.add_claim("HGNC:1", "upregulated_in", "MONDO:0005005", "EV:t:1")
    basis = {"started_utc": "x", "config_sha256": "x", "evidence": provenance.fingerprints(store, ["EV:t:1"])}
    v = Verifier(store=store, basis=basis, semantic_check=False)
    stmt = {"text": "GENEA log2FC 1.5", "evidence_ids": ["EV:t:1"]}
    assert v.check_statement(stmt)[0] is not None
    # rebuilding on another day changes only the date: not stale
    store.conn.execute("UPDATE evidence SET date = '2026-02-02' WHERE id = 'EV:t:1'")
    assert v.check_statement(stmt)[0] is not None
    # a changed value is stale
    store.conn.execute("""UPDATE evidence SET props = '{"gene": "GENEA", "log2fc": 2.5}' WHERE id = 'EV:t:1'""")
    ok, log = v.check_statement(stmt)
    assert ok is None and log["stale_evidence"] == ["EV:t:1"]


def test_id_missing_from_basis_fails_closed(tmp_path):
    store = GraphStore(tmp_path / "g.sqlite")
    with store.batch():
        store.add_node("HGNC:1", "Gene", "GENEA")
        store.add_node("MONDO:0005005", "Disease", "ccRCC")
        store.add_evidence("EV:t:1", "omics_de", "test", None, gene="GENEA", log2fc=1.5)
        store.add_evidence("EV:t:2", "omics_de", "test", None, gene="GENEA", log2fc=1.5)
        store.add_claim("HGNC:1", "upregulated_in", "MONDO:0005005", "EV:t:1")
        store.add_claim("HGNC:1", "upregulated_in", "MONDO:0005005", "EV:t:2")
    basis = {"started_utc": "x", "config_sha256": "x", "evidence": provenance.fingerprints(store, ["EV:t:1"])}
    v = Verifier(store=store, basis=basis, semantic_check=False)
    ok, log = v.check_statement({"text": "GENEA is up", "evidence_ids": ["EV:t:2"]})
    assert ok is None and log["stale_evidence"] == ["EV:t:2"]


def test_extend_basis_never_refreshes(tmp_path):
    store = GraphStore(tmp_path / "g.sqlite")
    with store.batch():
        store.add_evidence("EV:t:1", "x", "test", None, v=1)
        store.add_evidence("EV:t:2", "x", "test", None, v=2)
    provenance.start_basis(tmp_path, store, ["EV:t:1"])
    first = provenance.load_basis(tmp_path)["evidence"]["EV:t:1"]
    store.conn.execute("""UPDATE evidence SET props = '{"v": 9}' WHERE id = 'EV:t:1'""")
    b = provenance.extend_basis(tmp_path, store, ["EV:t:1", "EV:t:2"])
    assert b["evidence"]["EV:t:1"] == first and "EV:t:2" in b["evidence"]


@pytest.mark.skipif(not (PROCESSED / "discovery_candidates.parquet").exists(), reason="run `ccrcc run-all` first")
def test_discovery_v2_gates():
    d = pd.read_parquet(PROCESSED / "discovery_candidates.parquet")
    if "replicated" not in d.columns:
        pytest.skip("discovery not rebuilt with M14")
    fic = d[d.status == "first-in-class opportunity"]
    assert fic.replicated.all()
    assert fic[fic.track.isin(["tumor", "stroma"])].purity_robust.all()
    assert (d.rank_in_track >= 1).all() and d.groupby("track").rank_in_track.min().eq(1).all()
    summ = json.loads((PROCESSED / "discovery_summary.json").read_text())
    assert set(summ["top_by_track"]) == set(d.track.unique())
