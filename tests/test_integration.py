"""M6 checks: deterministic convergence and the evidence verifier."""
import pytest

from ccrcc.graph import GraphStore, evidence_id
from ccrcc.paths import GRAPH_DB, PROCESSED
from ccrcc.verify import Verifier

pytestmark = pytest.mark.skipif(not GRAPH_DB.exists() or not (PROCESSED / "gene_universe.parquet").exists(),
                                reason="build the graph first")


def test_convergence_is_deterministic_and_bounded():
    from ccrcc.config import load
    from ccrcc.integrate import compute
    a, b = compute(), compute()
    assert a["gene"].head(50).tolist() == b["gene"].head(50).tolist()
    sw = load()["convergence"]["stream_weights"]
    assert a["score"].max() <= sum(sw.values()) + 1e-9
    for k in sw:
        assert a[f"c_{k}"].between(0, 1).all(), k


def test_verifier_requires_existing_and_gene_linked_evidence():
    v = Verifier()
    ca9_de = evidence_id("omics_de", "protein", "CA9")
    vhl_cnv = evidence_id("omics_cnv", "VHL")
    assert GraphStore().check_evidence_ids([ca9_de, vhl_cnv])["invalid"] == []
    ok, _ = v.check_statement({"gene": "CA9", "statement": "x", "evidence_ids": [ca9_de]})
    assert ok and ok["evidence_ids"] == [ca9_de]
    # real id but about another gene -> rejected
    bad, log = v.check_statement({"gene": "CA9", "statement": "x", "evidence_ids": [vhl_cnv]})
    assert bad is None and log["not_linked_to_gene"] == [vhl_cnv]
    # invented id -> rejected
    bad, log = v.check_statement({"gene": "CA9", "statement": "x", "evidence_ids": ["EV:omics_de:000000000000"]})
    assert bad is None and log["nonexistent"]
    # mixed -> invalid ids stripped
    ok, _ = v.check_statement({"gene": "CA9", "statement": "x", "evidence_ids": [ca9_de, vhl_cnv]})
    assert ok["evidence_ids"] == [ca9_de]


def test_legacy_symbol_resolves_for_verification():
    v = Verifier()
    ev = evidence_id("omics_de", "protein", "COXFA4L2")
    ok, _ = v.check_statement({"gene": "NDUFA4L2", "statement": "x", "evidence_ids": [ev]})
    assert ok is not None
