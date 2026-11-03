"""M4 checks: every relation goes through a resolvable Evidence record."""
import sqlite3

import pytest

from ccrcc.graph import GraphStore, evidence_id
from ccrcc.paths import GRAPH_DB


def test_claim_requires_existing_evidence(tmp_path):
    s = GraphStore(tmp_path / "g.sqlite")
    s.add_node("HGNC:1", "Gene", "A")
    s.add_node("MONDO:1", "Disease", "D")
    with pytest.raises(sqlite3.IntegrityError):
        s.add_claim("HGNC:1", "upregulated_in", "MONDO:1", "EV:missing")
    ev = s.add_evidence(evidence_id("t", 1), "omics_de", "DATASET:X", log2fc=1.0)
    s.add_claim("HGNC:1", "upregulated_in", "MONDO:1", ev)
    s.conn.commit()
    assert s.neighbors("HGNC:1")[0]["evidence_ids"] == [ev]
    assert s.check_evidence_ids([ev, "EV:nope"]) == {"valid": [ev], "invalid": ["EV:nope"]}


def test_reloading_evidence_replaces_old_values(tmp_path):
    s = GraphStore(tmp_path / "g.sqlite")
    ev = s.add_evidence(evidence_id("t", 1), "omics_de", "DATASET:X", log2fc=1.0)
    s.add_evidence(ev, "omics_de", "DATASET:X", log2fc=2.0)
    assert s.get_evidence([ev])[0]["log2fc"] == 2.0


@pytest.mark.skipif(not GRAPH_DB.exists(), reason="run `ccrcc graph` first")
def test_built_graph_is_fully_evidenced():
    st = GraphStore().stats()
    assert st["claims_without_resolvable_evidence"] == 0
    for t in ("Gene", "Protein", "Pathway", "Drug", "Disease", "Dataset", "Publication", "Chunk"):
        assert st["nodes_by_type"].get(t, 0) > 0, t
    assert st["evidence_by_type"]["literature_relation"] > 0
    ca9 = GraphStore().resolve("CA9")
    preds = {n["predicate"] for n in GraphStore().neighbors(ca9["id"], limit=500)}
    assert {"upregulated_in", "targets", "associated_with"} <= preds
