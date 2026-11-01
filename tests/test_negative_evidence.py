"""A measured unfavorable result must remain citable, not disappear as missing data."""
import pytest
from ccrcc.graph import GraphStore, evidence_id
from ccrcc.paths import GRAPH_DB


@pytest.mark.skipif(not GRAPH_DB.exists(), reason="requires rebuilt graph")
def test_measured_non_dependency_and_null_survival_are_citable():
    store = GraphStore()
    ids = [evidence_id(kind, 'CA9') for kind in ('depmap_dependency','tcga_survival')]
    records = {e['type']:e for e in store.get_evidence(ids)}
    assert records['depmap_dependency']['mean_effect_ccrcc'] > 0
    assert records['tcga_survival']['padj'] > .05
    assert records['tcga_survival']['direction'] == 'association not established'
    for eid in ids:
        predicates = {r[0] for r in store.conn.execute('SELECT predicate FROM claims WHERE evidence_id=?',(eid,))}
        assert predicates & {'dependency_assessed_in','survival_assessed_in'}
        assert not predicates & {'selective_dependency_in','prognostic_in'}
