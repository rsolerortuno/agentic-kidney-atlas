import pytest
from ccrcc.graph import GraphStore, evidence_id
from ccrcc.paths import GRAPH_DB
from ccrcc.render import evidence_view


@pytest.mark.skipif(not GRAPH_DB.exists(), reason='requires rebuilt local graph')
def test_monomer_heuristic_is_not_labelled_a_degron():
    store=GraphStore()
    e=store.get_evidence([evidence_id('structural_context','Q16790','KAA')])[0]
    assert 'structural_context_score' in e and 'mgd_score' not in e
    assert 'exposed_loop' in e['components'] and 'degron' not in e['components']
    assert e['n_accessible_lysines']==len(e['accessible_lysines'])
    assert 'no degron or degradation prediction' in e['method']
    view=evidence_view(e)
    assert 'structural-context score' in str(view) and 'degron' not in str(view)
