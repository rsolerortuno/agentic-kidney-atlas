"""Offline tests; none launches a model or depends on downloaded omics files."""
import json
import pytest
from benchmarks.architecture.baseline import solve,combine
from benchmarks.architecture.evaluate import score,bootstrap_delta
from benchmarks.architecture.run import parse_answers,batches,verify_frozen
from benchmarks.architecture.judge import parse as judge_parse
from benchmarks.architecture.build import HERE

def test_frozen_files_have_not_drifted():
    assert verify_frozen()['n_cases']==36

def test_no_same_seed_variants_in_one_model_context():
    cases=json.loads((HERE/'frozen/cases.json').read_text())
    for rep in range(1,4):
        bb=batches(cases,rep)
        assert len(bb)==6
        assert {c['id'] for b in bb for c in b}=={c['id'] for c in cases}
        assert all(len({c['id'].split('v')[0] for c in b})==len(b) for b in bb)

def test_label_sabotage_and_false_abstention():
    gold={'decision':'insufficient','required_evidence':['e1'],'relevant_evidence':['e1','e2']}
    assert score({'decision':'supported','evidence_ids':['e1']},gold)['critical_error']
    assert score({'decision':'insufficient','evidence_ids':['e1']},gold)['correct']
    assert not score(None,gold)['correct']
    bad=score({'decision':'insufficient','evidence_ids':['WRONG_GENE']},gold)
    assert bad['required_recall']==0 and bad['citation_precision']==0
    assert not bad['citation_complete']

def test_rule_does_not_treat_nonsignificance_as_negative():
    case={'id':'test','family':'treatment','evidence':[{'id':'e1','kind':'difference_x','q':.9},{'id':'e2','kind':'within_arm','q':.0001}]}
    assert solve(case)['decision']=='insufficient'
    case['evidence'][0]['q']=.001
    assert solve(case)['decision']=='supported'

def test_missing_staining_not_safe():
    case={'id':'test','family':'surface','evidence':[{'id':'e1','kind':'topology','surface_score':.9},{'id':'e2','kind':'normal_tissue','available':False,'staining':None}]}
    assert solve(case)['decision']=='insufficient'

def test_vote_requires_three_and_does_not_copy_dissenting_citations():
    yes={'decision':'supported','evidence_ids':['YES'],'rationale':'supported'}
    no={'decision':'contradicted','evidence_ids':['NO'],'rationale':'contradicted'}
    assert combine({'id':'test'},[yes,yes,no,no,None])['decision']=='insufficient'
    result=combine({'id':'test'},[yes,yes,yes,no,no])
    assert result['decision']=='supported' and result['evidence_ids']==['YES']

def test_parser_rejects_missing_duplicate_and_foreign_cases():
    a={'id':'x','decision':'supported','evidence_ids':['e'],'rationale':'why'}
    assert parse_answers({'answers':[a]},['x'])==[a]
    for rows in [[],[a,a],[dict(a,id='other')]]:
        with pytest.raises(ValueError): parse_answers({'answers':rows},['x'])

def test_cluster_bootstrap_averages_repetitions_not_pseudoreplication():
    rows=[{'base':b,'arm':a,'correct':a=='D'} for b in ['x','y'] for a in ['B','D'] for _ in range(3)]
    out=bootstrap_delta(rows,'D','B',n=100)
    assert out['clusters']==2 and out['difference']==1 and out['interval_95']==[1,1]

def test_judge_invalid_score_rejected():
    with pytest.raises(ValueError): judge_parse({'answers':[{'id':'x','scores':{'accuracy':3}}]},['x'])

def test_model_packets_exclude_gene_names_sources_and_gold():
    cases=json.loads((HERE/'frozen/cases.json').read_text())
    assert all(set(c)=={'id','family','entity','question','evidence'} for c in cases)
    assert all(c['entity'].startswith('ENTITY_') for c in cases)
