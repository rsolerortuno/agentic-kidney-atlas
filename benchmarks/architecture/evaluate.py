"""Label evaluator independent of baseline A and the production claim verifier."""
import json
import random
from collections import defaultdict
from pathlib import Path

LABELS={'supported','contradicted','insufficient'}

def score(answer, expected):
    answer=answer if isinstance(answer,dict) else {}
    label=answer.get('decision'); cited=answer.get('evidence_ids',[])
    valid=label in LABELS and isinstance(cited,list) and all(isinstance(e,str) for e in cited)
    cited=set(cited) if isinstance(cited,list) and all(isinstance(e,str) for e in cited) else set()
    relevant=set(expected['relevant_evidence']); required=set(expected['required_evidence'])
    return {'valid':valid,'correct':bool(valid and label==expected['decision']),
            'critical_error':bool(valid and label in {'supported','contradicted'} and label!=expected['decision']),
            'required_recall':len(cited & required)/len(required),
            'citation_precision':len(cited & relevant)/len(cited) if cited else 0.,
            'wrong_or_irrelevant_citations':len(cited-relevant),
            'abstains':bool(valid and label=='insufficient'),
            'expected_abstention':expected['decision']=='insufficient',
            'citation_complete':bool(required<=cited and cited<=relevant)}

def bootstrap_delta(rows, a, b, field='correct', n=4000):
    """Average repetitions and variants within base; resample base-case clusters.

    Small pilot, bounded tasks. Does not model all shared-gene/batch dependencies.
    """
    grouped=defaultdict(lambda:defaultdict(list))
    for r in rows: grouped[r['base']][r['arm']].append(float(r[field]))
    delta=[]
    for v in grouped.values():
        if a in v and b in v: delta.append(sum(v[a])/len(v[a])-sum(v[b])/len(v[b]))
    if not delta: return None
    rng=random.Random(8217)
    draws=sorted(sum(rng.choices(delta,k=len(delta)))/len(delta) for _ in range(n))
    return {'difference':sum(delta)/len(delta),'interval_95':[draws[int(n*.025)],draws[int(n*.975)]], 'clusters':len(delta),'method':'paired percentile bootstrap over base-case means; exploratory, no multiplicity correction'}

def evaluate(root):
    root=Path(root); gold=json.loads((root/'frozen/gold.json').read_text(encoding='utf-8')); rows=[]
    for p in sorted((root/'runs').glob('*/answers_*.json')):
        run=p.parent.name; arm=p.stem.split('_')[-1]
        document=json.loads(p.read_text(encoding='utf-8'))
        answers={a['id']:a for a in document['answers']}
        for cid,g in gold.items():
            rows.append({'run':run,'arm':arm,'id':cid,'decision':answers.get(cid,{}).get('decision'),**{k:g[k] for k in ('base','family','variant','is_synthetic')},**score(answers.get(cid),g)})
    arms={}
    for arm in 'ABCD':
        rr=[r for r in rows if r['arm']==arm]
        if not rr: continue
        stats={'n':len(rr)}
        for k in ['correct','critical_error','required_recall','citation_precision','valid','citation_complete']:
            stats[k]=sum(r[k] for r in rr)/len(rr)
        pos=[r for r in rr if r['expected_abstention']]; neg=[r for r in rr if not r['expected_abstention']]
        stats['correct_abstention']=sum(r['abstains'] for r in pos)/len(pos) if pos else None
        stats['unnecessary_abstention']=sum(r['abstains'] for r in neg)/len(neg) if neg else None
        stats['families']={f:sum(r['correct'] for r in rr if r['family']==f)/sum(r['family']==f for r in rr) for f in {r['family'] for r in rr}}
        stats['real_accuracy']=sum(r['correct'] for r in rr if not r['is_synthetic'])/sum(not r['is_synthetic'] for r in rr)
        stats['synthetic_accuracy']=sum(r['correct'] for r in rr if r['is_synthetic'])/sum(r['is_synthetic'] for r in rr)
        # Require BOTH answers correct for a successful counterfactual pair.
        by={(r['run'],r['base'],r['variant']):r for r in rr}; pairs=[]
        for (run,base,v),r in by.items():
            if v!=0: continue
            for vv in (1,2):
                other=by.get((run,base,vv))
                if other: pairs.append(r['correct'] and other['correct'])
        stats['paired_variant_success']=sum(pairs)/len(pairs) if pairs else None
        arms[arm]=stats
    out={'scope':'frozen-packet interpretation pilot; not end-to-end MCP or therapeutic validation','arms':arms,'paired_differences':{f'{a}-{b}':bootstrap_delta(rows,a,b) for a,b in [('B','A'),('C','B'),('D','C'),('D','B')]},'rows':rows}
    (root/'results.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    return out
