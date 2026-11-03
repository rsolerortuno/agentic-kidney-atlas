"""A: explicit domain rules. Does not load gold labels or generated provenance.

This is deliberately a strong baseline for these bounded, rule-defined tasks.
It is NOT the existing target-ranking algorithm, nor a claim of general reasoning.
"""
from collections import Counter

def solve(case):
    rows=[r for r in case['evidence'] if r['kind']!='unrelated_gene']
    family=case['family']; label='insufficient'
    if family in ('replication','multiomics'):
        significant=[r.get('fc') is not None and r.get('q') is not None and r['q']<.05 for r in rows]
        if family=='replication' and all(significant):
            product=rows[0]['fc']*rows[1]['fc']
            label='supported' if product>0 else 'contradicted' if product<0 else 'insufficient'
        elif family=='multiomics':
            if any(ok and r['fc']<0 for ok,r in zip(significant,rows)): label='contradicted'
            elif all(significant) and all(r['fc']>0 for r in rows): label='supported'
    elif family=='cell_origin':
        if all(r.get('share') is not None for r in rows):
            largest=max(r['share'] for r in rows)
            winners=[r['kind'] for r in rows if r['share']==largest]
            if len(winners)==1: label='supported' if winners[0]=='RCC' else 'contradicted'
    elif family=='surface':
        a,b=rows
        if (a.get('surface_score') is not None and a['surface_score']<.7) or (b.get('available') and b.get('staining') is not None and b['staining']>0): label='contradicted'
        elif a.get('surface_score') is not None and a['surface_score']>=.7 and b.get('available') and b.get('staining')==0: label='supported'
    elif family=='treatment':
        r=next(r for r in rows if r['kind'].startswith('difference_'))
        if r.get('q') is not None and r['q']<.05: label='supported'
    elif family=='degradation':
        r=next(r for r in rows if r['kind']=='degradation_assay')
        if r.get('measured') and isinstance(r.get('positive'),bool): label='supported' if r['positive'] else 'contradicted'
    return {'id':case['id'],'decision':label,'evidence_ids':[r['id'] for r in rows], 'rationale':f'The {family} rule yields {label}. This is evidence classification, not therapeutic validation.'}

def combine(case, outputs):
    """C: majority of valid independent votes; abstain without >=3 matching votes.

    Citation union uses only agreeing votes. No LLM sees another worker output.
    """
    votes=[o for o in outputs if o and o.get('decision') in ('supported','contradicted','insufficient')]
    counts=Counter(o['decision'] for o in votes)
    top,n=counts.most_common(1)[0] if counts else ('insufficient',0)
    label=top if n>=3 else 'insufficient'
    agree=[o for o in votes if o['decision']==label]
    return {'id':case['id'],'decision':label,'evidence_ids':sorted({e for o in agree for e in o.get('evidence_ids',[]) if isinstance(e,str)}),
            'rationale':'Independent specialist vote counts: '+str(dict(counts))+'. '+(' '.join(o.get('rationale','') for o in agree))[:1800],
            'votes':dict(counts),'worker_failures':5-len(votes)}
