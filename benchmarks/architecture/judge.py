"""Secondary LLM assessment. An explicit separate command, never the main score.

Two blind presentation orders with the same model measure order sensitivity;
they are NOT two independent experts or independent model families.
"""
import argparse
import json
import random
from pathlib import Path
from .build import HERE, save
from .run import call, fingerprint, verify_frozen

DIMENSIONS=['accuracy','integration','contradictions','calibration','next_experiment']
RUBRIC='''Score each anonymous answer against the supplied evidence, not stylistic polish or length.
Use five dimensions, each 0,1,2: accuracy (major error / minor imprecision / correct);
integration (misses required layers / lists them / relates them);
contradictions (misses relevant conflict / notices / interprets correctly, including missing vs negative);
calibration (overclaims / partly qualified / conclusion proportional to evidence);
next_experiment (irrelevant or absent / generic / concrete and resolves uncertainty).
For each score give an exact quote or identify the omission supporting it. Flag critical biological overclaims separately. Deterministic answers may be terse: do not penalize style, but do not invent reasoning absent from the answer.
Return {"answers":[{"id":"anonymous answer id","scores":{"accuracy":0,"integration":0,"contradictions":0,"calibration":0,"next_experiment":0},"justification":"Evidence-specific explanation with quotes","critical_error":false}]}.
These are bounded classification tasks; next-experiment quality is secondary and may not be fully supported by these packets.'''

def parse(text,ids):
    if isinstance(text,str):
        text=text.strip().removeprefix('```json').removesuffix('```').strip()
        text=json.loads(text)
    aa=text['answers']
    if len(aa)!=len(ids) or {a['id'] for a in aa}!=set(ids): raise ValueError('Invalid judge IDs')
    for a in aa:
        if set(a['scores'])!=set(DIMENSIONS): raise ValueError('Invalid rubric')
        if any(type(s)!=int or s not in (0,1,2) for s in a['scores'].values()): raise ValueError('Invalid judge scores')
        if not isinstance(a.get('justification'),str) or type(a.get('critical_error'))!=bool: raise ValueError('Missing judge explanation')
    return aa

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--execute',action='store_true'); ap.add_argument('--claude',type=Path); a=ap.parse_args()
    if not a.execute:
        print('Dry run: 12 native cases, four arms, first replicate only, two reversed anonymous orders; 24 judge calls. No model called.'); return
    if not a.claude or not a.claude.is_file(): ap.error('Explicit --claude required')
    verify_frozen(); h=fingerprint(); cases=json.loads((HERE/'frozen/cases.json').read_text())
    source={arm:json.loads((HERE/'runs'/('deterministic' if arm=='A' else 'rep1')/f'answers_{arm}.json').read_text()) for arm in 'ABCD'}
    if any(x['study_hash']!=h for x in source.values()): raise ValueError('Protocol drift')
    mapping={}; rng=random.Random(7842)
    for case in [c for c in cases if c['id'].endswith('v0')]:
        items=[]
        for arm in 'ABCD':
            answer=next((r for r in source[arm]['answers'] if r['id']==case['id']),None)
            if answer is None: raise ValueError('Missing answer; do not silently omit a failed architecture')
            anon=f'answer_{rng.getrandbits(64):016x}'
            mapping[anon]={'arm':arm,'case':case['id']}
            items.append({'id':anon,'decision':answer['decision'],'evidence_ids':answer['evidence_ids'],'rationale':answer['rationale']})
        rng.shuffle(items)
        for order,sequence in enumerate([items,list(reversed(items))]):
            prompt=RUBRIC+'\nCASE:\n'+json.dumps(case)+'\nANONYMOUS ANSWERS:\n'+json.dumps(sequence)
            call(a.claude,HERE/'judges'/f'{case["id"]}_order{order}.json',prompt,[i['id'] for i in items],1.,4000,h,parser=parse)
        save(HERE/'judges/identity_map.json',mapping)

if __name__=='__main__': main()
