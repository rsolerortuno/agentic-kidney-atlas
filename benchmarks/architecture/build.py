"""Build auditable, de-identified evidence packets from a read-only graph snapshot.

Expected answers are explicit experimental assignments, not calls to baseline A.
The original source record is retained outside model-visible packets.
"""
import argparse
import copy
import hashlib
import json
import sqlite3
from pathlib import Path

HERE = Path(__file__).parent
OUT = HERE / 'frozen'

def save(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False)+'\n', encoding='utf-8', newline='\n')

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

QUESTIONS = {
 'replication': 'Claim: the discovery protein change replicates in direction in the independent cohort. Require q<0.05 in both cohorts and matching nonzero signs. An opposite significant sign contradicts the claim; a nonsignificant or unavailable result is insufficient.',
 'multiomics': 'Claim: both RNA and protein are significantly increased. Require positive log2 fold changes and q<0.05 in both layers. A significant decrease in either layer contradicts the claim; missing or nonsignificant measurements are insufficient.',
 'cell_origin': 'Claim: the largest share of this gene\'s transcripts in the sampled tumour compartment comes from RCC cells. Compare transcript shares, not the percentage of cells expressing the gene. A larger non-RCC share contradicts the claim; missing shares or tied maxima are insufficient. Do not infer expression in every tumour cell.',
 'surface': 'Claim: this candidate passes this exploratory surface-screen rule: surface_score>=0.7 AND available vital-organ staining with level=0. This is a screening rule, not demonstrated safety. A score below 0.7 or measured vital-organ staining above 0 contradicts eligibility; missing required data are insufficient.',
 'treatment': 'Claim: there is statistical evidence of a gene-by-treatment interaction at adjusted q<0.05. A within-arm association alone does not establish this. A reported interaction q>=0.05 or missing interaction is insufficient (not proof of no interaction).',
 'degradation': 'Claim: compound-dependent degradation has been experimentally demonstrated in this specified system. Require a measured degradation assay with positive result. A measured negative assay contradicts the claim; a pocket, ligase expression or recruitment alone is insufficient. Do not infer ccRCC efficacy.'
}

def build(db, dest=OUT):
    if (dest/'manifest.json').exists():
        raise FileExistsError('Frozen suite already exists; use a new destination to prevent retrospective edits.')
    con=sqlite3.connect(f'file:{Path(db).resolve().as_posix()}?mode=ro', uri=True)
    def source(kind, gene, layer=None):
        rows=con.execute('SELECT id,type,source,props FROM evidence WHERE type=? AND json_extract(props,\'$.gene\')=?', (kind,gene)).fetchall()
        for eid,typ,src,p in rows:
            d=json.loads(p)
            if layer is None or d.get('layer')==layer:
                return {'id':eid,'type':typ,'source':src,'props':d}
        raise ValueError(f'Missing source: {kind} {gene} {layer}')
    seeds=[]
    for gene in ['NNMT','SLC16A3']:
        r=source('confirmatory_replication',gene); p=r['props']
        assert p['protein_discovery_log2fc']>0 and p['protein_conf_log2fc']>0
        assert p['protein_discovery_padj']<.05 and p['protein_conf_padj']<.05
        seeds.append(('replication',gene,[r],[{'kind':'discovery','fc':p['protein_discovery_log2fc'],'q':p['protein_discovery_padj']},{'kind':'confirmation','fc':p['protein_conf_log2fc'],'q':p['protein_conf_padj']}],'supported'))
    for gene in ['CA9','NNMT']:
        rs=[source('omics_de',gene,l) for l in ['rna','protein']]
        assert all(r['props']['log2fc']>0 and r['props']['padj']<.05 for r in rs)
        seeds.append(('multiomics',gene,rs,[{'kind':r['props']['layer'],'fc':r['props']['log2fc'],'q':r['props']['padj']} for r in rs],'supported'))
    for gene in ['CA9','COL23A1']:
        r=source('sc_summary',gene); groups=r['props']['groups']
        groups=sorted(groups,key=lambda g:g['share_of_tumor_transcripts'],reverse=True)
        assert groups[0]['cell_group']=='RCC'
        seeds.append(('cell_origin',gene,[r],[{'kind':g['cell_group'],'share':g['share_of_tumor_transcripts'],'pct_expressing':g['pct_expressing']} for g in groups],'supported'))
    for gene in ['CA9','CD70']:
        r=source('adc_profile',gene); p=r['props']
        # Record the observed result explicitly; never claim the screen proves safety.
        if p['surface_score']<.7 or (p['staining_data_available'] and p['vital_organ_staining_level']>0): label='contradicted'
        elif p['staining_data_available'] and p['vital_organ_staining_level']==0: label='supported'
        else: label='insufficient'
        seeds.append(('surface',gene,[r],[{'kind':'topology','surface_score':p['surface_score']},{'kind':'normal_tissue','available':p['staining_data_available'],'staining':p.get('vital_organ_staining_level')}],label))
    for gene in ['SULF1','CA9']:
        r=source('treatment_response',gene); arms=r['props']['arms']
        k=next(k for k in arms if k.startswith('difference_'))
        assert arms[k]['padj']>=.05
        seeds.append(('treatment',gene,[r],[{'kind':k,'q':arms[k]['padj'],'log_hr':arms[k]['log_hr_per_sd']},{'kind':'within_arm','q':next(v['padj'] for k,v in arms.items() if not k.startswith('difference_'))}],'insufficient'))
    for gene in ['CA9','EPAS1']:
        r=source('structural_context',gene)
        seeds.append(('degradation',gene,[r],[{'kind':'predicted_monomer','pocket_volume':r['props']['volume_a3'],'mean_plddt':r['props']['mean_plddt']},{'kind':'degradation_assay','measured':False,'positive':None}],'insufficient'))
    cases=[]; gold={}; provenance={}
    for i,(family,gene,sources,records,real_label) in enumerate(seeds):
        base=f'case{i+1:02d}'
        for variant in range(3):
            rec=copy.deepcopy(records); label=real_label; mutation='none; normalized real record projection'
            if variant==1:
                if family in ['replication','multiomics']:
                    rec[1]['fc']=-abs(rec[1]['fc']); label='contradicted'; mutation='second fold-change sign reversed'
                elif family=='cell_origin':
                    rec[0]['share'],rec[1]['share']=rec[1]['share'],rec[0]['share']; label='contradicted'; mutation='RCC and leading other group shares exchanged'
                elif family=='surface':
                    rec[0]['surface_score']=.9; rec[1].update(available=True,staining=3); label='contradicted'; mutation='surface score .9 and observed vital-organ staining 3'
                elif family=='treatment':
                    rec[0]['q']=.001; label='supported'; mutation='interaction q set to .001; synthetic, not a clinical finding'
                else:
                    rec[1].update(measured=True,positive=True); label='supported'; mutation='synthetic positive degradation assay; not a real experiment'
            elif variant==2:
                if family in ['replication','multiomics']:
                    rec[1].update(fc=None,q=None); label='insufficient'; mutation='second layer withheld'
                elif family=='cell_origin':
                    for r in rec: r['share']=None
                    label='insufficient'; mutation='transcript shares withheld; expressing-cell percentages retained'
                elif family=='surface':
                    rec[0]['surface_score']=.9; rec[1].update(available=False,staining=None); label='insufficient'; mutation='surface score .9; staining withheld'
                elif family=='treatment':
                    rec[0]['q']=None; rec[1]['q']=.001; label='insufficient'; mutation='interaction withheld; significant within-arm association added'
                else:
                    rec[1].update(measured=True,positive=False); label='contradicted'; mutation='synthetic negative degradation assay; not a real experiment'
            cid=f'{base}v{variant}'
            for n,r in enumerate(rec): r['id']=f'{cid}:e{n+1}'
            # Every variant additionally includes an irrelevant record to measure citation discipline.
            rec.append({'id':f'{cid}:d1','kind':'unrelated_gene','entity':'OTHER','rna_fc':8,'q':.0001})
            c={'id':cid,'family':family,'entity':f'ENTITY_{i+1:02d}','question':QUESTIONS[family], 'evidence':rec}
            if family=='cell_origin': required=[r['id'] for r in rec if r['kind']!='unrelated_gene']
            elif family in ['treatment','degradation']: required=[rec[0 if family=='treatment' else 1]['id']]
            else: required=[rec[0]['id'],rec[1]['id']]
            relevant=[r['id'] for r in rec if r['kind']!='unrelated_gene']
            cases.append(c)
            gold[cid]={'decision':label,'required_evidence':required,'relevant_evidence':relevant,'base':base,'family':family,'variant':variant,'is_synthetic':variant!=0}
            provenance[cid]={'gene':gene,'sources':sources,'mutation':mutation,'normalized_evidence':rec}
    save(dest/'cases.json',cases); save(dest/'gold.json',gold); save(dest/'provenance.json',provenance)
    save(dest/'manifest.json',{'schema':1,'source_db_sha256':sha(Path(db)),'files':{p.name:sha(p) for p in dest.glob('*.json')},'n_cases':len(cases),'n_base_cases':len(seeds),'split':'evaluation; never used for prompt tuning','scope':'closed-world, masked-gene interpretation; not therapeutic validity or end-to-end MCP retrieval'})
    return cases

if __name__=='__main__':
    ap=argparse.ArgumentParser(); ap.add_argument('--db',required=True); ap.add_argument('--out',type=Path,default=OUT)
    a=ap.parse_args(); print(f'Built {len(build(a.db,a.out))} scenarios')
