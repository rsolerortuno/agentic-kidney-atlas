"""Small manually reviewed abstract annotations, tied to immutable source text.

Reported assay numbers are not newly computed measurements or independent
replications. Context and units remain part of each metric's name/annotation.
"""
import hashlib
import json
from pathlib import Path
import pandas as pd
from .graph import GraphStore, evidence_id
from .paths import LITERATURE

CATALOG = Path(__file__).parent / 'reference_data' / 'literature_review.json'


def validate(entry, abstract):
    if hashlib.sha256(abstract.encode()).hexdigest() != entry['abstract_sha256']:
        raise ValueError(f"Abstract changed for PMID {entry['pmid']}; review required")
    for quote in entry['source_quotes']:
        if quote not in abstract:
            raise ValueError('Reviewed source quotation is absent')


def load(store=None, publications=None):
    store = store or GraphStore()
    publications = publications if publications is not None else pd.read_parquet(LITERATURE / 'publications.parquet')
    lookup={str(r.pmid):str(r.abstract) for r in publications.itertuples()}
    ids=[]
    with store.batch():
        for entry in json.loads(CATALOG.read_text(encoding='utf-8')):
            abstract=lookup.get(str(entry['pmid']))
            if abstract is None:
                raise ValueError(f"Missing reviewed publication PMID {entry['pmid']}")
            validate(entry,abstract)
            gene=store.resolve(entry['gene'])
            if not gene or gene['type']!='Gene':
                raise ValueError(f"Missing reviewed gene {entry['gene']}")
            pub=f"PMID:{entry['pmid']}"
            store.add_node(pub,'Publication',entry['title'])
            ev=evidence_id('literature_review',entry['gene'],entry['pmid'],entry['abstract_sha256'])
            store.add_evidence(ev,'literature_review',pub,gene=entry['gene'],pmid=entry['pmid'],
                               url=f"https://pubmed.ncbi.nlm.nih.gov/{entry['pmid']}/",
                               abstract_sha256=entry['abstract_sha256'],source_quotes=entry['source_quotes'],
                               interpretation=entry['interpretation'],reported_metrics=entry.get('reported_metrics',{}))
            store.add_claim(gene['id'],'reviewed_in',pub,ev);ids.append(ev)
    return ids
