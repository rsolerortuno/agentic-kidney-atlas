"""M4: build the knowledge graph from omics results (M2), literature (M3) and the external cache (M5)."""
from __future__ import annotations

import json
import re
from datetime import date

import pandas as pd

from .config import load
from .graph import GraphStore, evidence_id
from .literature_entities import is_ccrcc, resolve_gene
from .paths import LITERATURE, PROCESSED, RAW

DISEASE = "MONDO:0005005"
TODAY = date.today().isoformat()
DATASETS = {
    "rna": ("DATASET:CPTAC-CCRCC-RNA", "CPTAC ccRCC discovery RNA-seq (log2 FPKM)"),
    "protein": ("DATASET:CPTAC-CCRCC-PROTEOME", "CPTAC ccRCC discovery proteome (TMT)"),
    "phospho": ("DATASET:CPTAC-CCRCC-PHOSPHO", "CPTAC ccRCC discovery phosphoproteome (gene level, TMT)"),
    "cnv": ("DATASET:CPTAC-CCRCC-CNV", "CPTAC ccRCC discovery copy number (gene level)"),
}


def gene_node_id(hgnc_id: str) -> str:
    return hgnc_id if hgnc_id.startswith("HGNC:") else f"HGNC:{hgnc_id}"


def _slug(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", text.lower()).strip("_")[:60]


def build(store: GraphStore | None = None) -> dict:
    cfg = load()
    store = store or GraphStore()
    store.reset()
    hg = pd.read_csv(RAW / "hgnc_complete_set.txt", sep="\t", dtype=str,
                     usecols=["hgnc_id", "symbol", "name", "location", "ensembl_gene_id", "uniprot_ids", "entrez_id"])
    by_symbol = hg.set_index("symbol")
    by_entrez = hg.dropna(subset=["entrez_id"]).set_index("entrez_id")["symbol"].to_dict()

    def gene(symbol: str) -> str | None:
        if symbol not in by_symbol.index:
            return None
        r = by_symbol.loc[symbol]
        nid = gene_node_id(r["hgnc_id"])
        store.add_node(nid, "Gene", symbol, hgnc_id=r["hgnc_id"], full_name=r["name"], location=r["location"],
                       ensembl=r["ensembl_gene_id"], uniprot=r["uniprot_ids"])
        return nid

    with store.batch():
        store.add_node(DISEASE, "Disease", cfg["disease"]["name"], mondo=DISEASE,
                       open_targets=cfg["disease"]["open_targets_id"])
        for _, (did, dname) in DATASETS.items():
            store.add_node(did, "Dataset", dname, source="LinkedOmics CPTAC-CCRCC discovery cohort")

        # ---------------------------------------------------------- omics (M2)
        thr = cfg["omics"]
        de = pd.read_parquet(PROCESSED / "de_results.parquet").dropna(subset=["hgnc_symbol"])
        de = de.sort_values("padj").drop_duplicates(["layer", "hgnc_symbol"])
        for r in de.itertuples():
            gid = gene(r.hgnc_symbol)
            if gid is None:
                continue
            sig = r.padj < thr["padj"] and abs(r.log2fc) >= thr["min_abs_log2fc"][r.layer]
            pred = ("upregulated_in" if r.log2fc > 0 else "downregulated_in") if sig else "not_changed_in"
            did = DATASETS[r.layer][0]
            ev = store.add_evidence(evidence_id("omics_de", r.layer, r.hgnc_symbol), "omics_de", did, TODAY,
                                    layer=r.layer, gene=r.hgnc_symbol, log2fc=round(r.log2fc, 4),
                                    padj=float(f"{r.padj:.3g}"), pvalue=float(f"{r.pvalue:.3g}"), method=r.method,
                                    n_pairs=int(r.n_pairs), n_tumor=int(r.n_tumor), n_normal=int(r.n_normal),
                                    unit=r.unit, significant=bool(sig))
            store.add_claim(gid, pred, DISEASE, ev)
            store.add_claim(gid, "measured_in", did, ev)

        cnv = pd.read_parquet(PROCESSED / "cnv_frequencies.parquet").dropna(subset=["hgnc_symbol"])
        for r in cnv.drop_duplicates("hgnc_symbol").itertuples():
            gid = gene(r.hgnc_symbol)
            if gid is None:
                continue
            recurrent = max(r.loss_freq, r.gain_freq) >= 0.3
            pred = ("lost_in" if r.loss_freq >= r.gain_freq else "gained_in") if recurrent else "copy_number_assessed_in"
            ev = store.add_evidence(evidence_id("omics_cnv", r.hgnc_symbol), "omics_cnv", DATASETS["cnv"][0], TODAY,
                                    gene=r.hgnc_symbol, loss_freq=round(r.loss_freq, 3), gain_freq=round(r.gain_freq, 3),
                                    n_tumor=int(r.n_tumor), arm=r.arm, threshold_abs_log2ratio=0.2)
            store.add_claim(gid, pred, DISEASE if recurrent else DATASETS["cnv"][0], ev)

        corr = pd.read_parquet(PROCESSED / "rna_protein_correlation.parquet").dropna(subset=["hgnc_symbol"])
        for r in corr.drop_duplicates("hgnc_symbol").itertuples():
            gid = gene(r.hgnc_symbol)
            if gid is None or r.padj >= 0.05:
                continue
            ev = store.add_evidence(evidence_id("omics_rna_protein", r.hgnc_symbol), "omics_rna_protein_correlation",
                                    DATASETS["protein"][0], TODAY, gene=r.hgnc_symbol, spearman_rho=round(r.spearman_rho, 3),
                                    padj=float(f"{r.padj:.3g}"), n=int(r.n))
            store.add_claim(gid, "rna_protein_correlated_in", DATASETS["protein"][0], ev)

        # ---------------------------------------------------------- literature (M3)
        pubs = pd.read_parquet(LITERATURE / "publications.parquet").set_index("pmid")

        def publication(pmid: int) -> str:
            nid = f"PMID:{pmid}"
            p = pubs.loc[pmid]
            store.add_node(nid, "Publication", p["title"][:300] or nid, pmid=int(pmid), year=p["year"],
                           journal=p["journal"][:200], url=f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/")
            return nid

        def entity(text: str) -> str:
            g = resolve_gene(text)
            if g and (gid := gene(g)):
                return gid
            if is_ccrcc(text):
                return DISEASE
            return store.add_node(f"CONCEPT:{_slug(text)}", "Concept", text.strip()[:120])

        rel = pd.read_parquet(LITERATURE / "relations.parquet")
        for r in rel.itertuples():
            pid = publication(int(r.pmid))
            s, o = entity(r.subject), entity(r.object)
            ev = store.add_evidence(evidence_id("literature", r.relation_id), "literature_relation", pid, TODAY,
                                    relation_id=r.relation_id, pmid=int(r.pmid), sentence=r.literal_sentence,
                                    confidence=float(r.confidence), subject_text=r.subject, object_text=r.object,
                                    extractor="extractor-relaciones (Claude Sonnet 5.5), verbatim-validated")
            store.add_claim(s, r.relation, o, ev)
            store.add_claim(pid, "reports", s, ev)

        ann = pd.read_parquet(LITERATURE / "annotations.parquet")
        genes_ann = ann[(ann["type"] == "Gene") & ann["identifier"].notna()].drop_duplicates(["pmid", "identifier"])
        for r in genes_ann.itertuples():
            sym = None
            for ident in str(r.identifier).split(";"):
                sym = by_entrez.get(ident.strip())
                if sym:
                    break
            if not sym or (gid := gene(sym)) is None:
                continue
            pid = publication(int(r.pmid))
            ev = store.add_evidence(evidence_id("pubtator", r.pmid, r.identifier), "pubtator_annotation", "PubTator3",
                                    TODAY, pmid=int(r.pmid), ncbi_gene=str(r.identifier), mention=r.text)
            store.add_claim(pid, "mentions", gid, ev)

        # gene mentions per publication (PubTator-normalized): mention text -> gene node id
        mentions_by_pmid: dict[int, dict[str, str]] = {}
        for r in genes_ann.itertuples():
            sym = next((by_entrez.get(i.strip()) for i in str(r.identifier).split(";") if by_entrez.get(i.strip())), None)
            if sym and sym in by_symbol.index and isinstance(r.text, str) and len(r.text) >= 2:
                mentions_by_pmid.setdefault(int(r.pmid), {})[r.text] = gene_node_id(by_symbol.loc[sym, "hgnc_id"])

        chunks = pd.read_parquet(LITERATURE / "chunks.parquet")
        for r in chunks.itertuples():
            pid = publication(int(r.pmid))
            cid = f"CHUNK:{r.chunk_id}"
            store.add_node(cid, "Chunk", r.text[:120], text=r.text, pmid=int(r.pmid))
            ev = store.add_evidence(evidence_id("chunk", r.chunk_id), "provenance", "PubTator3 abstract", TODAY,
                                    pmid=int(r.pmid))
            store.add_claim(cid, "part_of", pid, ev)
            for mention, gid in mentions_by_pmid.get(int(r.pmid), {}).items():
                if re.search(rf"(?<![A-Za-z0-9]){re.escape(mention)}(?![A-Za-z0-9])", r.text):
                    store.add_claim(cid, "mentions", gid, ev)

        # ---------------------------------------------------------- external (M5 cache)
        ot = pd.read_parquet(PROCESSED / "external_ot_associations.parquet")
        for r in ot.itertuples():
            gid = gene(r.symbol)
            if gid is None:
                continue
            ev = store.add_evidence(evidence_id("open_targets_assoc", r.ensembl), "open_targets_association",
                                    "Open Targets Platform", TODAY, score=round(r.score, 4), ensembl=r.ensembl,
                                    gene=r.symbol, genetic_association=round(r.dt_genetic_association, 4),
                                    somatic_mutation=round(r.dt_somatic_mutation, 4),
                                    literature=round(r.dt_literature, 4), known_drug=round(r.dt_known_drug, 4),
                                    disease=cfg["disease"]["open_targets_id"],
                                    url=f"https://platform.opentargets.org/evidence/{r.ensembl}/{cfg['disease']['open_targets_id']}")
            store.add_claim(gid, "associated_with", DISEASE, ev)

        details = json.loads((PROCESSED / "external_targets.json").read_text())
        for sym, d in details.items():
            gid = gene(sym)
            if gid is None:
                continue
            t = d.get("open_targets")
            if t:
                ev_t = store.add_evidence(evidence_id("open_targets_tract", sym), "open_targets_tractability",
                                          "Open Targets Platform", TODAY, gene=sym,
                                          tractable_modalities=t["tractable_modalities"],
                                          has_approved_or_clinical_drug=t["has_approved_or_clinical_drug"],
                                          labels=t["tractability_labels"][:20])
                store.add_claim(gid, "tractability_assessed_in", DISEASE, ev_t)
                for drug in t["drugs"]:
                    did = f"CHEMBL:{drug['chembl']}"
                    store.add_node(did, "Drug", drug["name"], chembl=drug["chembl"])
                    ev = store.add_evidence(evidence_id("open_targets_drug", drug["chembl"], sym), "open_targets_drug",
                                            "Open Targets Platform", TODAY, gene=sym, drug=drug["name"],
                                            max_clinical_stage=drug["max_stage"],
                                            url=f"https://platform.opentargets.org/drug/{drug['chembl']}")
                    store.add_claim(did, "targets", gid, ev)
                for pw in t["pathways"]:
                    pid = f"REACTOME:{pw['reactome']}"
                    store.add_node(pid, "Pathway", pw["name"], reactome=pw["reactome"], top_level=pw["top_level"])
                    ev = store.add_evidence(evidence_id("reactome", pw["reactome"], sym), "pathway_membership",
                                            "Reactome via Open Targets", TODAY, gene=sym, reactome=pw["reactome"])
                    store.add_claim(gid, "in_pathway", pid, ev)
            for p in d.get("string") or []:
                pg = gene(p["partner"])
                if pg is None:
                    continue
                ev = store.add_evidence(evidence_id("string", *sorted([sym, p["partner"]])), "string_interaction",
                                        "STRING", TODAY, a=sym, b=p["partner"], score=p["score"],
                                        experimental=p["experimental"])
                store.add_claim(gid, "interacts_with", pg, ev)
            u = d.get("uniprot")
            if u:
                prot = f"UNIPROT:{u['accession']}"
                store.add_node(prot, "Protein", u["name"] or u["accession"], accession=u["accession"],
                               function=u["function"], subcellular_location=u["subcellular_location"])
                ev = store.add_evidence(evidence_id("uniprot", u["accession"]), "uniprot_entry", "UniProtKB/Swiss-Prot",
                                        TODAY, accession=u["accession"],
                                        url=f"https://www.uniprot.org/uniprotkb/{u['accession']}")
                store.add_claim(gid, "encodes", prot, ev)

    from .literature_review import load as load_reviewed_literature
    load_reviewed_literature(store)
    stats = store.stats()
    (PROCESSED / "graph_stats.json").write_text(json.dumps(stats, indent=1))
    return stats


if __name__ == "__main__":
    print(json.dumps(build(), indent=1))
