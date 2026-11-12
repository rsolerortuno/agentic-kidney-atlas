"""external-mcp: Open Targets, STRING and UniProt (disk-cached), with graph evidence ids when available."""
from __future__ import annotations

import pandas as pd
from mcp.server.fastmcp import FastMCP

from .. import external as ext
from ..config import load
from ..graph import evidence_id
from ..paths import AGENT_CACHE, PROCESSED, RAW
from .common import canonical_gene, clean, existing

mcp = FastMCP("external-mcp")
ext.WRITE_ROOT = AGENT_CACHE   # agent lookups never write into the frozen input tree (data/cache)


def _ensembl(symbol: str) -> str | None:
    hg = pd.read_csv(RAW / "hgnc_complete_set.txt", sep="\t", usecols=["symbol", "ensembl_gene_id"], dtype=str)
    m = hg[hg["symbol"] == symbol]["ensembl_gene_id"].dropna()
    return m.iloc[0] if len(m) else None


@mcp.tool()
def open_targets_association(gene: str) -> dict:
    """Open Targets Platform overall association score between a gene and clear cell renal carcinoma
    (MONDO_0005005), its rank among associated targets, and the evidence_id."""
    g = canonical_gene(gene)
    ot = pd.read_parquet(PROCESSED / "external_ot_associations.parquet").reset_index(drop=True)
    sub = ot[ot["symbol"] == g]
    if sub.empty:
        return {"gene": g, "message": f"not among the top {len(ot)} Open Targets associations for ccRCC"}
    r = sub.iloc[0]
    return clean({"gene": g, "score": round(r["score"], 4), "rank": int(sub.index[0]) + 1, "of": len(ot),
                  "evidence_id": existing(evidence_id("open_targets_assoc", r["ensembl"]))})


@mcp.tool()
def open_targets_target(gene: str) -> dict:
    """Open Targets target profile: tractability (small molecule, antibody, PROTAC), known drugs with
    maximum clinical stage, and Reactome pathways. Includes evidence ids for each drug/tractability record."""
    g = canonical_gene(gene)
    ens = _ensembl(g)
    if not ens:
        return {"gene": g, "message": "no Ensembl id"}
    t = ext.ot_target(ens)
    if not t:
        return {"gene": g, "message": "not in Open Targets"}
    t["tractability_evidence_id"] = existing(evidence_id("open_targets_tract", g))
    for d in t["drugs"]:
        d["evidence_id"] = existing(evidence_id("open_targets_drug", d["chembl"], g))
    t["pathways"] = t["pathways"][:15]
    return clean(t)


@mcp.tool()
def string_interactions(gene: str, limit: int = 10) -> list[dict]:
    """High-confidence (score >= 0.7) STRING protein interaction partners for a gene."""
    g = canonical_gene(gene)
    out = ext.string_partners(g, limit=limit)
    for p in out:
        p["evidence_id"] = existing(evidence_id("string", *sorted([g, p["partner"]])))
    return clean(out)


@mcp.tool()
def uniprot_entry(gene: str) -> dict:
    """UniProtKB/Swiss-Prot entry for the human gene: accession, protein name, function, subcellular location."""
    g = canonical_gene(gene)
    u = ext.uniprot(g)
    if not u:
        return {"gene": g, "message": "no reviewed UniProt entry"}
    u["evidence_id"] = existing(evidence_id("uniprot", u["accession"]))
    return clean(u)


@mcp.tool()
def disease_top_targets(n: int = 25) -> list[dict]:
    """Top Open Targets associated targets for ccRCC (MONDO_0005005), with scores."""
    ot = pd.read_parquet(PROCESSED / "external_ot_associations.parquet").head(min(n, 200))
    return clean([{"gene": r.symbol, "score": round(r.score, 4),
                   "evidence_id": existing(evidence_id("open_targets_assoc", r.ensembl))} for r in ot.itertuples()])


if __name__ == "__main__":
    load()
    mcp.run()
