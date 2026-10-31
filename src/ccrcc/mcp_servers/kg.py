"""kg-mcp: query the knowledge graph (GraphStore) and the literature chunk index."""
from __future__ import annotations

import json

from mcp.server.fastmcp import FastMCP

from ..graph import evidence_id
from ..paths import RUNS
from ..search import search_chunks
from .common import canonical_gene, clean, store

mcp = FastMCP("kg-mcp")


def _resolve(entity: str) -> dict | None:
    s = store()
    return s.resolve(entity) or s.resolve(canonical_gene(entity))


@mcp.tool()
def neighbors(entity: str, predicate: str | None = None, limit: int = 50) -> dict:
    """Neighbors of a graph entity (gene symbol, node id such as 'MONDO:0005005', 'PMID:123', drug name).
    Each neighbor lists the predicate, direction and the evidence_ids supporting the link.
    Optional predicate filter, e.g. 'targets', 'interacts_with', 'in_pathway', 'upregulated_in'."""
    n = _resolve(entity)
    if not n:
        return {"error": f"entity '{entity}' not found"}
    nb = store().neighbors(n["id"], limit=500)
    if predicate:
        nb = [x for x in nb if x["predicate"] == predicate]
    for x in nb:
        x["evidence_ids"] = x["evidence_ids"][:5]
    return clean({"entity": {k: n[k] for k in ("id", "type", "name")}, "neighbors": nb[:limit],
                  "total": len(nb)})


@mcp.tool()
def evidence_for(entity: str, evidence_type: str | None = None, limit: int = 50) -> dict:
    """All evidence records attached to an entity's relations (omics statistics, literature sentences with
    PMID, Open Targets scores, drugs, STRING, UniProt). Optional evidence_type filter, e.g.
    'literature_relation', 'omics_de', 'open_targets_drug', 'open_targets_association'."""
    n = _resolve(entity)
    if not n:
        return {"error": f"entity '{entity}' not found"}
    ev = store().evidence_for(n["id"], limit=2000)
    ev = [e for e in ev if e.get("type") not in ("provenance", "pubtator_annotation")]
    if evidence_type:
        ev = [e for e in ev if e.get("type") == evidence_type]
    seen, out = set(), []
    for e in ev:
        if e["evidence_id"] not in seen:
            seen.add(e["evidence_id"])
            out.append(e)
    return clean({"entity": {k: n[k] for k in ("id", "type", "name")}, "evidence": out[:limit], "total": len(out)})


@mcp.tool()
def search_literature(query: str, k: int = 10) -> list[dict]:
    """BM25 search over sentences of the ingested ccRCC abstracts. Returns sentence text, PMID and the
    evidence_id of the sentence chunk (cite it, or cite literature_relation evidence from evidence_for)."""
    hits = search_chunks(query, k=min(k, 30))
    for h in hits:
        ev = evidence_id("chunk", h["chunk_id"])
        h["evidence_id"] = ev if store().check_evidence_ids([ev])["valid"] else None
        h["url"] = f"https://pubmed.ncbi.nlm.nih.gov/{h['pmid']}/"
    return hits


@mcp.tool()
def check_evidence_ids(ids: list[str]) -> dict:
    """Check which evidence ids exist in the graph. Use before returning findings."""
    return store().check_evidence_ids(ids)


@mcp.tool()
def get_evidence(ids: list[str]) -> list[dict]:
    """Full records for a list of evidence ids (max 50)."""
    return clean(store().get_evidence(ids[:50]))


@mcp.tool()
def convergence_ranking(n: int = 25) -> list[dict]:
    """Deterministic convergence ranking from the latest scripts/integrate.py run (if any)."""
    runs = sorted(p for p in RUNS.glob("*/convergence.json"))
    if not runs:
        return [{"message": "no convergence run yet; the orchestrator runs scripts/integrate.py"}]
    return json.loads(runs[-1].read_text())["ranking"][:n]


if __name__ == "__main__":
    mcp.run()
