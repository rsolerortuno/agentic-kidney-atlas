"""kg-mcp: query the knowledge graph (GraphStore) and the literature chunk index."""
from __future__ import annotations

import json
import os

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
def numeric_annotation(evidence_id: str, field: list, display: str, unit: str = "native") -> dict:
    """Build the exact annotation a number needs to pass verification (docs/numeric-claims.md). field = path of keys
    into the cited record, e.g. ["log2fc"] or ["arms", "sunitinib", "log_hr_per_sd"]; display = the number as you will
    print it, e.g. "2.27" (sign and decimals must match the record); unit = "native", "percent" (registered fractions)
    or "native_percent" (fields already 0-100). Returns `annotation` (paste it once, unchanged, into the statement where
    the number goes) and `numeric_claim` (put it in the finding's "numeric_claims" list), or `error` explaining why the
    number does not match the record. Use it for every number you write."""
    from ..semantic import numeric_claim_text
    recs = store().get_evidence([evidence_id])
    if not recs:
        return {"error": f"evidence id {evidence_id} does not exist"}
    claim = {"evidence_id": evidence_id, "field": field, "display": display, "unit": unit}
    try:
        return {"annotation": numeric_claim_text(claim, recs[0]), "numeric_claim": claim}
    except (ValueError, TypeError) as exc:
        return {"error": str(exc), "record": clean(recs[0])}


def _latest_convergence() -> dict | None:
    from ..integrate import run_dir
    if os.environ.get("CCRCC_RUN_DIR") and (run_dir() / "convergence.json").exists():
        return json.loads((run_dir() / "convergence.json").read_text())
    runs = sorted(p for p in RUNS.glob("[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]/convergence.json"))
    return json.loads(runs[-1].read_text()) if runs else None


@mcp.tool()
def convergence_ranking(n: int = 25, profile: str = "convergence") -> list[dict]:
    """Deterministic ranking from the latest scripts/integrate.py run.

    profile: 'convergence' (all evidence streams, including literature and drugs), 'target' (tumour-up expression,
    purity-robust tumour-intrinsic signal, selective DepMap dependency, driver status, kinase activity, predicted
    tractability) or 'biomarker' (change in either direction, TCGA replication, prognosis, grade, RNA-protein
    agreement). Profiles leave out literature and drug history. Each row carries its components, the M13 values
    (mutation, DepMap, TCGA, purity, kinase) and `rank_interval_90`: 5th-95th percentile rank and probability of
    staying in the top 10/25/100 when the weights are randomly perturbed (robustness)."""
    conv = _latest_convergence()
    if not conv:
        return [{"message": "no convergence run yet; the orchestrator runs scripts/integrate.py"}]
    if profile == "convergence":
        return conv["ranking"][:n]
    rows = (conv.get("profiles") or {}).get(profile)
    if rows is None:
        return [{"message": f"unknown profile '{profile}'; use convergence, target or biomarker"}]
    return rows[:n]


@mcp.tool()
def gene_scores(gene: str) -> dict:
    """Where one gene sits in the convergence ranking and in the target / biomarker profiles (rank, score,
    components, robustness interval). Returns only what the latest run computed; absent = outside the stored top."""
    conv = _latest_convergence()
    if not conv:
        return {"message": "no convergence run yet"}
    g = canonical_gene(gene)
    out = {"gene": g}
    for name, rows in [("convergence", conv["ranking"])] + list((conv.get("profiles") or {}).items()):
        hit = next((r for r in rows if r["gene"] == g), None)
        out[name] = hit if hit else "not in the stored top of this ranking"
    return clean(out)


if __name__ == "__main__":
    mcp.run()
