"""Helpers shared by the MCP servers."""
from __future__ import annotations

import math
from functools import lru_cache

from ..graph import GraphStore
from ..literature_entities import resolve_gene


@lru_cache(maxsize=1)
def store() -> GraphStore:
    return GraphStore()


def canonical_gene(symbol: str) -> str:
    """Accept approved, previous or dataset-era symbols (e.g. NDUFA4L2 -> COXFA4L2)."""
    return resolve_gene(symbol) or symbol.strip().upper()


def existing(ev_id: str) -> str | None:
    """Return the evidence id only if it exists in the graph (agents may only cite existing ids)."""
    return ev_id if store().check_evidence_ids([ev_id])["valid"] else None


def clean(obj):
    """Make values JSON-safe (NaN -> None, numpy scalars -> Python)."""
    if isinstance(obj, dict):
        return {k: clean(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [clean(v) for v in obj]
    if hasattr(obj, "item"):
        obj = obj.item()
    if isinstance(obj, float) and math.isnan(obj):
        return None
    return obj


def applied_record(gene: str, ev_type: str, missing: str) -> dict:
    """The evidence record of one gene for a result loaded by kg_applied.py (treatment, subtype, ADC, immune,
    cell replication, confirmatory cohort). The record itself is returned with its evidence_id: cite that id."""
    g = canonical_gene(gene)
    from ..graph import evidence_id
    eid = evidence_id(ev_type, g)
    rec = store().get_evidence([eid])
    if not rec:
        return {"gene": g, "message": missing}
    out = {k: v for k, v in rec[0].items() if k not in ("source", "retrieved", "date", "type")}
    return clean({**out, "evidence_id": eid})
