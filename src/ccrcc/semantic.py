"""Semantic verification of statements against the evidence they cite (deterministic, no LLM).

Existence and gene linkage (verify.py) prove that a statement points at real evidence about the right gene. They do
not prove that the statement says what the evidence says. This module adds two checks:

1. Numbers. Every number written in a statement must match a numeric value of the evidence it cites, of the
   gene's row in the convergence ranking, or a count from the project's own summary files (cohort sizes, cell counts),
   within the rounding the statement uses (6.3 matches 6.312; 97.6 % matches 0.976; 1.7e-49 matches 1.72e-49).
   Small integers (0-10), years and numbers after "phase", "rank", "top" or "#" are not checked.
2. Direction. If a gene statement uses only "up" words (or only "down" words) and every significant omics_de
   evidence it cites for that gene points the other way, the statement contradicts its evidence.
"""
from __future__ import annotations

import json
import math
import re
from functools import lru_cache

from .paths import PROCESSED

NUM = re.compile(r"(?<![\w.])([-−])?(\d+(?:\.\d+)?(?:[eE][-−+]?\d+)?)(\s?%)?(?!\w)")
# Words after which a number is a label, not a measurement. Whole words only: single letters ("p", "q") would exempt
# p- and q-values and any word ending in that letter ("from 0.5"). Numbers glued to letters ("M13", "CD70", "3p") are
# never extracted by NUM in the first place.
EXEMPT_BEFORE = re.compile(r"(\b(phase|rank|ranked|top|grade|stage|tier|version|step|chromosome|chr|pmid:?|uniprot|chembl)"
                           r"|#|[A-Za-z]-)\s*$", re.I)
THRESHOLD_BEFORE = re.compile(r"(\bat|<|≤|<=|>|≥|>=|threshold|cut-?off|alpha|α|fdr|bh)\s*$", re.I)
THRESHOLDS = {0.001, 0.01, 0.05, 0.1, 0.25}   # significance levels written as such ("not significant at 0.05")
DRUG_AFTER = re.compile(r"^\s+[A-Z][A-Z0-9-]{3,}")   # "IODINE I 124 GIRENTUXIMAB": part of a drug name
UP = re.compile(r"\b(up|upregulated|up-regulated|increased|elevated|higher|overexpress\w*|induced|gained)\b", re.I)
DOWN = re.compile(r"\b(down|downregulated|down-regulated|decreased|reduced|lower|depleted|suppressed)\b", re.I)


def _numeric_leaves(obj, out: set[float]) -> set[float]:
    if isinstance(obj, bool):
        return out
    if isinstance(obj, (int, float)):
        if math.isfinite(obj):
            out.add(float(obj))
    elif isinstance(obj, dict):
        for v in obj.values():
            _numeric_leaves(v, out)
    elif isinstance(obj, (list, tuple)):
        for v in obj:
            _numeric_leaves(v, out)
    elif isinstance(obj, str):
        # numbers inside text values (a literature sentence quoting "P = 0.001", a drug called "XL-820")
        for m in re.finditer(r"\d+(?:\.\d+)?(?:[eE][-+]?\d+)?", obj.replace("\u2212", "-")):
            try:
                out.add(float(m.group(0)))
            except ValueError:
                pass
    return out


@lru_cache(maxsize=1)
def project_facts() -> frozenset:
    """Numbers stated in the project's own summary files (cohort sizes, counts) that statements may quote."""
    vals: set[float] = set()
    for name in ("omics_summary.json", "sc_meta.json", "layers_summary.json", "graph_stats.json",
                 "spatial_summary.json", "layer_mutations_meta.json", "layer_dependency_lines.json",
                 "discovery_summary.json", "structure_summary.json", "safety_clinical_summary.json"):
        p = PROCESSED / name
        if p.exists():
            try:
                _numeric_leaves(json.loads(p.read_text(encoding="utf-8")), vals)
            except json.JSONDecodeError:
                pass
    # counts only (cohort sizes, cells, pairs, lines): a pool of every decimal in the summaries would match almost any
    # number written in a sentence and make the check toothless; decimals must come from the cited evidence.
    return frozenset(v for v in vals if float(v).is_integer() and abs(v) > 10)


def numbers_in(text: str) -> list[dict]:
    out = []
    for m in NUM.finditer(text.replace("−", "-")):
        raw = m.group(2)
        before = text[max(0, m.start() - 12):m.start()]
        val = float(raw)
        decimals = len(raw.split(".")[1].split("e")[0].split("E")[0]) if "." in raw else 0
        sci = "e" in raw.lower()
        is_int = "." not in raw and not sci
        if is_int and (val <= 10 or 1900 <= val <= 2100):
            continue
        if EXEMPT_BEFORE.search(before) or DRUG_AFTER.match(text[m.end():m.end() + 20]):
            continue
        if val in THRESHOLDS and THRESHOLD_BEFORE.search(before):
            continue
        out.append({"text": m.group(0).strip(), "value": -val if m.group(1) else val, "decimals": decimals,
                    "sci": sci, "percent": bool(m.group(3))})
    return out


def _matches(n: dict, candidates) -> bool:
    v = n["value"]
    for c in candidates:
        for target in (c, abs(c), c * 100, abs(c) * 100) if not n["sci"] else (c, abs(c)):
            if n["sci"]:
                if target != 0 and abs(abs(v) - abs(target)) <= 0.051 * abs(target):
                    return True
                continue
            tol = 0.5 * 10 ** (-n["decimals"]) + 1e-9
            if abs(abs(v) - abs(target)) <= tol:
                return True
            # integers written for large values ("178,817 cells" is split by the regex; "0.98" vs 0.976)
    return False


def check(text: str, gene: str | None, evidence: list[dict], extra_values: set[float] | None = None) -> dict:
    """Return unmatched numbers and direction contradictions for one statement."""
    cand: set[float] = set()
    for e in evidence:
        _numeric_leaves(e, cand)
    cand |= set(extra_values or ())
    cand |= set(project_facts())
    unmatched = [n["text"] for n in numbers_in(text) if not _matches(n, cand)]
    contradiction = None
    if gene:
        de = [e for e in evidence if e.get("type") == "omics_de" and e.get("significant")
              and (e.get("gene") or "").upper() == gene.upper()]
        signs = {1 if e["log2fc"] > 0 else -1 for e in de if e.get("log2fc") is not None}
        up, down = bool(UP.search(text)), bool(DOWN.search(text))
        if len(signs) == 1 and up != down:
            said = 1 if up else -1
            if said != next(iter(signs)):
                contradiction = f"text says {'up' if up else 'down'}, cited omics evidence says " \
                                f"{'up' if next(iter(signs)) > 0 else 'down'}"
    return {"unmatched_numbers": unmatched, "direction_contradiction": contradiction}
