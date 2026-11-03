"""How often does the verifier catch a wrong statement? Deliberately corrupted copies of verified statements.

Every evidence-bearing statement of a run (agent findings, narrative paragraphs, targets) that passes verification is
copied and corrupted in four ways, one at a time:

- number:    the first printed number gets a different last digit (inside or outside a canonical annotation)
- direction: the first direction word is flipped (higher <-> lower, up <-> down, increased <-> decreased, ...)
- gene:      a gene-specific statement is attributed to another gene of the same run
- evidence:  a gene-specific statement cites the evidence of a statement about another gene instead of its own

A corruption counts as caught when the verifier rejects the statement or removes any of its evidence ids, the same rule
`verify_final` applies. Rates are reported per kind and only over statements where the corruption could be applied.
The direction check is narrow by design (significant omics_de evidence only), so its rate is reported twice: over all
flipped statements and over those that cite significant omics_de evidence for the gene concerned.
"""
from __future__ import annotations

import json
import random
import re
from pathlib import Path

from . import semantic
from .verify import Verifier, _statements

FLIP = [("higher", "lower"), ("increased", "decreased"), ("elevated", "reduced"), ("upregulated", "downregulated"),
        ("up-regulated", "down-regulated"), ("overexpressed", "underexpressed"), ("up", "down"), ("gained", "depleted")]
FLIP_MAP = {**{a: b for a, b in FLIP}, **{b: a for a, b in FLIP}}
WORD = re.compile(r"\b(" + "|".join(re.escape(w) for w in sorted(FLIP_MAP, key=len, reverse=True)) + r")\b", re.I)


def _text_key(stmt: dict) -> str:
    return next(k for k in ("statement", "text", "rationale") if stmt.get(k))


def corrupt_number(stmt: dict) -> dict | None:
    key = _text_key(stmt)
    nums = semantic.numbers_in(stmt[key], with_spans=True)
    if not nums:
        return None
    n = nums[0]
    s, e = n["start"], n["end"]
    digits = [i for i in range(s, e) if stmt[key][i].isdigit()]
    i = digits[-1]
    new = str((int(stmt[key][i]) + 1) % 10)
    text = stmt[key][:i] + new + stmt[key][i + 1:]
    return {**stmt, key: text}


def corrupt_direction(stmt: dict) -> dict | None:
    key = _text_key(stmt)
    m = WORD.search(stmt[key])
    if not m:
        return None
    w = m.group(0)
    rep = FLIP_MAP[w.lower()]
    rep = rep.capitalize() if w[0].isupper() else rep
    return {**stmt, key: stmt[key][:m.start()] + rep + stmt[key][m.end():]}


def _de_gene(stmt: dict, store) -> bool:
    """Does the statement cite significant omics_de evidence for its gene, or for the one gene its sentence names?"""
    ev = store.get_evidence([i for i in stmt.get("evidence_ids", []) if isinstance(i, str)])
    genes = {(e.get("gene") or "").upper() for e in ev if e.get("type") == "omics_de" and e.get("significant")}
    if stmt.get("gene"):
        return stmt["gene"].upper() in genes
    return bool(genes)


def run(run_dir: Path, seed: int = 7, out_path: Path | None = None) -> dict:
    findings = json.loads((Path(run_dir) / "findings.json").read_text(encoding="utf-8"))
    from .provenance import load_basis
    v = Verifier(basis=load_basis(Path(run_dir)))
    rng = random.Random(seed)

    def caught(stmt) -> bool:
        ok, _ = v.check_statement(stmt)
        return ok is None or len(ok["evidence_ids"]) != len(stmt.get("evidence_ids", []))

    stmts = [(loc, s) for loc, s in _statements(findings) if not caught(s)]
    gene_stmts = [(loc, s) for loc, s in stmts if s.get("gene")]
    out = {k: {"applied": 0, "caught": 0, "missed_examples": []} for k in ("number", "direction", "gene", "evidence")}
    out["direction"].update({"applied_with_de_evidence": 0, "caught_with_de_evidence": 0})

    def record(kind, loc, mutant, extra=None):
        hit = caught(mutant)
        o = out[kind]
        o["applied"] += 1
        o["caught"] += hit
        if extra is not None:
            o["applied_with_de_evidence"] += extra
            o["caught_with_de_evidence"] += extra and hit
        if not hit and len(o["missed_examples"]) < 5:
            o["missed_examples"].append({"location": loc, "text": mutant[_text_key(mutant)][:240]})

    for loc, s in stmts:
        m = corrupt_number(s)
        if m:
            record("number", loc, m)
        m = corrupt_direction(s)
        if m:
            record("direction", loc, m, extra=_de_gene(s, v.store))
    for loc, s in gene_stmts:
        others = sorted({o["gene"] for _, o in gene_stmts if o["gene"].upper() != s["gene"].upper()})
        if others:
            record("gene", loc, {**s, "gene": rng.choice(others)})
        donors = [o for _, o in gene_stmts if o["gene"].upper() != s["gene"].upper()]
        if donors:
            d = rng.choice(donors)
            record("evidence", loc, {**s, "evidence_ids": d["evidence_ids"]})
    for o in out.values():
        o["rate"] = round(o["caught"] / o["applied"], 3) if o["applied"] else None
    d = out["direction"]
    d["rate_with_de_evidence"] = (round(d["caught_with_de_evidence"] / d["applied_with_de_evidence"], 3)
                                  if d["applied_with_de_evidence"] else None)
    res = {"run": Path(run_dir).name, "statements_passing": len(stmts), "seed": seed, "corruptions": out,
           "scope": "Detection of deliberately corrupted copies of statements that pass verification; a missed "
                    "corruption is a statement the automatic checks would accept although it is wrong."}
    from .paths import ROOT
    (out_path or ROOT / "evals" / "verifier_sabotage.json").write_text(json.dumps(res, indent=1), encoding="utf-8")
    return res
