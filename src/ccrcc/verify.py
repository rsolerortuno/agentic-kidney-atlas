"""Verification of agent outputs and final findings (deterministic, no LLM).

A statement is supported only if it cites >= 1 evidence id that (a) exists in the graph and (b) when the
statement is about a gene, is attached to a claim involving that gene, and (c) its numbers and its direction of
change agree with the evidence it cites (semantic.py). Invalid ids are removed; statements left without support,
or contradicted by their own evidence, are dropped and logged.
"""
from __future__ import annotations

import json
from pathlib import Path

from . import semantic
from .graph import GraphStore
from .literature_entities import resolve_gene

AGENT_FILES = ("agent_omics.json", "agent_literature.json", "agent_external.json", "agent_cells.json")


def _ranking_values(conv_path: Path | None) -> dict[str, set[float]]:
    """Numeric values of each gene's convergence row (rank, score, components, fold changes) that statements
    may quote in addition to their evidence."""
    if not conv_path or not Path(conv_path).exists():
        return {}
    conv = json.loads(Path(conv_path).read_text(encoding="utf-8"))
    out = {}
    for r in conv.get("ranking", []):
        row = {k: v for k, v in r.items() if k != "evidence_ids"}
        out[r["gene"]] = semantic._numeric_leaves(row, set())
    for name in ("profiles",):
        for prof, rows in (conv.get(name) or {}).items():
            for r in rows:
                out.setdefault(r["gene"], set()).update(semantic._numeric_leaves(r, set()))
    return out


class Verifier:
    def __init__(self, store: GraphStore | None = None, conv_path: Path | None = None, semantic_check: bool = True):
        self.store = store or GraphStore()
        self._gene_ev: dict[str, set[str]] = {}
        self.semantic_check = semantic_check
        self.ranking = _ranking_values(conv_path)

    def _gene_evidence(self, gene: str) -> set[str] | None:
        sym = resolve_gene(gene) or gene
        node = self.store.resolve(sym)
        if not node or node["type"] != "Gene":
            return None
        if node["id"] not in self._gene_ev:
            rows = self.store.conn.execute(
                "SELECT evidence_id FROM claims WHERE subject = ? OR object = ?", (node["id"], node["id"])).fetchall()
            self._gene_ev[node["id"]] = {r[0] for r in rows}
        return self._gene_ev[node["id"]]

    def check_statement(self, stmt: dict) -> tuple[dict | None, dict]:
        ids = [i for i in stmt.get("evidence_ids", []) if isinstance(i, str)]
        existing = set(self.store.check_evidence_ids(ids)["valid"]) if ids else set()
        gene = stmt.get("gene")
        linked = existing
        if gene:
            gene_ev = self._gene_evidence(gene)
            linked = existing & gene_ev if gene_ev is not None else set()
        valid = [i for i in ids if i in linked]
        text = stmt.get("statement") or stmt.get("text") or stmt.get("rationale") or ""
        log = {"text": text[:160], "gene": gene,
               "cited": len(ids), "nonexistent": sorted(set(ids) - existing),
               "not_linked_to_gene": sorted(existing - linked), "kept": bool(valid)}
        if valid and self.semantic_check:
            ev = self.store.get_evidence(valid)
            extra: set[float] = set()
            for g in {gene, *[e.get("gene") for e in ev]} - {None}:
                extra |= self.ranking.get(resolve_gene(g) or g, set())
            sem = semantic.check(text, resolve_gene(gene) or gene if gene else None, ev, extra)
            log.update(sem)
            if sem["unmatched_numbers"] or sem["direction_contradiction"]:
                log["kept"] = False
                log["semantic_failure"] = True
                return None, log
        if not valid:
            return None, log
        return {**stmt, "evidence_ids": valid}, log


def verify_agents(run_dir: Path) -> dict:
    v = Verifier(conv_path=run_dir / "convergence.json")
    report = {"agents": {}}
    for name in AGENT_FILES:
        path = run_dir / name
        if not path.exists():
            report["agents"][name] = {"missing": True}
            continue
        data = json.loads(path.read_text(encoding="utf-8"))
        kept, logs = [], []
        for f in data.get("findings", []):
            ok, log = v.check_statement(f)
            logs.append(log)
            if ok:
                kept.append(ok)
        data["findings"] = kept
        (run_dir / name.replace(".json", ".verified.json")).write_text(json.dumps(data, indent=1), encoding="utf-8")
        report["agents"][name] = {"submitted": len(logs), "kept": len(kept), "dropped": len(logs) - len(kept),
                                  "dropped_detail": [l for l in logs if not l["kept"]],
                                  "semantic_failures": sum(1 for l in logs if l.get("semantic_failure")),
                                  "ids_removed": sum(len(l["nonexistent"]) + len(l["not_linked_to_gene"]) for l in logs)}
    (run_dir / "verification_agents.json").write_text(json.dumps(report, indent=1), encoding="utf-8")
    return report


def _statements(findings: dict):
    """Yield (location, statement) for every evidence-bearing statement in findings.json."""
    for agent, block in findings.get("agent_findings", {}).items():
        for i, f in enumerate(block.get("findings", [])):
            yield f"agent_findings.{agent}[{i}]", f
    narr = findings.get("narrative", {})
    for i, p in enumerate(narr.get("summary", [])):
        yield f"narrative.summary[{i}]", p
    for sec, paras in narr.get("sections", {}).items():
        for i, p in enumerate(paras):
            yield f"narrative.sections.{sec}[{i}]", p
    for t in findings.get("targets", []):
        yield f"targets.{t.get('gene')}", t


def verify_final(path: Path) -> dict:
    """Final check on findings.json: count statements without valid, gene-linked evidence."""
    v = Verifier(conv_path=path.parent / "convergence.json")
    findings = json.loads(path.read_text(encoding="utf-8"))
    unsupported, total = [], 0
    for loc, stmt in _statements(findings):
        total += 1
        ok, log = v.check_statement(stmt)
        if not ok or len(ok["evidence_ids"]) != len(stmt.get("evidence_ids", [])):
            unsupported.append({"location": loc, **log})
    report = {"file": str(path), "statements": total, "unsupported": len(unsupported),
              "semantic_failures": sum(1 for u in unsupported if u.get("semantic_failure")), "detail": unsupported}
    (path.parent / "verification_final.json").write_text(json.dumps(report, indent=1), encoding="utf-8")
    return report
