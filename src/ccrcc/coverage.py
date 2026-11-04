"""Integration coverage: which verified specialist results did the narrative leave out?

The orchestrator writes the narrative from the agents' verified findings. The pipeline benchmark showed that the
specialists found decisive results (a protein that does not replicate, a second trial, a drug screen) that never reached
the report. This check lists, before assembly:

- key findings (marked ``"key": true`` by the agent that made them) none of whose evidence ids the narrative cites;
- evidence types that some verified finding cites but the narrative never does.

It is advisory: the orchestrator either integrates each item or says in ``limitations`` why it was left out.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

EV_TYPE = re.compile(r"^EV:([a-z0-9_]+):")


def _narrative_ids(n: dict) -> set[str]:
    paras = list(n.get("summary", [])) + list(n.get("targets", []))
    for v in (n.get("sections") or {}).values():
        paras += list(v)
    return {i for p in paras for i in p.get("evidence_ids", [])}


def _types(ids) -> set[str]:
    return {m.group(1) for i in ids if (m := EV_TYPE.match(i))}


def coverage(run: Path) -> dict:
    n = json.loads((run / "narrative.json").read_text(encoding="utf-8"))
    cited = _narrative_ids(n)
    limitations = " ".join(x.get("text", "") for x in n.get("limitations", []))
    missing_key, agent_ids, n_key = [], set(), 0
    for p in sorted(run.glob("agent_*.verified.json")):
        for f in json.loads(p.read_text(encoding="utf-8")).get("findings", []):
            ids = set(f.get("evidence_ids", []))
            agent_ids |= ids
            if f.get("key"):
                n_key += 1
                if not ids & cited:
                    missing_key.append({"agent": p.name.replace(".verified.json", ""), "gene": f.get("gene"),
                                        "statement": f.get("statement", "")[:400],
                                        "explained_in_limitations": bool(f.get("gene")) and f["gene"] in limitations})
    out = {"key_findings": n_key, "key_findings_not_cited": missing_key,
           "evidence_types_not_cited": sorted(_types(agent_ids) - _types(cited)),
           "share_of_finding_evidence_cited": round(len(agent_ids & cited) / len(agent_ids), 3) if agent_ids else None}
    (run / "coverage.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    return out
