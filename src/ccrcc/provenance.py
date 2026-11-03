"""Stale-narrative guard (M14).

The agent step writes text against the graph as it is at that moment. If the deterministic pipeline is re-run later
(new data, a code change), evidence ids stay the same but their values can change, and a saved narrative would then
quote numbers the graph no longer contains. That happened once: the M12 report carried agent text written during M11.

`runs/<run>/basis.json` records, for every evidence id the agents and the orchestrator could cite, a fingerprint of the
record's values (everything except the build date). It is written when the agent step starts (`scripts/integrate.py`:
the convergence ids) and completed after the agents finish (`verify.py agents`: the ids they cite). Ids already in the
basis are never overwritten, so a later `run-all` cannot refresh it. `verify.py final` then fails every statement that
cites an id whose current fingerprint differs from the basis, or that is missing from the basis (fail closed): the
narrative must be regenerated with `/informe`.
"""
from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

from .paths import CONFIG

IGNORED = {"evidence_id", "date"}


def fingerprint(record: dict) -> str:
    body = {k: v for k, v in record.items() if k not in IGNORED}
    return hashlib.sha256(json.dumps(body, sort_keys=True, default=str).encode()).hexdigest()[:16]


def fingerprints(store, ids) -> dict[str, str]:
    ids = sorted(set(i for i in ids if isinstance(i, str)))
    out = {}
    for k in range(0, len(ids), 500):
        for e in store.get_evidence(ids[k:k + 500]):
            out[e["evidence_id"]] = fingerprint(e)
    return out


def basis_path(run_dir: Path) -> Path:
    return Path(run_dir) / "basis.json"


def load_basis(run_dir: Path) -> dict | None:
    p = basis_path(run_dir)
    return json.loads(p.read_text(encoding="utf-8")) if p.exists() else None


def start_basis(run_dir: Path, store, ids) -> dict:
    """New agent run: replace any previous basis with the current fingerprints of `ids`."""
    basis = {"started_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
             "config_sha256": hashlib.sha256(CONFIG.read_bytes()).hexdigest(),
             "evidence": fingerprints(store, ids)}
    basis_path(run_dir).write_text(json.dumps(basis, indent=0), encoding="utf-8")
    return basis


def extend_basis(run_dir: Path, store, ids) -> dict | None:
    """Add ids that are not in the basis yet (existing entries keep their agent-time fingerprint)."""
    basis = load_basis(run_dir)
    if basis is None:
        return None
    new = [i for i in set(ids) if i not in basis["evidence"]]
    if new:
        basis["evidence"].update(fingerprints(store, new))
        basis_path(run_dir).write_text(json.dumps(basis, indent=0), encoding="utf-8")
    return basis


def convergence_ids(conv: dict) -> list[str]:
    ids = []
    for rows in [conv.get("ranking", [])] + list((conv.get("profiles") or {}).values()):
        for r in rows:
            for v in r.get("evidence_ids", {}).values():
                ids += v
    return ids
