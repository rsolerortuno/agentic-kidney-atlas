"""Deterministic convergence scoring -> runs/<date>/convergence.json (weights in config.yaml)."""
import json
import sys

from pathlib import Path

from ccrcc.graph import GraphStore
from ccrcc.integrate import write
from ccrcc.provenance import convergence_ids, start_basis

if __name__ == "__main__":
    # This script is the first step of the agent run (/informe). It also starts runs/<run>/basis.json, the record of
    # the evidence values the agents will see (stale-narrative guard, src/ccrcc/provenance.py).
    day = sys.argv[1] if len(sys.argv) > 1 else None
    # A fresh agent run must not inherit anything from an earlier one: old agent files, narrative or findings in the same
    # run folder could be read by the agents or merged into the new report (this made replicate runs less independent).
    from ccrcc.integrate import run_dir
    for pat in ("agent_*.json", "narrative.json", "findings.json", "verification_*.json", "basis.json"):
        for f in run_dir(day).glob(pat):
            f.unlink()
    out = write(day=day)
    conv = json.loads(Path(out["path"]).read_text())
    basis = start_basis(Path(out["path"]).parent, GraphStore(), convergence_ids(conv))
    out["basis_evidence"] = len(basis["evidence"])
    print(json.dumps(out, indent=1))
