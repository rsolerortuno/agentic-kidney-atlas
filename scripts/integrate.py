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
    out = write(day=day)
    conv = json.loads(Path(out["path"]).read_text())
    basis = start_basis(Path(out["path"]).parent, GraphStore(), convergence_ids(conv))
    out["basis_evidence"] = len(basis["evidence"])
    print(json.dumps(out, indent=1))
