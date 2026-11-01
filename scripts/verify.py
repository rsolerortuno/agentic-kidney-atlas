"""Verify evidence ids.

  uv run python scripts/verify.py agents [YYYY-MM-DD]   -> drops agent findings without valid, gene-linked evidence
  uv run python scripts/verify.py final  [YYYY-MM-DD]   -> checks findings.json; exit 1 if any statement is unsupported
"""
import json
import sys

from ccrcc.integrate import run_dir
from ccrcc.verify import verify_agents, verify_final

if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "final"
    d = run_dir(sys.argv[2] if len(sys.argv) > 2 else None)
    if mode == "agents":
        rep = verify_agents(d)
        for name, r in rep["agents"].items():
            print(f"{name}: {r}" if r.get("missing") else
                  f"{name}: submitted={r['submitted']} kept={r['kept']} dropped={r['dropped']} ids_removed={r['ids_removed']} "
                  f"semantic_failures={r.get('semantic_failures', 0)}")
    else:
        rep = verify_final(d / "findings.json")
        print(f"statements checked: {rep['statements']}")
        print(f"unsupported claims: {rep['unsupported']} (of which semantic: numbers or direction not matching the "
              f"cited evidence: {rep.get('semantic_failures', 0)}; stale: cited evidence changed since the agent run: "
              f"{rep.get('stale', 0)})")
        if rep.get("basis") is None:
            print("note: no basis.json for this run (agent step predates M14): staleness not checked")
        if rep.get("stale"):
            print("the graph changed after the agents wrote this run: re-run /informe (scripts/run_informe.*)")
        if rep["unsupported"]:
            print(json.dumps(rep["detail"][:10], indent=1))
            sys.exit(1)
