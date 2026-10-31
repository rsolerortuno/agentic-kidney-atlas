"""Deterministic convergence scoring -> runs/<date>/convergence.json (weights in config.yaml)."""
import json
import sys

from ccrcc.integrate import write

if __name__ == "__main__":
    day = sys.argv[1] if len(sys.argv) > 1 else None
    print(json.dumps(write(day=day), indent=1))
