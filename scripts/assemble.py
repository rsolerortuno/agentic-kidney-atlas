"""Assemble runs/<date>/findings.json from convergence.json, verified agent outputs and narrative.json."""
import json
import sys

from ccrcc.assemble import assemble
from ccrcc.integrate import run_dir

if __name__ == "__main__":
    print(json.dumps(assemble(run_dir(sys.argv[1] if len(sys.argv) > 1 else None)), indent=1))
