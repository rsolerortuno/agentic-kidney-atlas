"""Assemble runs/<date>/findings.json from convergence.json, verified agent outputs and narrative.json."""
import json
import sys

from ccrcc.assemble import assemble
from ccrcc.integrate import current_run, run_dir

if __name__ == "__main__":
    print(json.dumps(assemble(run_dir(sys.argv[1]) if len(sys.argv) > 1 else current_run()), indent=1))
