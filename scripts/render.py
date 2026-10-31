"""Render reports/ccrcc_report.html from the latest runs/<date>/findings.json (or a given date).

With CCRCC_RUN_DIR set (replicate runs, see src/ccrcc/baselines.py) the report is written next to that run's
findings.json instead, so a replicate never overwrites the published report."""
import os
import sys

from ccrcc.integrate import run_dir
from ccrcc.paths import RUNS
from ccrcc.render import render

if __name__ == "__main__":
    if len(sys.argv) > 1:
        print(render(run_dir(sys.argv[1]) / "findings.json"))
    elif os.environ.get("CCRCC_RUN_DIR"):
        d = run_dir()
        print(render(d / "findings.json", out=d / "report.html"))
    else:
        print(render(sorted(RUNS.glob("[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]/findings.json"))[-1]))
