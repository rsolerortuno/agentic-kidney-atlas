"""Render reports/ccrcc_report.html from the latest runs/<date>/findings.json (or a given date)."""
import sys

from ccrcc.integrate import run_dir
from ccrcc.paths import RUNS
from ccrcc.render import render

if __name__ == "__main__":
    if len(sys.argv) > 1:
        path = run_dir(sys.argv[1]) / "findings.json"
    else:
        path = sorted(RUNS.glob("*/findings.json"))[-1]
    print(render(path))
