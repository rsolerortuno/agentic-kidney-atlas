"""List verified specialist results the narrative left out (advisory; writes RUN/coverage.json).

  uv run python scripts/coverage.py [YYYY-MM-DD]
"""
import sys

from ccrcc.coverage import coverage
from ccrcc.integrate import current_run, run_dir

if __name__ == "__main__":
    d = run_dir(sys.argv[1]) if len(sys.argv) > 1 else current_run()
    rep = coverage(d)
    miss = rep["key_findings_not_cited"]
    print(f"key findings: {rep['key_findings']}; not cited in the narrative: {len(miss)}")
    for m in miss:
        note = " (gene named in limitations)" if m["explained_in_limitations"] else ""
        print(f"  - [{m['agent']}] {m['gene']}{note}: {m['statement']}")
    print("evidence types cited by agents but not by the narrative:", ", ".join(rep["evidence_types_not_cited"]) or "none")
