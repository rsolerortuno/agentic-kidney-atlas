"""M6/M7 checks: final findings have 0 unsupported statements; the HTML report validates."""
import importlib.util

import pytest

from ccrcc.paths import REPORTS, ROOT, RUNS
from ccrcc.verify import verify_final

FINDINGS = sorted(RUNS.glob("[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]/findings.json"))
REPORT = REPORTS / "ccrcc_report.html"


@pytest.mark.skipif(not FINDINGS, reason="no orchestrated run yet")
def test_findings_have_no_unsupported_statements():
    rep = verify_final(FINDINGS[-1])
    assert rep["statements"] > 0
    assert rep["unsupported"] == 0, rep["detail"][:3]


@pytest.mark.skipif(not REPORT.exists(), reason="report not rendered yet")
def test_report_validates():
    spec = importlib.util.spec_from_file_location("validate_report", ROOT / "scripts" / "validate_report.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    assert mod.main(REPORT) == 0
