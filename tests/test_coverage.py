import json

from ccrcc.coverage import coverage


def _write(p, obj):
    p.write_text(json.dumps(obj), encoding="utf-8")


def test_key_finding_not_cited_is_listed(tmp_path):
    _write(tmp_path / "agent_omics.verified.json", {"findings": [
        {"gene": "HK2", "key": True, "statement": "HK2 does not replicate.", "evidence_ids": ["EV:confirmatory_protein:a"]},
        {"gene": "CA9", "key": True, "statement": "CA9 up.", "evidence_ids": ["EV:omics_de:b"]},
        {"gene": "NNMT", "statement": "NNMT up.", "evidence_ids": ["EV:prism_response:c"]}]})
    _write(tmp_path / "narrative.json", {"summary": [{"text": "CA9", "evidence_ids": ["EV:omics_de:b"]}],
                                         "sections": {}, "targets": [], "limitations": [{"text": "HK2 left out"}]})
    c = coverage(tmp_path)
    assert c["key_findings"] == 2
    assert [m["gene"] for m in c["key_findings_not_cited"]] == ["HK2"]
    assert c["key_findings_not_cited"][0]["explained_in_limitations"]
    assert c["evidence_types_not_cited"] == ["confirmatory_protein", "prism_response"]
    assert c["share_of_finding_evidence_cited"] == round(1 / 3, 3)
    assert (tmp_path / "coverage.json").exists()
