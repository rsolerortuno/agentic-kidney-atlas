"""M3 checks: every kept relation is verbatim-grounded with a PMID; the validator rejects bad extractions."""
import json

import pandas as pd
import pytest

from ccrcc import literature as lit
from ccrcc.paths import LITERATURE

REL = LITERATURE / "relations.parquet"


@pytest.mark.skipif(not REL.exists(), reason="run literature extraction first")
def test_kept_relations_are_verbatim_and_have_pmid():
    rel = pd.read_parquet(REL)
    pubs = pd.read_parquet(LITERATURE / "publications.parquet").set_index("pmid")
    assert len(rel) > 0
    assert rel["pmid"].notna().all()
    assert set(rel["relation"]) <= set(lit.RELATIONS)
    for _, r in rel.iterrows():
        text = lit._norm(pubs.loc[r["pmid"], "title"] + " " + pubs.loc[r["pmid"], "abstract"])
        assert lit._norm(r["literal_sentence"]) in text, r["relation_id"]


@pytest.mark.skipif(not REL.exists(), reason="run literature extraction first")
def test_extraction_cap_respected():
    pmids = set()
    for b in lit.BATCH_DIR.glob("batch_*.jsonl"):
        pmids |= {json.loads(x)["pmid"] for x in b.read_text(encoding="utf-8").splitlines() if x}
    assert len(pmids) <= 100


def test_validator_rejects_bad_relations(tmp_path, monkeypatch):
    lit_dir = tmp_path / "lit"
    (lit_dir / "extraction_batches").mkdir(parents=True)
    (lit_dir / "extracted").mkdir()
    pd.DataFrame([{"pmid": 1, "title": "CA9 in ccRCC.", "abstract": "CA9 is overexpressed in ccRCC tumors.",
                   "year": 2020, "journal": "J"}]).to_parquet(lit_dir / "publications.parquet")
    (lit_dir / "extraction_batches" / "batch_01.jsonl").write_text(json.dumps({"pmid": 1}) + "\n")
    good = {"pmid": 1, "subject": "CA9", "relation": "upregulated_in", "object": "ccRCC",
            "literal_sentence": "CA9 is overexpressed in ccRCC tumors.", "confidence": 0.9}
    rows = [good,
            {**good, "literal_sentence": "CA9 is massively overexpressed in ccRCC."},  # paraphrase
            {**good, "relation": "cures"},                                          # not in schema
            {**good, "pmid": 999},                                                  # pmid not in batch
            {**good, "confidence": 3}]
    (lit_dir / "extracted" / "x.jsonl").write_text("\n".join(json.dumps(r) for r in rows) + "\nnot json\n")
    monkeypatch.setattr(lit, "LITERATURE", lit_dir)
    monkeypatch.setattr(lit, "BATCH_DIR", lit_dir / "extraction_batches")
    monkeypatch.setattr(lit, "EXTRACTED_DIR", lit_dir / "extracted")
    res = lit.validate_extractions()
    assert res["kept"] == 1
    assert res["reasons"] == {"sentence_not_in_abstract": 1, "relation_not_allowed": 1, "pmid_not_in_batch": 1,
                              "bad_confidence": 1, "invalid_json": 1}
