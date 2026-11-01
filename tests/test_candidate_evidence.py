import pytest
from ccrcc.candidate_evidence import build, expression
from ccrcc.paths import PROCESSED


def test_nonsignificant_opposite_is_not_a_significant_contradiction():
    row = {"effect":-1.2,"q":.07}
    assert expression("replication",row,"effect","q","test",reference=1)["status"] == "non_favorable"
    row["q"] = .01
    assert expression("replication",row,"effect","q","test",reference=1)["status"] == "contradictory"
    row["effect"] = None
    assert expression("replication",row,"effect","q","test",reference=1)["status"] == "missing"


def test_missing_discovery_cannot_establish_replication_or_biomarker(tmp_path):
    assert expression("replication", {"effect":2, "q":.001}, "effect", "q", "test", replication=True)["status"] == "missing"
    card = build(["UNMEASURED"], tmp_path)[0]
    assert card["measured"] == 0
    assert all(m["status"] == "not_established" for m in card["modalities"])


@pytest.mark.skipif(not (PROCESSED / "confirmatory_protein.parquet").exists(), reason="requires real data")
def test_missing_ca9_is_not_failed_hk2_replication():
    cards = {c["gene"]: c for c in build(["CA9","HK2","CD70"])}
    obs = lambda gene,name: next(o for o in cards[gene]["observations"] if o["name"] == name)
    assert obs("CA9","Confirmatory protein")["status"] == "missing"
    assert obs("HK2","Confirmatory protein")["status"] == "non_favorable"
    assert obs("CD70","Cell attribution")["status"] == "contradictory"
    assert obs("CD70","Vital-organ staining")["status"] == "missing"
    assert all(next(m for m in c["modalities"] if m["name"] == "Degradation / MGD")["status"] == "not_established" for c in cards.values())
