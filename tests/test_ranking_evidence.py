from ccrcc.graph import GraphStore
from ccrcc.ranking_evidence import load
from ccrcc.semantic import numeric_claim_text, check


def test_ranking_metrics_are_bound_to_profile_and_content(tmp_path):
    store = GraphStore(tmp_path / "graph.sqlite")
    store.add_node("HGNC:1", "Gene", "CA9")
    conv = {"ranking":[{"gene":"CA9", "rank":2, "score":3.5}],
            "profiles":{"target":[{"gene":"CA9", "rank":1, "score":4.0}]}}
    ids = load(conv, store)
    ev = next(e for e in store.get_evidence(ids) if e["profile"] == "target")
    claim = {"evidence_id":ev["evidence_id"], "field":["rank"], "display":"1", "unit":"native"}
    text = numeric_claim_text(claim, ev)
    assert not check(text, "CA9", [ev], numeric_claims=[claim])["numeric_binding_errors"]
    assert check(text.replace("profile=target", "profile=convergence"), "CA9", [ev], numeric_claims=[claim])["numeric_binding_errors"]
    conv["profiles"]["target"][0]["rank"] = 3
    assert load(conv, store) != ids
