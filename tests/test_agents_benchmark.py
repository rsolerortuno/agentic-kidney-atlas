from ccrcc.agents_benchmark import ARMS, CHECK_LABELS, build


def test_agents_benchmark_reads_result_files():
    d = build()
    assert d, "benchmark result files are missing"
    keys = {c["key"] for c in d["conditions"]}
    assert {"runs", "runs_first_attempt", "runs_structured"} <= keys
    for c in d["conditions"]:
        assert set(c["arms"]) == set(ARMS)
        for a in c["arms"].values():
            assert 0 <= a["accuracy"] <= 100 and a["failed"] <= a["calls"]
    arms = {x["arm"] for x in d["pipeline"]}
    assert {"multi", "multi_improved", "single"} <= arms
    assert set(d["checks"]) <= set(CHECK_LABELS)
