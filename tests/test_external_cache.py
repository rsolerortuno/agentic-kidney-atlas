"""Agent-time lookups stay out of the frozen input tree; the highest drug stage is not an alphabetical maximum."""
import json

from ccrcc import external as ext
from ccrcc.discovery import max_stage


def test_max_stage_orders_by_clinical_progress():
    assert max_stage(["PHASE_4", "APPROVAL", "PHASE_1"]) == "APPROVAL"       # alphabetical maximum was PHASE_4
    assert max_stage(["PHASE_2", "PHASE_2_3", "PHASE_1_2"]) == "PHASE_2_3"
    assert max_stage(["WITHDRAWAL", "PHASE_1"]) == "PHASE_1"                  # a clinical stage wins over other labels
    assert max_stage(["WITHDRAWAL"]) == "WITHDRAWAL"
    assert max_stage([]) is None


def test_agent_lookups_do_not_write_into_the_frozen_cache(tmp_path, monkeypatch):
    frozen, agent = tmp_path / "cache", tmp_path / "cache_agent"
    (frozen / "string").mkdir(parents=True)
    monkeypatch.setattr(ext, "CACHE", frozen)
    monkeypatch.setattr(ext, "WRITE_ROOT", agent)
    calls = []

    def fetch():
        calls.append(1)
        return {"v": 1}

    assert ext._cached("string", "new-key", fetch)["data"] == {"v": 1}
    assert not list(frozen.rglob("*.json"))                 # nothing added to the frozen tree
    assert len(list(agent.rglob("*.json"))) == 1
    assert ext._cached("string", "new-key", fetch)["data"] == {"v": 1} and len(calls) == 1   # served from the agent cache

    # an answer already in the frozen cache is read from there, never fetched again or copied
    name = next(agent.rglob("*.json")).name
    (frozen / "string" / name).write_text(json.dumps({"key": "new-key", "data": {"v": 2}}), encoding="utf-8")
    assert ext._cached("string", "new-key", fetch)["data"] == {"v": 2} and len(calls) == 1


def test_mcp_server_redirects_writes():
    from ccrcc.mcp_servers import external as ext_mcp   # noqa: F401  (importing the server sets the write root)
    from ccrcc.paths import AGENT_CACHE, CACHE
    assert ext.WRITE_ROOT == AGENT_CACHE != CACHE
