import json
import pytest
from ccrcc.reproducibility import sha256, verify, require_verified


def test_detects_same_size_corruption_and_missing_file(tmp_path):
    data = tmp_path / "data"
    data.mkdir()
    f = data / "sample.bin"
    f.write_bytes(b"original")
    manifest = tmp_path / "inputs.json"
    manifest.write_text(json.dumps({"files": [{"path": "sample.bin", "bytes": 8, "sha256": sha256(f)}]}))
    before = manifest.read_bytes()
    assert verify(manifest, data)["ok"]
    f.write_bytes(b"modified")
    assert verify(manifest, data)["errors"][0]["error"] == "checksum mismatch"
    f.unlink()
    assert verify(manifest, data)["errors"][0]["error"] == "missing"
    assert manifest.read_bytes() == before


def test_configuration_drift_and_preflight_failure_preserve_outputs(tmp_path, monkeypatch):
    from ccrcc import cli, reproducibility
    config = tmp_path / "config.yaml"
    config.write_text("weights: original")
    manifest = tmp_path / "inputs.json"
    manifest.write_text(json.dumps({"files":[], "config_sha256":sha256(config)}))
    assert verify(manifest, tmp_path, tmp_path)["ok"]
    config.write_text("weights: modified")
    with pytest.raises(RuntimeError, match="configuration"):
        require_verified(manifest, tmp_path, tmp_path)
    processed = tmp_path / "processed"
    processed.mkdir()
    sentinel = processed / "valuable-result.txt"
    sentinel.write_text("preserve me")
    monkeypatch.setattr(cli, "PROCESSED", processed)
    def fail():
        raise RuntimeError("frozen input failed")
    monkeypatch.setattr(reproducibility, "require_verified", fail)
    with pytest.raises(RuntimeError, match="frozen input failed"):
        cli.run_all(fresh=True)
    assert sentinel.read_text() == "preserve me"


def test_frozen_mode_blocks_http_and_restores_transport():
    import requests
    from ccrcc.reproducibility import offline_inputs
    original = requests.sessions.Session.request
    with pytest.raises(RuntimeError, match="Frozen reproduction attempted network"):
        with offline_inputs():
            requests.get("https://example.invalid/missing-input")
    assert requests.sessions.Session.request is original


def test_run_all_is_offline_by_default(monkeypatch):
    import requests
    from ccrcc import cli
    def pipeline(**kwargs):
        assert kwargs == {"fresh": False, "frozen": True}
        with pytest.raises(RuntimeError, match="Frozen reproduction attempted network"):
            requests.post("https://example.invalid/api", json={})
        return 0
    monkeypatch.setattr(cli, "_run_all", pipeline)
    assert cli.run_all(fresh=False) == 0


def test_extra_cached_input_and_changed_saved_analysis_are_rejected(tmp_path):
    data = tmp_path / "data"
    cache = data / "cache"
    cache.mkdir(parents=True)
    narrative = tmp_path / "narrative.json"
    narrative.write_text('{"claim":"original"}')
    manifest = tmp_path / "inputs.json"
    manifest.write_text(json.dumps({
        "files": [], "input_trees": ["cache"],
        "project_files": [{"path":"narrative.json", "sha256":sha256(narrative)}]
    }))
    assert verify(manifest, data, tmp_path)["ok"]
    (cache / "unexpected.json").write_text('{}')
    narrative.write_text('{"claim":"changed"}')
    errors = verify(manifest, data, tmp_path)["errors"]
    assert {e["error"] for e in errors} == {
        "unlisted input in frozen tree", "saved analysis input checksum mismatch"}


def test_frozen_baseline_never_invokes_model_for_missing_response(tmp_path, monkeypatch):
    from ccrcc import baselines
    monkeypatch.setattr(baselines, "OUT", tmp_path)
    monkeypatch.setattr(baselines, "_claude", lambda: "claude")
    def unexpected(*args, **kwargs):
        pytest.fail("Frozen reconstruction must not regenerate model responses")
    monkeypatch.setattr(baselines.subprocess, "run", unexpected)
    with pytest.raises(FileNotFoundError, match="Frozen model response missing"):
        baselines.closed_book(runs=1, cached_only=True)


def test_selected_run_requires_frozen_narrative_and_extra_replicates_fail(tmp_path):
    from ccrcc import baselines
    run = tmp_path / "runs" / "2026-11-01"
    run.mkdir(parents=True)
    narrative = run / "narrative.json"
    narrative.write_text('{}')
    manifest = tmp_path / "inputs.json"
    expected = {"files":[], "analysis_run":"runs/2026-11-01", "project_files":[],
                "project_input_globs":["runs/replicates/*/findings.json"]}
    manifest.write_text(json.dumps(expected))
    assert verify(manifest, tmp_path, tmp_path)["errors"][0]["error"] == "selected narrative is not frozen"
    expected["project_files"] = [{"path":"runs/2026-11-01/narrative.json", "sha256":sha256(narrative)}]
    manifest.write_text(json.dumps(expected))
    assert verify(manifest, tmp_path, tmp_path)["analysis_run"] == "runs/2026-11-01"
    extra = tmp_path / "runs" / "replicates" / "new" / "findings.json"
    extra.parent.mkdir(parents=True)
    extra.write_text('{}')
    assert verify(manifest, tmp_path, tmp_path)["errors"][0]["error"] == "unlisted saved analysis input"


def test_baseline_uses_selected_run_even_if_newer_directory_exists(tmp_path, monkeypatch):
    from ccrcc import baselines
    for name in ["2026-11-01", "2026-11-02"]:
        p=tmp_path / name
        p.mkdir()
        (p / "findings.json").write_text('{}')
    monkeypatch.setattr(baselines, "RUNS", tmp_path)
    chosen=tmp_path / "2026-11-01"
    assert baselines.run_dirs(official=chosen) == [chosen]
