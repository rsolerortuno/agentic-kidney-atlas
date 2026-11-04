"""Offline checks of the v2 architecture benchmark (no model calls; only the frozen suite in the repository)."""
import json

import pytest

from benchmarks.architecture_v2.build import FROZEN, QUESTIONS, gold, verdict
from benchmarks.architecture_v2.evaluate import bootstrap
from benchmarks.architecture_v2.run import SPECIALISTS, batches, parse, prompt_full, prompt_specialist

CASES = json.loads((FROZEN / "cases.json").read_text(encoding="utf-8"))
GOLD = json.loads((FROZEN / "gold.json").read_text(encoding="utf-8"))


def test_frozen_suite_matches_its_manifest():
    import hashlib
    m = json.loads((FROZEN / "manifest.json").read_text(encoding="utf-8"))
    for name, digest in m["files"].items():
        assert hashlib.sha256((FROZEN / name).read_bytes()).hexdigest() == digest, name
    assert m["n_cases"] == len(CASES) == 24


def test_gold_is_the_rubric_applied_to_the_target_records():
    for c in CASES:
        target = [{k: v for k, v in r.items() if k not in ("id", "entity")} for r in c["records"]
                  if r["entity"] == c["target_entity"]]
        g = gold(target)
        assert all(g[q] == GOLD[c["id"]][q] for q in QUESTIONS + ["verdict"]), c["id"]


def test_dossiers_hide_gene_names_and_mix_in_other_entities():
    for c in CASES:
        assert len({r["entity"] for r in c["records"]}) == 3
        assert "unrelated" not in json.dumps(c)
        assert all(r["entity"].startswith("ENTITY_") for r in c["records"])


def test_verdict_classes_are_balanced_enough_to_discriminate():
    counts = {}
    for g in GOLD.values():
        counts[g["verdict"]] = counts.get(g["verdict"], 0) + 1
    assert max(counts.values()) / len(GOLD) < 0.35 and len(counts) == 4


def test_missing_data_is_insufficient_not_negative():
    assert gold([])["protein_up"] == "insufficient"
    assert gold([{"type": "normal_tissue_ihc", "measured": False, "vital_max_level": None}])["vital_safety_concern"] == "insufficient"
    shares = {"type": "single_cell_tumour", "groups": [{"cell_group": "RCC", "share_of_tumour_transcripts": 0.3, "pct_expressing": 80},
                                                       {"cell_group": "CD8_T", "share_of_tumour_transcripts": 0.5, "pct_expressing": 10}]}
    assert gold([shares])["tumour_intrinsic"] == "no"          # transcript share, not pct_expressing, decides


def test_verdict_rule_order():
    assert verdict(dict(protein_up="yes", replicated="yes", tumour_intrinsic="yes", surface_accessible="yes",
                        vital_safety_concern="no")) == "surface_target_candidate"
    assert verdict(dict(protein_up="no", replicated="yes", tumour_intrinsic="no")) == "loss_biomarker"
    assert verdict(dict(protein_up="yes", replicated="yes", tumour_intrinsic="no")) == "infiltration_signal"


def test_no_context_holds_two_variants_of_one_gene():
    for size in (4, 12):
        for rep in (1, 2, 3):
            bb = batches(CASES, rep, size)
            assert sorted(c["id"] for b in bb for c in b) == sorted(c["id"] for c in CASES)
            assert all(len({c["id"][:3] for c in b}) == len(b) for b in bb)


def test_specialists_see_only_their_record_types():
    batch = CASES[:2]
    for name, (types, qs) in SPECIALISTS.items():
        p = prompt_specialist(batch, name)
        payload = json.loads(p[p.index("CASES\n") + 6:])
        assert {r["type"] for c in payload for r in c["records"]} <= set(types)
    assert sorted(q for _, qs in SPECIALISTS.values() for q in qs) == sorted(QUESTIONS)
    assert "verdict" in prompt_full(batch, "role")


def test_parser_takes_first_json_object_and_rejects_missing_cases():
    a = {"id": "x", **{q: "yes" for q in QUESTIONS}, "verdict": "loss_biomarker"}
    assert parse('{"answers": [' + json.dumps(a) + ']} trailing text', ["x"], QUESTIONS, True)[0]["id"] == "x"
    with pytest.raises(ValueError):
        parse('{"answers": []}', ["x"], QUESTIONS, True)
    with pytest.raises(ValueError):
        parse('{"answers": [' + json.dumps({**a, "verdict": "maybe"}) + ']}', ["x"], QUESTIONS, True)


def test_bootstrap_pairs_by_gene():
    rows = [{"gene": g, "arm": arm, "v": float(arm == "B")} for g in range(4) for arm in "BC" for _ in range(3)]
    out = bootstrap(rows, "C", "B", "v", n=200)
    assert out["genes"] == 4 and out["diff"] == -1.0 and out["ci95"] == [-1.0, -1.0]


def test_structured_schema_fixes_ids_labels_and_count():
    from benchmarks.architecture_v2.run import json_schema
    s = json_schema(["a", "b"], ["protein_up"], True)
    arr = s["properties"]["answers"]
    assert arr["minItems"] == arr["maxItems"] == 2
    item = arr["items"]
    assert item["properties"]["id"]["enum"] == ["a", "b"]
    assert item["properties"]["protein_up"]["enum"] == ["insufficient", "no", "yes"]
    assert "verdict" in item["required"] and item["additionalProperties"] is False
