import pytest
from ccrcc import semantic

EV = [{"gene": "CA9", "log2fc": 2.27, "padj": 0.001, "n_pairs": 75}]


@pytest.mark.parametrize("text", ["log2FC -2.27", "log2FC 227", "log2FC 0.001",
                                  "padj 2.27", "log2FC -2", "unknown metric 2.27"])
def test_reject_sign_unit_wrong_field_and_unbound_values(text):
    assert semantic.check(text, "CA9", EV)["unmatched_numbers"]


def test_exact_field_rounding_and_no_other_gene():
    assert not semantic.check("log2FC 2.27, padj 0.001, 75 pairs", "CA9", EV)["unmatched_numbers"]
    assert semantic.check("log2FC 2.27", "VHL", EV)["unmatched_numbers"]


def test_extra_ranking_numbers_cannot_rescue_a_wrong_metric():
    assert semantic.check("log2FC 9.99", "CA9", EV, {9.99})["unmatched_numbers"]


def test_ambiguous_layers_and_invalid_percent_are_rejected():
    assert semantic.check("RNA log2FC 2.27", "CA9", EV + [{**EV[0], "layer":"protein"}])["unmatched_numbers"]
    assert semantic.check("log2FC 227%", "CA9", EV)["unmatched_numbers"]
    assert semantic.check("padj 0.1%", "CA9", EV)["unmatched_numbers"]
    assert semantic.check("Measured in 3 samples", "CA9", EV)["unmatched_numbers"]
    assert semantic.check("RNA log2FC 2.27", "CA9", [{**EV[0], "layer":"protein"}])["unmatched_numbers"]


NESTED = {"evidence_id":"EV:test", "gene":"CA9", "type":"treatment", "source":"trial1",
          "cohort":"validation", "arms":{"atezo":{"pvalue":.0317, "log_hr":-.31}}, "freq":.127}
CLAIM = {"evidence_id":"EV:test", "field":["arms","atezo","log_hr"], "display":"-0.31", "unit":"native"}


def test_nested_metric_binding_preserves_context_and_sign():
    text = semantic.numeric_claim_text(CLAIM, NESTED)
    assert "trial1" in text and "validation" in text and "arms.atezo.log_hr=-0.31" in text
    result = semantic.check(text, "CA9", [NESTED], numeric_claims=[CLAIM])
    assert not result["unmatched_numbers"] and not result["numeric_binding_errors"]
    for bad in (dict(CLAIM, display="0.31"), dict(CLAIM, field=["arms","atezo","pvalue"]),
                dict(CLAIM, evidence_id="EV:uncited"), dict(CLAIM, unit="percent")):
        assert semantic.check(text, "CA9", [NESTED], numeric_claims=[bad])["numeric_binding_errors"]
    for changed in (text.replace("trial1", "trial2"), text.replace("validation", "discovery"), text + text):
        assert semantic.check(changed, "CA9", [NESTED], numeric_claims=[CLAIM])["numeric_binding_errors"]


def test_percent_is_explicit_and_scientific_precision_is_supported():
    claim = dict(CLAIM, field=["freq"], display="1.27e1", unit="percent")
    text = semantic.numeric_claim_text(claim, NESTED)
    assert not semantic.check(text, "CA9", [NESTED], numeric_claims=[claim])["numeric_binding_errors"]
    assert semantic.check(text, "VHL", [NESTED], numeric_claims=[claim])["numeric_binding_errors"]


def test_singlecell_percentages_are_not_multiplied_twice():
    evidence = {**NESTED, "pct_expressing":54.7, "cell_group":"RCC"}
    claim = dict(CLAIM, field=["pct_expressing"], display="54.7", unit="native_percent")
    annotation = semantic.numeric_claim_text(claim, evidence)
    assert 'cell_group=RCC' in annotation and '54.7%' in annotation
    assert not semantic.check(annotation, 'CA9', [evidence], numeric_claims=[claim])["numeric_binding_errors"]
    with pytest.raises(ValueError, match="fraction"):
        semantic.numeric_claim_text(dict(claim, display='5470', unit='percent'), evidence)


def test_survival_direction_is_not_tumour_normal_expression_direction():
    evidence=[{'gene':'ALDOB','type':'omics_de','significant':True,'log2fc':-2.27}]
    assert not semantic.check('Higher expression is associated with better survival.', 'ALDOB', evidence)['direction_contradiction']
    assert semantic.check('ALDOB is upregulated in tumour. Higher expression is associated with better survival.', 'ALDOB', evidence)['direction_contradiction']


@pytest.mark.parametrize("field", [["arms","missing"], ["gene"], [], [True], ["arms",0]])
def test_invalid_metric_paths_fail_closed(field):
    with pytest.raises(ValueError):
        semantic.numeric_claim_text(dict(CLAIM, field=field), NESTED)


def test_verifier_passes_bindings_and_final_does_not_hide_filtered_failures(tmp_path, monkeypatch):
    import json
    from ccrcc import verify
    class Store:
        def check_evidence_ids(self, ids):
            return {"valid":["EV:test"]}
        def get_evidence(self, ids):
            return [NESTED]
    v = verify.Verifier(store=Store())
    monkeypatch.setattr(v, "_gene_evidence", lambda gene: {"EV:test"})
    statement = {"gene":"CA9", "evidence_ids":["EV:test"],
                 "statement":semantic.numeric_claim_text(CLAIM, NESTED), "numeric_claims":[CLAIM]}
    assert v.check_statement(statement)[0] is not None
    statement["numeric_claims"] = [dict(CLAIM, unit="percent")]
    assert v.check_statement(statement)[1]["semantic_failure"]
    monkeypatch.setattr(verify, "Verifier", lambda **kwargs: v)
    monkeypatch.setattr(verify.provenance, "load_basis", lambda path: None)
    path = tmp_path / "findings.json"
    path.write_text(json.dumps({"dropped_during_assembly":[{"where":"summary", "kept":False,
        "semantic_failure":True, "text":"Incorrect signed effect"}]}))
    report = verify.verify_final(path)
    assert report["unsupported"] == 1
    assert report["interpretation"] == "Automatic checks detected failures"
