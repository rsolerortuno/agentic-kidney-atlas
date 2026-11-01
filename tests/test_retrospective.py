"""Clinical eligibility and first observed trial date answer different questions."""
import pandas as pd

from ccrcc.outcomes import clinical_targets_from_records


def test_first_trial_includes_phase_one_and_unknown_dates_are_not_late():
    rows = []
    for gene, early in [("A", "2010-01-01"), ("B", "2017-01-01"), ("C", None)]:
        for i, (phase, date) in enumerate([("PHASE_1", early), ("PHASE_2", "2018-01-01"),
                                           ("PHASE_3", "2019-01-01")]):
            rows.append(dict(targetId=gene, drugId=gene, clinicalReportId=f"{gene}{i}",
                             clinicalStage=phase, studyStartDate=date, diseaseFromSource="renal cell carcinoma"))
    result = clinical_targets_from_records(pd.DataFrame(rows), {g: g for g in "ABC"}).set_index("gene")
    assert result.loc["A", "first_trial"] == pd.Timestamp("2010-01-01")
    assert not result.loc["A", "holdout"]
    assert result.loc["B", "holdout"]
    assert not result.loc["C", "holdout"]
    assert result.loc["C", "date_incomplete"]
    assert (result["trials"] == 2).all()
