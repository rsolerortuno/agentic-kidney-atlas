"""M2 checks: statistics behave correctly and known ccRCC biology is recovered."""
import numpy as np
import pandas as pd
import pytest

from ccrcc.omics import cnv_frequencies, differential
from ccrcc.paths import PROCESSED

DE = PROCESSED / "de_results.parquet"
HIF_TARGETS = ["CA9", "NDUFA4L2", "ANGPTL4", "EGLN3", "VEGFA"]


def test_differential_paired_detects_shift():
    rng = np.random.default_rng(0)
    cases = [f"C{i}" for i in range(20)]
    normal = pd.DataFrame(rng.normal(0, 1, (2, 20)), index=["UP", "NULL"], columns=cases)
    tumor = normal + rng.normal(0, 0.1, (2, 20))
    tumor.loc["UP"] += 2.0
    de = differential(tumor, normal).set_index("gene")
    assert de.loc["UP", "method"] == "paired_t"
    assert de.loc["UP", "padj"] < 1e-6 and de.loc["UP", "log2fc"] == pytest.approx(2.0, abs=0.1)
    assert de.loc["NULL", "padj"] > 0.01


def test_cnv_frequencies():
    cnv = pd.DataFrame([[-0.5, -0.3, 0.0, 0.5]], index=["G"], columns=list("abcd"))
    f = cnv_frequencies(cnv).iloc[0]
    assert f["loss_freq"] == 0.5 and f["gain_freq"] == 0.25


@pytest.mark.skipif(not DE.exists(), reason="run `ccrcc omics` first")
def test_hif_targets_up_in_rna_and_protein_when_measured():
    de = pd.read_parquet(DE)
    for layer in ("rna", "protein"):
        sub = de[(de["layer"] == layer) & de["gene"].isin(HIF_TARGETS)].set_index("gene")
        if layer == "rna":
            assert set(sub.index) == set(HIF_TARGETS)
        else:
            # EGLN3 is absent from the CPTAC proteome; VEGFA has 14 tumor / 4 normal values (not testable)
            assert {"CA9", "NDUFA4L2", "ANGPTL4"} <= set(sub.index)
        assert (sub["padj"] < 0.05).all(), sub
        assert (sub["log2fc"] > 0).all(), sub


@pytest.mark.skipif(not (PROCESSED / "cnv_arm_summary.parquet").exists(), reason="run `ccrcc omics` first")
def test_3p_loss_is_most_frequent_arm_loss():
    arm = pd.read_parquet(PROCESSED / "cnv_arm_summary.parquet").set_index("arm")
    assert arm["median_loss_freq"].idxmax() == "3p"
    assert arm.loc["3p", "median_loss_freq"] > 0.5
