"""M1 check: downloaded files match the LinkedOmics CPTAC-CCRCC dimensions."""
import json

import pytest

from ccrcc.paths import RAW

MANIFEST = RAW / "manifest.json"
pytestmark = pytest.mark.skipif(not MANIFEST.exists(), reason="run `ccrcc download` first")

EXPECTED = {
    "HS_CPTAC_CCRCC_CNV_gene_Tumor.cct": 110,
    "HS_CPTAC_CCRCC_proteome_Tumor.cct": 110,
    "HS_CPTAC_CCRCC_proteome_Normal.cct": 84,
    "HS_CPTAC_CCRCC_phosphoproteome_gene_Tumor.cct": 110,
    "HS_CPTAC_CCRCC_phosphoproteome_gene_Normal.cct": 84,
    "HS_CPTAC_CCRCC_RNAseq_fpkm_log2_Tumor.cct": 110,
    "HS_CPTAC_CCRCC_RNAseq_fpkm_log2_Normal.cct": 75,
}


def test_manifest_dimensions_and_hashes():
    files = {e["file"]: e for e in json.loads(MANIFEST.read_text())["files"]}
    assert len(files) == 9
    for name, n_samples in EXPECTED.items():
        assert files[name]["columns"] == n_samples, name
        assert len(files[name]["sha256"]) == 64
    assert files["hgnc_complete_set.txt"]["rows"] > 40000
