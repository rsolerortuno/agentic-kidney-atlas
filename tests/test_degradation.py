from copy import deepcopy
import pytest
from ccrcc.degradation import references, validate_reference
from ccrcc.paths import CACHE


def test_catalog_never_claims_candidate_transfer():
    for ref in references():
        assert ref["gene"] and ref["ligase"] and ref["compound"]
        assert ref["evidence_type"] == "experimental_reference"
        assert "no ccRCC" in ref["scope"]
        assert ref["numbering"].startswith("PDB author")


@pytest.mark.skipif(not (CACHE / "ternary/5HXB.pdb").exists(), reason="requires downloaded references")
def test_experimental_chain_compound_and_residue_identity():
    for ref in references():
        path = CACHE / "ternary" / (ref["pdb_id"]+".pdb")
        pdb = validate_reference(ref,path)
        assert ref["ligand_resname"] in pdb
        if ref["key_residues"]:
            wrong = deepcopy(ref)
            wrong["key_residues"][next(iter(wrong["key_residues"]))] = "XXX"
            with pytest.raises(ValueError,match="identity mismatch"):
                validate_reference(wrong,path)
