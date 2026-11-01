import json
import pytest
from ccrcc.protein_viewer import build, script_json
from ccrcc.paths import PROCESSED


def test_script_json_cannot_close_script_tag():
    text = script_json({"value": "</script><script>alert(1)</script>"})
    assert "<" not in text
    assert json.loads(text)["value"].startswith("</script>")


@pytest.mark.skipif(not (PROCESSED / "structure_summary.json").exists(), reason="local structures required")
def test_actual_coordinates_cover_every_highlighted_residue():
    data = build()
    assert data["proteins"]
    for protein in data["proteins"]:
        residues = {int(line[22:26]) for line in protein["pdb"].splitlines() if line.startswith("ATOM")}
        assert len(protein["sha256"]) == 64
        for pocket in protein["pockets"]:
            assert set(pocket["residues"]) <= residues
            assert pocket["volume_a3"] > 0
