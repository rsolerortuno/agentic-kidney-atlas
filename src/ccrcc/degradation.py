"""Curated experimental MGD references; no generic degradability prediction."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd
import requests

from .paths import CACHE, PROCESSED

CATALOG = Path(__file__).parent / "reference_data" / "degradation.json"


def references() -> list[dict]:
    return json.loads(CATALOG.read_text(encoding="utf-8"))["references"]


def download_references():
    folder = CACHE / "ternary"
    folder.mkdir(parents=True, exist_ok=True)
    for ref in references():
        path = folder / (ref["pdb_id"] + ".pdb")
        if not path.exists():
            response = requests.get(ref["coordinates_url"], timeout=120)
            response.raise_for_status()
            # The frozen reference was normalized as UTF-8 text with LF newlines.
            content = response.text.replace("\r\n", "\n").encode("utf-8")
            if hashlib.sha256(content).hexdigest() != ref["coordinates_sha256"]:
                raise ValueError(f"Reference version changed: {ref['pdb_id']}; review before updating the catalog")
            path.write_bytes(content)
        validate_reference(ref, path)


def validate_reference(ref: dict, path: Path) -> str:
    if hashlib.sha256(path.read_bytes()).hexdigest() != ref["coordinates_sha256"]:
        raise ValueError(f"Reference checksum mismatch: {path.name}")
    lines = path.read_text(encoding="utf-8").splitlines()
    atom_lines = [l for l in lines if l.startswith(("ATOM  ", "HETATM"))]
    chains = {l[21] for l in atom_lines if l.startswith("ATOM")}
    if not {ref["target_chain"], ref["ligase_chain"], ref["adaptor_chain"]} <= chains:
        raise ValueError("Reference lacks an annotated protein chain")
    residues = {int(l[22:26]):l[17:20] for l in atom_lines if l.startswith("ATOM") and l[21] == ref["target_chain"]}
    for position, identity in ref["key_residues"].items():
        if residues.get(int(position)) != identity:
            raise ValueError(f"Degron residue identity mismatch: {ref['id']} {position}")
    if not any(l[17:20] == ref["ligand_resname"] and l[21] == ref["ligase_chain"] for l in atom_lines):
        raise ValueError("Reference lacks the annotated compound")
    selected_chains = {ref["target_chain"], ref["ligase_chain"], ref["adaptor_chain"]}
    selected = [l for l in atom_lines if l[21] in selected_chains and l[17:20] != "HOH"]
    return "\n".join(selected) + "\nEND\n"


def ligase_context() -> list[dict]:
    """Expression and alteration observations, never a claim of functional E3 activity."""
    genes = ["CRBN", "DDB1", "CUL4A", "RBX1", "VHL"]
    universe = pd.read_parquet(PROCESSED / "gene_universe.parquet").set_index("gene")
    mutations = pd.read_parquet(PROCESSED / "layer_mutations.parquet").set_index("gene")
    sc = pd.read_parquet(PROCESSED / "sc_celltype_expression.parquet")
    out = []
    for gene in genes:
        u = universe.loc[gene].to_dict() if gene in universe.index else {}
        m = mutations.loc[gene].to_dict() if gene in mutations.index else {}
        cells = sc[(sc["hgnc_symbol"] == gene) & (sc["tissue"] == "tumor")]
        out.append({"gene":gene, "protein_log2fc":u.get("protein_log2fc"), "protein_padj":u.get("protein_padj"),
                    "mutation_frequency":m.get("freq"), "copy_loss_frequency":u.get("loss_freq"),
                    "cell_expression":cells[["cell_group", "pct_expressing", "n_donors"]].to_dict("records"),
                    "functionality":"not established", "subcellular_colocalization":"not measured"})
    return out


def build() -> dict:
    refs = []
    for ref in references():
        path = CACHE / "ternary" / (ref["pdb_id"] + ".pdb")
        refs.append({**ref, "coordinates_available":path.is_file(),
                     "pdb": validate_reference(ref, path) if path.is_file() else None})
    return {"references":refs, "ligases":ligase_context(),
            "prediction_status":"No new MGD or ternary prediction is made. Reference matching is annotation, not a classifier."}


if __name__ == "__main__":
    download_references()
    print(json.dumps({"validated_references":[r["id"] for r in references()]}))
