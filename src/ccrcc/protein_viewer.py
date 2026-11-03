"""Package real coordinates and pocket-wall residues for the offline 3D explorer."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

from .paths import CACHE, PROCESSED


def build(summary: Path | None = None, models: Path | None = None, references: list | None = None) -> dict:
    summary = summary or PROCESSED / "structure_summary.json"
    models = models or CACHE / "alphafold"
    result = {"proteins": [], "unavailable": []}
    summary_data = json.loads(summary.read_text(encoding="utf-8")) if summary.exists() else {"results":[]}
    for entry in summary_data["results"]:
        accession = entry.get("uniprot", "")
        path = models / f"{accession}.pdb"
        if not entry.get("modeled") or not path.is_file():
            result["unavailable"].append(entry["gene"])
            continue
        pdb = path.read_text(encoding="utf-8")
        atoms = [line for line in pdb.splitlines() if line.startswith("ATOM  ")]
        # AlphaFold monomer files: preserve chain and UniProt numbering, fail on mismatch.
        chains = sorted({line[21:22].strip() for line in atoms})
        if len(chains) != 1:
            raise ValueError(f"Expected one chain in {path.name}; found {chains}")
        residues = {int(line[22:26]) for line in atoms}
        pockets = []
        for pocket in entry.get("pockets", []):
            missing = set(pocket["residues"]) - residues
            if missing:
                raise ValueError(f"Pocket residue numbering mismatch for {entry['gene']}: {sorted(missing)}")
            pockets.append({k: pocket[k] for k in ("rank", "cavity", "residues", "volume_a3",
                                                   "max_depth_a", "mean_plddt")})
        result["proteins"].append({"gene": entry["gene"], "accession": accession,
            "chain": chains[0], "pdb": "\n".join(atoms) + "\nEND\n", "pockets": pockets,
            "source_url": entry.get("url"), "model_url": entry.get("model_url"),
            "version": entry.get("version"), "model": entry.get("model_entity"),
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "detection_min_plddt": summary_data["config"]["min_plddt_keep"],
            "evidence_type": "Predicted monomer; computational cavity detection"})
    for ref in references or []:
        if ref.get("pdb"):
            result["proteins"].append({**ref, "experimental":True, "pockets":[],
                "chain":ref["target_chain"], "accession":ref["pdb_id"],
                "source_url":ref["source"], "sha256":ref["coordinates_sha256"]})
    return result


def gzip_b64(text: str) -> str:
    """gzip + base64 text for an inline <script> block; the page inflates it with the browser's DecompressionStream.

    The AlphaFold coordinates and the tissue viewer made the single-file report 36 MB; compressed it stays one
    self-contained file (opened from disk, no external requests) at a fraction of the size. base64 has no '<'."""
    import base64
    import gzip
    return base64.b64encode(gzip.compress(text.encode("utf-8"), compresslevel=9, mtime=0)).decode("ascii")


def script_json(data: dict) -> str:
    """Safe in an inline script, including untrusted structure metadata."""
    return json.dumps(data, separators=(",", ":"), ensure_ascii=True).replace("<", "\\u003c").replace(
        ">", "\\u003e").replace("&", "\\u0026")
