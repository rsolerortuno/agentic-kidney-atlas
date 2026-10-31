"""Structure and pocket context: AlphaFold DB models, cavity detection and molecular-glue-degrader (MGD) scoring.

Deterministic, no LLM, no local docking software beyond the pyKVFinder wheel:

1. The AlphaFold DB model of each candidate's UniProt accession is downloaded (cached on disk) with its
   per-residue pLDDT (the PDB B-factor column).
2. Residues below `structure.min_plddt_keep` are removed: in a model of a mostly disordered protein they are
   ribbons of spaghetti that create large, meaningless "cavities".
3. pyKVFinder detects cavities on the trimmed model (volume, area, maximum depth, average hydropathy, lining
   residues). Cavities are reported in UniProt residue numbering.
4. Each cavity gets an MGD site score in [0, 1] from four components (weights in config.yaml):
   - confidence: mean pLDDT of the lining residues (50 -> 0, 90 -> 1). A pocket in disordered model regions
     cannot be designed against.
   - geometry: a volume inside `pocket_volume` and a buried, deep cavity (molecular glues bind shallow but
     well-defined clefts; this prefers real clefts over surface dents and over huge interior voids).
   - lysines: solvent-accessible lysines within `lysine_radius` of the cavity centroid. A glue that recruits an
     E3 must leave a lysine the E2~ubiquitin can reach; with none, ubiquitin transfer has no acceptor.
   - degron: an exposed, non-helical stretch (loop or beta-hairpin) within `hairpin_radius`. Known glue degrons
     (the IMiD zinc-finger G-loop, the CK1alpha / GSPT1 beta-hairpin) are of that shape.
   Secondary structure is assigned geometrically from CA distances (P-SEA-style), and relative solvent
   accessibility with Biopython's Shrake-Rupley, so no DSSP binary is required.
5. Validation: PDBe-KB ligand-binding residues (experimental holo structures, UniProt numbering) are used as a
   positive control — how many of them fall in a detected cavity, and in which rank.

This is a prioritization of *sites* for a chemistry campaign, not a prediction that a degrader exists.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import requests
from Bio.PDB import PDBParser, ShrakeRupley

from .config import load
from .graph import GraphStore
from .paths import CACHE, PROCESSED

AF_API = "https://alphafold.ebi.ac.uk/api/prediction/"
PDBE_LIGAND_SITES = "https://www.ebi.ac.uk/pdbe/graph-api/uniprot/ligand_sites/"
AF_DIR = CACHE / "alphafold"
UA = {"User-Agent": "ccrcc-graphrag/0.1 (research; multi-omics GraphRAG case study)"}
# Tien et al. 2013 theoretical maximum side-chain + backbone SASA per residue (A^2)
MAX_ASA = {"ALA": 129, "ARG": 274, "ASN": 195, "ASP": 193, "CYS": 167, "GLN": 225, "GLU": 223, "GLY": 104,
           "HIS": 224, "ILE": 197, "LEU": 201, "LYS": 236, "MET": 224, "PHE": 240, "PRO": 159, "SER": 155,
           "THR": 172, "TRP": 285, "TYR": 263, "VAL": 174}


def _get_json(url: str, cache_name: str) -> dict | list | None:
    AF_DIR.mkdir(parents=True, exist_ok=True)
    path = AF_DIR / cache_name
    if path.exists():
        return json.loads(path.read_text(encoding="utf-8"))
    r = requests.get(url, timeout=120, headers=UA)
    if r.status_code == 404:
        path.write_text("null", encoding="utf-8")
        return None
    r.raise_for_status()
    path.write_text(r.text, encoding="utf-8")
    return r.json()


def alphafold_model(accession: str) -> dict | None:
    """Download (once) the AlphaFold DB model for an accession. Returns metadata + local PDB path."""
    meta = _get_json(AF_API + accession, f"{accession}.api.json")
    if not meta:
        return None
    m = meta[0]
    pdb = AF_DIR / f"{accession}.pdb"
    if not pdb.exists():
        r = requests.get(m["pdbUrl"], timeout=300, headers=UA)
        r.raise_for_status()
        pdb.write_bytes(r.content)
    return {"accession": accession, "path": pdb, "model_entity": m["modelEntityId"], "version": m["latestVersion"],
            "mean_plddt": m["globalMetricValue"], "created": m["modelCreatedDate"][:10],
            "fraction_plddt_confident": round(m["fractionPlddtConfident"] + m["fractionPlddtVeryHigh"], 3),
            "length": m["uniprotEnd"], "url": f"https://alphafold.ebi.ac.uk/entry/{accession}",
            "model_url": m["pdbUrl"]}


# ------------------------------------------------------------------ structure parsing
def read_structure(path: Path) -> dict:
    s = PDBParser(QUIET=True).get_structure("af", str(path))
    chain = next(iter(next(iter(s))))
    res = [r for r in chain if r.id[0] == " " and "CA" in r]
    ShrakeRupley().compute(chain, level="R")
    out = {
        "numbers": np.array([r.id[1] for r in res]),
        "names": [r.get_resname() for r in res],
        "ca": np.array([r["CA"].coord for r in res], dtype=float),
        "plddt": np.array([float(np.mean([a.bfactor for a in r])) for r in res]),
        "rsa": np.array([r.sasa / MAX_ASA.get(r.get_resname(), 200) for r in res]),
    }
    out["ss"] = secondary_structure(out["ca"])
    out["residues"] = res
    return out


def secondary_structure(ca: np.ndarray) -> np.ndarray:
    """P-SEA-style geometric assignment from CA distances: 'H' helix, 'E' strand, 'C' coil."""
    n = len(ca)
    ss = np.full(n, "C", dtype="<U1")

    def d(i, j):
        return float(np.linalg.norm(ca[i] - ca[j]))

    for i in range(n - 4):
        d2, d3, d4 = d(i, i + 2), d(i, i + 3), d(i, i + 4)
        if abs(d3 - 5.3) < 0.5 and abs(d4 - 6.4) < 0.7 and abs(d2 - 5.5) < 0.6:
            ss[i:i + 5] = "H"
        elif abs(d2 - 6.7) < 0.7 and abs(d3 - 9.9) < 1.0 and abs(d4 - 12.4) < 1.2 and not (ss[i:i + 5] == "H").any():
            ss[i:i + 5] = "E"
    return ss


def trimmed_pdb(path: Path, st: dict, min_plddt: float, out: Path) -> int:
    keep = {int(nr) for nr, p in zip(st["numbers"], st["plddt"]) if p >= min_plddt}
    lines = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(("ATOM", "HETATM")):
            if int(line[22:26]) in keep:
                lines.append(line)
        elif line.startswith(("HEADER", "TITLE", "SEQRES", "TER", "END")):
            lines.append(line)
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return len(keep)


# ------------------------------------------------------------------ cavities and MGD scoring
def _fit(x: float, lo: float, hi: float) -> float:
    """1 inside [lo, hi], decaying smoothly outside (log scale)."""
    if x <= 0:
        return 0.0
    if lo <= x <= hi:
        return 1.0
    ref, span = (lo, math.log(lo / max(x, 1e-6))) if x < lo else (hi, math.log(x / hi))
    return float(max(0.0, 1 - span / math.log(3)))


def cavity_table(pdb: Path, st: dict) -> list[dict]:
    import pyKVFinder

    cfg = load()["structure"]
    r = pyKVFinder.run_workflow(str(pdb), include_depth=True, include_hydropathy=True)
    if not r or not r.ncav:
        return []
    by_number = {int(n): i for i, n in enumerate(st["numbers"])}
    lys = np.array([i for i, nm in enumerate(st["names"]) if nm == "LYS"])
    out = []
    for cav, vol in r.volume.items():
        nums = sorted({int(x[0]) for x in r.residues[cav] if int(x[0]) in by_number})
        if not nums:
            continue
        idx = np.array([by_number[n] for n in nums])
        centroid = st["ca"][idx].mean(axis=0)
        dist = np.linalg.norm(st["ca"] - centroid, axis=1)
        near_lys = [int(st["numbers"][i]) for i in lys
                    if dist[i] <= cfg["lysine_radius"] and st["rsa"][i] >= cfg["lysine_min_rsa"]]
        near = np.flatnonzero(dist <= cfg["hairpin_radius"])
        loop = [i for i in near if st["ss"][i] != "H" and st["rsa"][i] >= 0.3]
        runs, run = [], []
        for i in sorted(loop):
            if run and i == run[-1] + 1:
                run.append(i)
            else:
                runs.append(run)
                run = [i]
        runs.append(run)
        best_run = max((len(x) for x in runs), default=0)
        comp = {
            "confidence": float(np.clip((st["plddt"][idx].mean() - 50) / 40, 0, 1)),
            "geometry": _fit(vol, *cfg["pocket_volume"]) * float(np.clip(r.max_depth[cav] / 10, 0, 1)),
            "lysines": float(min(len(near_lys), cfg["lysine_saturation"]) / cfg["lysine_saturation"]),
            "degron": float(min(best_run, 5) / 5),
        }
        w = cfg["mgd_weights"]
        out.append({
            "cavity": cav, "volume_a3": round(float(vol), 1), "area_a2": round(float(r.area[cav]), 1),
            "max_depth_a": round(float(r.max_depth[cav]), 2), "avg_hydropathy": round(float(r.avg_hydropathy[cav]), 3),
            "n_residues": len(nums), "residues": nums, "centroid": [round(float(c), 2) for c in centroid],
            "mean_plddt": round(float(st["plddt"][idx].mean()), 1),
            "accessible_lysines": near_lys, "exposed_loop_run": int(best_run),
            "components": {k: round(v, 3) for k, v in comp.items()},
            "mgd_score": round(sum(w[k] * comp[k] for k in w), 4),
        })
    out.sort(key=lambda c: -c["mgd_score"])
    for i, c in enumerate(out, 1):
        c["rank"] = i
    return out


def ligand_sites(accession: str) -> dict:
    """Experimental ligand-binding residues from PDBe-KB (UniProt numbering) for validation."""
    data = _get_json(PDBE_LIGAND_SITES + accession, f"{accession}.ligand_sites.json")
    if not data or accession not in data:
        return {"ligands": 0, "residues": []}
    res, ligs = set(), []
    for lig in data[accession]["data"]:
        ligs.append(lig["accession"])
        for r in lig["residues"]:
            res.update(range(int(r["startIndex"]), int(r["endIndex"]) + 1))
    return {"ligands": len(ligs), "ligand_ids": sorted(set(ligs))[:10], "residues": sorted(res)}


def validate(cavities: list[dict], sites: dict, top: int) -> dict:
    """Two separate questions, using experimentally ligand-bound residues (PDBe-KB) as ground truth.

    detection: do the detected cavities (all of them) cover the known ligand residues? (cavity detector quality)
    ranking: where does the MGD score place the best-overlapping cavity, and how much is in the top `top`?
    (score quality: the score is not trained on this, so a poor rank is a real weakness, not a bug)."""
    if not sites["residues"] or not cavities:
        return {"checked": False}
    known = set(sites["residues"])
    overlaps = {c["rank"]: len(known & set(c["residues"])) for c in cavities}
    best = min((r for r, n in overlaps.items() if n), default=None)
    covered = set().union(*[set(c["residues"]) for c in cavities]) & known
    covered_top = set().union(*[set(c["residues"]) for c in cavities[:top]]) & known
    return {"checked": True, "n_ligands": sites["ligands"], "n_known_residues": len(known),
            "n_cavities": len(cavities), "recovered_residues": len(covered),
            "recall": round(len(covered) / len(known), 3),
            "recall_in_top_pockets": round(len(covered_top) / len(known), 3),
            "best_cavity_rank": best,
            "best_cavity_rank_fraction": None if best is None else round(best / len(cavities), 3),
            "overlap_by_rank": {r: n for r, n in overlaps.items() if n}}


def run() -> dict:
    cfg = load()["structure"]
    targets = load()["cell_context"]["targets"]
    store = GraphStore()
    AF_DIR.mkdir(parents=True, exist_ok=True)
    results, rows = [], []
    for gene in targets:
        node = store.resolve(gene)
        acc = (node or {}).get("uniprot")
        acc = acc.split("|")[0].strip() if isinstance(acc, str) and acc else None
        if not acc:
            results.append({"gene": gene, "modeled": False, "reason": "no UniProt accession in the graph"})
            continue
        model = alphafold_model(acc)
        if not model:
            results.append({"gene": gene, "modeled": False, "reason": "no AlphaFold DB entry", "uniprot": acc})
            continue
        st = read_structure(model["path"])
        trimmed = AF_DIR / f"{acc}.plddt{int(cfg['min_plddt_keep'])}.pdb"
        kept = trimmed_pdb(model["path"], st, cfg["min_plddt_keep"], trimmed)
        all_cav = cavity_table(trimmed, st) if kept >= 30 else []
        cav = all_cav[: cfg["top_pockets"]]
        sites = ligand_sites(acc)
        entry = {"gene": gene, "modeled": True, "uniprot": acc,
                 **{k: v for k, v in model.items() if k not in ("path",)},
                 "residues_kept_above_plddt": kept,
                 "fraction_confident_residues": round(float((st["plddt"] >= cfg["confident_plddt"]).mean()), 3),
                 "pockets": cav, "n_cavities": len(all_cav),
                 "ligand_site_validation": validate(all_cav, sites, cfg["top_pockets"])}
        entry["best_mgd_score"] = cav[0]["mgd_score"] if cav else None
        results.append(entry)
        for c in cav:
            rows.append({"gene": gene, "uniprot": acc, "rank": c["rank"], "cavity": c["cavity"],
                         "volume_a3": c["volume_a3"], "area_a2": c["area_a2"], "max_depth_a": c["max_depth_a"],
                         "avg_hydropathy": c["avg_hydropathy"], "mean_plddt": c["mean_plddt"],
                         "n_accessible_lysines": len(c["accessible_lysines"]),
                         "exposed_loop_run": c["exposed_loop_run"], "mgd_score": c["mgd_score"],
                         **{f"comp_{k}": v for k, v in c["components"].items()},
                         "residues": ",".join(map(str, c["residues"]))})
    import pandas as pd
    pd.DataFrame(rows).to_parquet(PROCESSED / "structure_pockets.parquet", index=False)
    summary = {"modeled": sum(1 for r in results if r["modeled"]), "targets": len(targets),
               "config": {k: cfg[k] for k in ("min_plddt_keep", "pocket_volume", "lysine_radius", "mgd_weights")},
               "tool": "pyKVFinder (SES surface, probe_in 1.4 A, probe_out 4.0 A, step 0.6 A)",
               "results": results}
    (PROCESSED / "structure_summary.json").write_text(json.dumps(summary, indent=1))
    return summary


if __name__ == "__main__":
    s = run()
    print(f"{s['modeled']}/{s['targets']} modeled")
    for r in s["results"]:
        if not r["modeled"]:
            print(f"{r['gene']:10s} -- {r['reason']}")
            continue
        v = r["ligand_site_validation"]
        print(f"{r['gene']:10s} {r['uniprot']:8s} plddt {r['mean_plddt']:5.1f} kept {r['residues_kept_above_plddt']:4d} "
              f"best MGD {r['best_mgd_score']} | pockets " +
              ", ".join(f"{c['volume_a3']}A3/d{c['max_depth_a']}/K{len(c['accessible_lysines'])}/{c['mgd_score']}"
                        for c in r["pockets"]) +
              (f" | ligand-site detection recall {v['recall']} top3 {v['recall_in_top_pockets']} (best rank {v['best_cavity_rank']}/{v['n_cavities']})" if v.get("checked") else ""))
