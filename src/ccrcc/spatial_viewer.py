"""Data for the interactive tissue viewer in the report (H&E image + Visium spots + neighbours).

Deterministic, no LLM. Four Visium sections of Li et al. 2022 are published as single-section datasets that carry the
H&E image (two tumour-core and two tumour-normal interface sections). For each, this module writes
data/processed/spatial_viewer.json with:
- the H&E image cropped to the tissue and downscaled (JPEG, base64),
- every in-tissue spot: position on that image, deconvolved cell-group shares (from spatial_deconvolution.parquet,
  computed on the same spots) and log1p CPM of the analysed genes (quantized),
- ring-1 neighbours of every spot (the same adjacency used for the neighbourhood statistics),
- the permutation-tested neighbour enrichment per gene (spatial_neighbors.parquet), so the viewer can show next to the
  live view which neighbour cell types are significantly enriched across all 23 sections.
The browser only recomputes descriptive quantities (top-decile spots of the picked gene, mean neighbour composition);
the statistics shown as p-values always come from the Python analysis.
"""
from __future__ import annotations

import base64
import io
import json

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from PIL import Image

from .config import analysis_targets, load
from .paths import PROCESSED, RAW
from .spatial import SPATIAL_GROUPS, neighbor_matrix, reference

SC = RAW / "singlecell"
SECTIONS = {"6800STDY12499408": "core", "6800STDY12499410": "core", "6800STDY12499406": "interface",
            "6800STDY12499508": "interface"}
MAX_PX = 1300
JPEG_QUALITY = 72
EXPR_SCALE = 20.0         # stored value = round(log1p(CPM) * EXPR_SCALE)


def viewer_genes() -> list[str]:
    cfg = load()["cell_context"]["spatial"]
    nb = pd.read_parquet(PROCESSED / "spatial_neighbors.parquet")
    genes = list(dict.fromkeys(cfg["map_genes"] + analysis_targets()))
    return [g for g in genes if g in set(nb["gene"])]


def section_payload(sample: str, region: str, dec: pd.DataFrame, genes: list[str], sym2ens: dict) -> dict | None:
    path = SC / f"li2022_visium_section_{sample}.h5ad"
    if not path.exists():
        return None
    a = ad.read_h5ad(path)
    key = [k for k in a.uns["spatial"] if k != "is_single"][0]
    lib = a.uns["spatial"][key]
    scale = float(lib["scalefactors"]["tissue_hires_scalef"])
    diam = float(lib["scalefactors"]["spot_diameter_fullres"]) * scale
    img = np.asarray(lib["images"]["hires"])
    ids = [f"{sample}_{b}" for b in a.obs_names]
    sec = dec[dec["sample"] == sample].set_index("spot")
    keep = np.array([i in sec.index for i in ids])
    a = a[keep].copy()
    ids = [i for i, k in zip(ids, keep) if k]
    sec = sec.loc[ids]
    xy = a.obsm["spatial"].astype(float) * scale
    pad = diam * 2
    x0, y0 = max(int(xy[:, 0].min() - pad), 0), max(int(xy[:, 1].min() - pad), 0)
    x1, y1 = min(int(xy[:, 0].max() + pad), img.shape[1]), min(int(xy[:, 1].max() + pad), img.shape[0])
    crop = Image.fromarray(img[y0:y1, x0:x1])
    f = min(1.0, MAX_PX / max(crop.size))
    crop = crop.resize((max(int(crop.size[0] * f), 1), max(int(crop.size[1] * f), 1)), Image.LANCZOS)
    buf = io.BytesIO()
    crop.convert("RGB").save(buf, format="JPEG", quality=JPEG_QUALITY, optimize=True)
    pos = np.column_stack([(xy[:, 0] - x0) * f, (xy[:, 1] - y0) * f])
    X = sp.csr_matrix(a.X, dtype=np.float64)
    libsize = np.asarray(X.sum(axis=1)).ravel()
    var = pd.Index(a.var_names)
    expr = {}
    for g in genes:
        ens = sym2ens.get(g)
        if ens in var:
            cpm = np.asarray(X[:, var.get_loc(ens)].todense()).ravel() / np.maximum(libsize, 1) * 1e6
            expr[g] = np.round(np.log1p(cpm) * EXPR_SCALE).astype(int).tolist()
    # ring-1 neighbours with the same rule as the statistics (computed on this section's spots only)
    spots = pd.DataFrame({"x": xy[:, 0], "y": xy[:, 1], "section": sample})
    A = neighbor_matrix(spots).tocsr()
    nbrs = [A.indices[A.indptr[i]:A.indptr[i + 1]].tolist() for i in range(len(spots))]
    groups = list(SPATIAL_GROUPS)
    return {"id": sample, "region": region, "donor": str(a.obs["donor_id"].iloc[0]), "n_spots": len(ids),
            "width": crop.size[0], "height": crop.size[1], "spot_px": round(diam * f, 2),
            "image": "data:image/jpeg;base64," + base64.b64encode(buf.getvalue()).decode(),
            "x": np.round(pos[:, 0], 1).tolist(), "y": np.round(pos[:, 1], 1).tolist(),
            "shares": {g: np.round(sec[g].to_numpy() * 100).astype(int).tolist() for g in groups},
            "expr": expr, "neighbors": nbrs}


def run() -> dict:
    dec = pd.read_parquet(PROCESSED / "spatial_deconvolution.parquet")
    genes = viewer_genes()
    ref = reference()
    sym2ens = ref["hgnc_symbol"].dropna().reset_index().drop_duplicates("hgnc_symbol").set_index("hgnc_symbol").iloc[:, 0].to_dict()
    sections = [p for s, r in SECTIONS.items() if (p := section_payload(s, r, dec, genes, sym2ens))]
    nb = pd.read_parquet(PROCESSED / "spatial_neighbors.parquet")
    nb = nb[nb["gene"].isin(genes)]
    stats = {g: [{"cell_group": r.cell_group, "log2": round(float(r.log2_enrichment), 3), "padj": float(f"{r.padj:.3g}"),
                  "observed": round(float(r.observed_neighbor_share), 4), "expected": round(float(r.expected_neighbor_share), 4)}
              for r in d.itertuples()] for g, d in nb.groupby("gene")}
    n_high = nb.groupby("gene")[["n_high_spots", "n_sections"]].first().to_dict("index")
    out = {"groups": list(SPATIAL_GROUPS), "genes": [g for g in genes if any(g in s["expr"] for s in sections)],
           "expr_scale": EXPR_SCALE, "sections": sections, "stats": stats, "stats_n": n_high,
           "validation": json.loads((PROCESSED / "spatial_summary.json").read_text())["validation"]}
    p = PROCESSED / "spatial_viewer.json"
    p.write_text(json.dumps(out, separators=(",", ":")), encoding="utf-8")
    return {"sections": [s["id"] for s in sections], "spots": sum(s["n_spots"] for s in sections),
            "genes": len(out["genes"]), "megabytes": round(p.stat().st_size / 1e6, 1)}


if __name__ == "__main__":
    print(json.dumps(run(), indent=1))
