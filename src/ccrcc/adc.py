"""Suitability for antibody-drug conjugates (ADC) and radioligand therapy (RLT).

Deterministic, no LLM. An ADC or a radioligand needs more than over-expression: the protein must sit on the cell
surface with an extracellular part to bind, be absent (as protein) from vital normal tissues, and be present on most
tumour cells, not on a subset, or the payload misses the rest. Each candidate gets four components in [0, 1]:

- surface:      UniProt topology and location. 1.0 = transmembrane with an extracellular topological domain of
                >= 20 residues and a cell-membrane / cell-surface location; 0.7 = transmembrane + membrane location
                without an annotated extracellular domain; 0.4 = membrane location only (GPI-anchored or peripheral);
                0 = intracellular or secreted. A protein observed at the cell surface by mass spectrometry (Cell Surface
                Protein Atlas, high confidence) gets at least 0.7 (M16).
- tumour_up:    min(1, protein log2FC / 2) when the CPTAC protein change is significant, else 0.5 x the RNA version.
- normal_ihc:   Human Protein Atlas normal-tissue immunohistochemistry. 1 - (max staining level over vital tissues) / 3,
                levels Not detected 0, Low 1, Medium 2, High 3; antibodies rated "Uncertain" count half. Normal kidney
                is reported separately (radioligands are cleared by the kidney) but is not part of the score. A protein
                enriched in a vital organ in the GTEx proteome (tissue-specificity score >= 4) is capped at 0.33 (M16).
- homogeneity:  share of ccRCC tumour cells with >= 1 UMI in scRNA-seq (Li 2022), taken as the 25th percentile over
                donors, so a gene expressed in every cell of a few tumours does not look homogeneous.

suitability = 0.30 surface + 0.25 tumour_up + 0.25 normal_ihc + 0.20 homogeneity (weights fixed before looking at the
ranking). Known ccRCC ADC / RLT targets serve as positive controls: CA9 (girentuximab, 177Lu/89Zr-girentuximab),
CD70 (several ADCs), ENPP3 (AGS-16C3F), HAVCR1 / TIM-1 (CDX-014) and FOLH1 / PSMA (tumour neovasculature).
"""
from __future__ import annotations

import json
import re
import zipfile

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

from .external import UNIPROT_URL, _cached, _request
from .paths import PROCESSED, RAW
from .singlecell import SC_FILE, cell_metadata, gene_table

EXTRA = RAW / "extra"
WEIGHTS = {"surface": 0.30, "tumour_up": 0.25, "normal_ihc": 0.25, "homogeneity": 0.20}
CONTROLS = ["CA9", "CD70", "ENPP3", "HAVCR1", "FOLH1"]
VITAL = {"heart muscle", "lung", "liver", "cerebral cortex", "cerebellum", "bone marrow", "pancreas", "stomach",
         "colon", "small intestine", "duodenum", "skeletal muscle", "esophagus", "spleen", "adrenal gland"}
LEVEL = {"Not detected": 0, "Low": 1, "Medium": 2, "High": 3}
SURFACE_LOC = re.compile(r"cell membrane|cell surface|apical cell membrane|basolateral cell membrane", re.I)
MEMBRANE_LOC = re.compile(r"membrane", re.I)


def uniprot_topology(symbol: str) -> dict | None:
    params = {"query": f"gene_exact:{symbol} AND organism_id:9606 AND reviewed:true",
              "fields": "accession,ft_transmem,ft_topo_dom,cc_subcellular_location,ft_lipid", "format": "json", "size": 1}
    data = _cached("uniprot_topology", json.dumps(params, sort_keys=True),
                   lambda: _request("GET", UNIPROT_URL, params=params).json())["data"]
    if not data.get("results"):
        return None
    r = data["results"][0]
    feats = r.get("features", [])
    tm = [f for f in feats if f["type"] == "Transmembrane"]
    extra = [f for f in feats if f["type"] == "Topological domain"
             and "extracellular" in (f.get("description") or "").lower()]
    ext_len = sum(f["location"]["end"]["value"] - f["location"]["start"]["value"] + 1 for f in extra
                  if f["location"]["end"].get("value") and f["location"]["start"].get("value"))
    gpi = any(f["type"] == "Lipidation" and "GPI" in (f.get("description") or "") for f in feats)
    locs = []
    for c in r.get("comments", []):
        if c.get("commentType") == "SUBCELLULAR LOCATION":
            locs += [s["location"]["value"] for s in c.get("subcellularLocations", []) if "location" in s]
    return {"accession": r["primaryAccession"], "n_transmembrane": len(tm), "extracellular_residues": int(ext_len),
            "gpi_anchor": gpi, "locations": sorted(set(locs))}


def surface_score(t: dict | None) -> float:
    if not t:
        return 0.0
    loc = " | ".join(t["locations"])
    surf = bool(SURFACE_LOC.search(loc))
    if t["n_transmembrane"] > 0 and t["extracellular_residues"] >= 20 and surf:
        return 1.0
    if (t["n_transmembrane"] > 0 or t["gpi_anchor"]) and surf:
        return 0.7
    if MEMBRANE_LOC.search(loc) and (surf or t["gpi_anchor"]):
        return 0.4
    return 0.0


def normal_ihc() -> pd.DataFrame:
    with zipfile.ZipFile(EXTRA / "hpa_normal_ihc.tsv.zip") as z:
        d = pd.read_csv(z.open(z.namelist()[0]), sep="\t")
    d["level"] = d["Level"].map(LEVEL)
    d = d.dropna(subset=["level"])
    d["weight"] = np.where(d["Reliability"].eq("Uncertain"), 0.5, 1.0)
    d["eff"] = d["level"] * d["weight"]
    tissue = d["Tissue"].str.lower()
    vit = d[tissue.isin(VITAL)].sort_values("eff", ascending=False).drop_duplicates("Gene name")
    kid = d[tissue.eq("kidney")].groupby("Gene name")["level"].max()
    out = pd.DataFrame({"ihc_vital_max_level": vit.set_index("Gene name")["eff"],
                        "ihc_vital_max_tissue": vit.set_index("Gene name")["Tissue"],
                        "ihc_vital_max_cell": vit.set_index("Gene name")["Cell type"]})
    out["ihc_kidney_max_level"] = kid
    out["ihc_reliability"] = d.drop_duplicates("Gene name").set_index("Gene name")["Reliability"]
    return out


def tumour_cell_homogeneity(genes: list[str]) -> pd.DataFrame:
    """Per donor, the share of ccRCC tumour cells (tumour tissue) with >= 1 UMI of each gene; streamed from the h5ad."""
    meta = cell_metadata()
    meta = meta[(meta["group"] == "RCC") & (meta["tissue"] == "tumor")].reset_index(drop=True)
    with h5py.File(SC_FILE, "r") as h:
        gt = gene_table(h)
        cols = gt.index[gt["hgnc_symbol"].isin(genes)].to_numpy()
        names = gt.loc[cols, "hgnc_symbol"].to_numpy()
        g = h["raw/X"]
        indptr = g["indptr"][:]
        donors = meta["donor"].unique()
        di = {d: i for i, d in enumerate(donors)}
        pos = np.zeros((len(donors), len(cols)))
        n = np.zeros(len(donors))
        rows = meta["row"].to_numpy()
        for start in range(0, len(rows), 20000):
            sel = rows[start:start + 20000]
            lo, hi = indptr[sel[0]], indptr[sel[-1] + 1]
            block = sp.csr_matrix((g["data"][lo:hi], g["indices"][lo:hi], indptr[sel[0]:sel[-1] + 2] - lo),
                                  shape=(sel[-1] - sel[0] + 1, len(gt)))[sel - sel[0]][:, cols]
            d_idx = meta["donor"].to_numpy()[start:start + 20000]
            onehot = sp.csr_matrix((np.ones(len(sel)), ([di[x] for x in d_idx], np.arange(len(sel)))),
                                   shape=(len(donors), len(sel)))
            pos += (onehot @ (block > 0).astype(np.float64)).toarray()
            n += np.asarray(onehot.sum(axis=1)).ravel()
    frac = pos / np.maximum(n, 1)[:, None]
    ok = n >= 50                       # donors with too few tumour cells are not used
    f = frac[ok]
    return pd.DataFrame({"gene": names, "rcc_pct_donor_p25": np.percentile(f, 25, axis=0) * 100,
                         "rcc_pct_donor_median": np.median(f, axis=0) * 100,
                         "rcc_pct_donor_min": f.min(axis=0) * 100, "n_donors": int(ok.sum())}).drop_duplicates("gene")


def candidates() -> list[str]:
    """Tumour-up genes (CPTAC RNA or protein) plus the positive controls."""
    u = pd.read_parquet(PROCESSED / "gene_universe.parquet")
    up = u[((u["rna_padj"] < 0.05) & (u["rna_log2fc"] >= 1)) | ((u["protein_padj"] < 0.05) & (u["protein_log2fc"] >= 0.5))]
    return sorted(set(up["gene"]) | set(CONTROLS))


def run() -> dict:
    genes = candidates()
    u = pd.read_parquet(PROCESSED / "gene_universe.parquet").set_index("gene")
    from .m16 import load_surfaceome, surface_call
    surf_tab = load_surfaceome()
    rows = []
    for g in genes:
        t = uniprot_topology(g)
        acc = (t or {}).get("accession")
        call, why = surface_call((t or {}).get("locations", []), acc, surf_tab)
        cspa = surf_tab is not None and acc in surf_tab.index and str(surf_tab.loc[acc, "cspa_category"]).startswith("1")
        rows.append({"gene": g, "uniprot": acc, "n_transmembrane": (t or {}).get("n_transmembrane"),
                     "extracellular_residues": (t or {}).get("extracellular_residues"),
                     "gpi_anchor": (t or {}).get("gpi_anchor"), "locations": "; ".join((t or {}).get("locations", [])),
                     "surface_call": call, "surface_evidence": why, "cspa_high_confidence": bool(cspa),
                     "c_surface": max(surface_score(t), 0.7 if cspa else 0.0)})
    d = pd.DataFrame(rows).set_index("gene")
    surf = d[d["c_surface"] > 0].index.tolist()
    d = d.join(normal_ihc(), how="left")
    d = d.join(tumour_cell_homogeneity(sorted(set(surf) | set(CONTROLS))).set_index("gene"), how="left")
    plfc, ppadj = u["protein_log2fc"].reindex(d.index), u["protein_padj"].reindex(d.index)
    rlfc, rpadj = u["rna_log2fc"].reindex(d.index), u["rna_padj"].reindex(d.index)
    d["protein_log2fc"], d["rna_log2fc"] = plfc, rlfc
    prot = (plfc / 2).clip(0, 1).where(ppadj < 0.05)
    rna = 0.5 * (rlfc / 4).clip(0, 1).where(rpadj < 0.05)
    d["c_tumour_up"] = prot.fillna(rna).fillna(0.0)
    d["c_normal_ihc"] = (1 - d["ihc_vital_max_level"] / 3).clip(0, 1)
    d["ihc_measured"] = d["c_normal_ihc"].notna()
    d["c_normal_ihc"] = d["c_normal_ihc"].fillna(0.5)         # no HPA antibody: neutral, flagged
    np_path = PROCESSED / "normal_proteome.parquet"
    if np_path.exists():
        npr = pd.read_parquet(np_path).drop_duplicates("gene").set_index("gene")
        for col in ("vital_protein_ts_max", "vital_protein_ts_tissue", "vital_protein_enriched"):
            d[col] = npr[col].reindex(d.index)
        enriched = d["vital_protein_enriched"].astype("boolean").fillna(False).astype(bool)
        d["c_normal_ihc"] = np.where(enriched, np.minimum(d["c_normal_ihc"], 0.33), d["c_normal_ihc"])
    d["c_homogeneity"] = (d["rcc_pct_donor_p25"] / 100 / 0.5).clip(0, 1).fillna(0.0)   # 50 % of cells saturates
    d["suitability"] = sum(w * d[f"c_{k}"] for k, w in WEIGHTS.items())
    d["control"] = d.index.isin(CONTROLS)
    d = d.sort_values("suitability", ascending=False)
    d["rank"] = np.arange(1, len(d) + 1)
    d.reset_index().to_parquet(PROCESSED / "adc_suitability.parquet", index=False)
    surf_df = d[d["c_surface"] > 0]
    summary = {"candidates": int(len(d)), "surface": int((d["c_surface"] > 0).sum()),
               "surface_with_ihc": int(surf_df["ihc_measured"].sum()), "weights": WEIGHTS,
               "controls": {g: {"rank": int(d.loc[g, "rank"]), "suitability": round(float(d.loc[g, "suitability"]), 3),
                                "of": int(len(d))} for g in CONTROLS if g in d.index},
               "control_median_percentile": float(np.median([d.loc[g, "rank"] / len(d) for g in CONTROLS if g in d.index])),
               "top": d.head(15).index.tolist()}
    (PROCESSED / "adc_summary.json").write_text(json.dumps(summary, indent=1, default=float))
    return summary


if __name__ == "__main__":
    print(json.dumps(run(), indent=1, default=float))
