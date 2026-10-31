"""Two annotation layers for target discovery (deterministic):

- normal-tissue safety: Human Protein Atlas consensus RNA (nTPM per tissue). A target expressed in vital organs risks on-target,
  off-tumor toxicity; kidney expression is not penalized (the tumor arises there, and the kidney is the tissue of interest).
- clinical association: Spearman correlation of tumor RNA / protein with histologic grade and pathologic stage across the CPTAC
  discovery tumors (the clinical file was downloaded in M1 but unused until now), BH-adjusted per layer and variable.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

from .config import load
from .omics import LAYERS, load_hgnc, read_matrix
from .paths import PROCESSED, RAW

HPA = RAW / "hpa" / "rna_tissue_consensus.tsv"
HPA_URL = "https://www.proteinatlas.org/download/tsv/rna_tissue_consensus.tsv.zip"
VITAL = ["heart muscle", "liver", "lung", "cerebral cortex", "bone marrow", "pancreas", "small intestine", "colon"]
GRADE = {"G1": 1, "G2": 2, "G3": 3, "G4": 4}
STAGE = {"Stage I": 1, "Stage II": 2, "Stage III": 3, "Stage IV": 4}


def download_hpa(force: bool = False) -> None:
    import io
    import zipfile

    import requests
    if HPA.exists() and not force:
        return
    import time
    HPA.parent.mkdir(parents=True, exist_ok=True)
    for attempt in range(6):   # proteinatlas.org occasionally resets the connection
        try:
            r = requests.get(HPA_URL, timeout=300)
            r.raise_for_status()
            break
        except (requests.ConnectionError, requests.Timeout, requests.HTTPError):
            if attempt == 5:
                raise
            time.sleep(min(2 ** attempt, 30))
    zipfile.ZipFile(io.BytesIO(r.content)).extract("rna_tissue_consensus.tsv", HPA.parent)


def normal_tissue() -> pd.DataFrame:
    download_hpa()
    t = pd.read_csv(HPA, sep="\t")
    wide = t.pivot_table(index="Gene", columns="Tissue", values="nTPM", aggfunc="max")
    tau = wide.apply(lambda r: float((1 - r / r.max()).sum() / (len(r) - 1)) if r.max() > 0 else np.nan, axis=1)
    vital = wide[[c for c in VITAL if c in wide.columns]]
    out = pd.DataFrame({
        "vital_max_ntpm": vital.max(axis=1), "vital_max_tissue": vital.idxmax(axis=1),
        "kidney_ntpm": wide["kidney"] if "kidney" in wide.columns else np.nan,
        "tissue_tau": tau, "n_tissues_expressed": (wide >= 1).sum(axis=1), "max_ntpm": wide.max(axis=1)})
    out.index.name = "ensembl"
    return out.reset_index()


def clinical_associations() -> pd.DataFrame:
    cli = pd.read_csv(RAW / "HS_CPTAC_CCRCC_CLI.tsi", sep="\t", index_col=0).iloc[1:]
    grade, stage = cli["Histologic_Grade"].map(GRADE), cli["Tumor_Stage_Pathological"].map(STAGE)
    hg = load_hgnc().drop_duplicates("symbol").set_index("symbol")
    rows = []
    for layer in ("rna", "protein"):
        m = read_matrix(LAYERS[layer][0])
        for var, y in (("grade", grade), ("stage", stage)):
            common = [c for c in m.columns if c in y.index and pd.notna(y[c])]
            yy = y[common].to_numpy(float)
            X = m[common]
            ok = X.notna().sum(axis=1) >= 40
            X = X[ok]
            rho, p = [], []
            for g, v in X.iterrows():
                mask = v.notna().to_numpy()
                r, pv = stats.spearmanr(v.to_numpy()[mask], yy[mask])
                rho.append(r)
                p.append(pv)
            p = np.asarray(p, float)
            good = np.isfinite(p)
            padj = np.full(len(p), np.nan)
            padj[good] = multipletests(p[good], method="fdr_bh")[1]
            rows.append(pd.DataFrame({"dataset_symbol": X.index, "layer": layer, "variable": var, "rho": rho,
                                      "pvalue": p, "padj": padj, "n": len(common)}))
    d = pd.concat(rows, ignore_index=True)
    d["hgnc_symbol"] = d["dataset_symbol"].map(lambda s: s if s in hg.index else None)
    return d.dropna(subset=["hgnc_symbol"])


def run() -> dict:
    ntis = normal_tissue()
    ntis.to_parquet(PROCESSED / "normal_tissue.parquet", index=False)
    clin = clinical_associations()
    clin.to_parquet(PROCESSED / "clinical_associations.parquet", index=False)
    sig = clin[clin.padj < 0.05]
    out = {"hpa_genes": int(len(ntis)), "clinical_tests": int(len(clin)),
           "significant_by_layer_variable": sig.groupby(["layer", "variable"]).size().to_dict().__repr__(),
           "tumors_with_clinical": int(clin["n"].max())}
    (PROCESSED / "safety_clinical_summary.json").write_text(json.dumps(out, indent=1))
    return out


if __name__ == "__main__":
    print(json.dumps(run(), indent=1))
