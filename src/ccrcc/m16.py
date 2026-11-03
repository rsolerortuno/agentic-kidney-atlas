"""M16: evidence layers added after the 2026-11-03 review (deterministic, no LLM).

Each source fills a gap the report left open:

- Surfaceome. The discovery screen called any "membrane" location antibody-accessible, so Golgi, ER and mitochondrial
  proteins scored as surface targets. Two curated surfaceome resources now decide, next to UniProt:
  SURFY (Bausch-Fluck et al. 2018 PNAS, in-silico human surfaceome, table S3) and the Cell Surface Protein Atlas (CSPA,
  Bausch-Fluck et al. 2015 PLoS One, S2 table A: proteins observed at the cell surface by mass spectrometry).
- Normal-tissue proteome. Safety rested on HPA RNA alone. The GTEx proteome (Jiang et al. 2020 Cell, table S2:
  32 normal tissues, 12,627 proteins) adds protein abundance and tissue-specificity (TS) scores in vital organs.
- CheckMate 025 (Braun et al. 2020 Nat Med, supplementary tables S1 and S4A): pre-treatment RNA of 250 patients randomised
  to nivolumab or everolimus, a second immunotherapy-versus-other-drug contrast next to IMmotion150.
- PRISM Repurposing 20Q2 primary screen (Corsello et al. 2020 Nat Cancer; DepMap figshare): drug response of the same cell
  lines DepMap uses, so CRISPR dependencies can be compared with what compounds against a gene's product do.

Files go to data/raw/m16/ with their own manifest (URL, UTC date, sha256, bytes). Nothing here is redistributed.
"""
from __future__ import annotations

import hashlib
import json
import re
from datetime import datetime, timezone

import numpy as np
import pandas as pd
from scipy import stats

from .paths import PROCESSED, RAW

M16 = RAW / "m16"
MANIFEST = M16 / "manifest.json"
PRISM = "https://ndownloader.figshare.com/files/"
FILES = {
    "surfy_table_S3.xlsx": "https://wollscheidlab.org/SURFY/table_S3_surfaceome.xlsx",
    "cspa_S2_File.xlsx": "https://wollscheidlab.org/CSPA/data/S2_File.xlsx",
    "gtex_proteome_table_S2.xlsx": "https://ars.els-cdn.com/content/image/1-s2.0-S0092867420310783-mmc3.xlsx",
    "braun2020_supplementary_tables.xlsx": ("https://static-content.springer.com/esm/art%3A10.1038%2Fs41591-020-0839-y/"
                                            "MediaObjects/41591_2020_839_MOESM2_ESM.xlsx"),
    "prism_20q2_primary_lfc.csv": PRISM + "36794580",
    "prism_20q2_primary_treatment_info.csv": PRISM + "36794583",
    "prism_20q2_primary_cell_line_info.csv": PRISM + "36794565",
}
# vital organs, as in the HPA-based safety score (kidney is the tumour's own organ and is not in the GTEx proteome)
VITAL_PROTEOME = ["Adrenal Gland", "Brain Cerebellum", "Brain Cortex", "Colon Sigmoid", "Colon Transverse",
                  "Esophagus Mucosa", "Esophagus Muscle", "Heart Atrial", "Heart Ventricle", "Liver", "Lung",
                  "Muscle Skeletal", "Pancreas", "Small Intestine", "Spleen", "Stomach"]
TS_ENRICHED = 4.0          # Jiang et al. 2020: TS >= 4 = tissue-enriched protein
MIN_CCRCC_LINES = 3


def _sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def download() -> dict:
    from .extra_download import _get
    M16.mkdir(parents=True, exist_ok=True)
    prev = {e["file"]: e for e in (json.loads(MANIFEST.read_text())["files"] if MANIFEST.exists() else [])}
    entries = []
    for name, url in FILES.items():
        dest = M16 / name
        # the publisher CDNs reject requests without a browser user agent
        _get(url, dest, headers={"User-Agent": "Mozilla/5.0 (ccrcc-graphrag data download)"})
        p = prev.get(name)
        same = p and p.get("bytes") == dest.stat().st_size
        entries.append({"file": name, "url": url, "bytes": dest.stat().st_size,
                        "downloaded_utc": p["downloaded_utc"] if same else datetime.now(timezone.utc).isoformat(timespec="seconds"),
                        "sha256": p["sha256"] if same else _sha256(dest)})
    MANIFEST.write_text(json.dumps({"source": "M16 layers (SURFY, CSPA, GTEx proteome, CheckMate 025, PRISM)",
                                    "files": entries}, indent=1))
    return {"files": len(entries), "bytes": sum(e["bytes"] for e in entries)}


# ---------------------------------------------------------------- surfaceome
CELL_SURFACE = re.compile(r"cell membrane|cell surface|plasma membrane|cell projection[^;]*membrane|microvillus membrane|"
                          r"cilium membrane|synaptic[^;]*membrane|postsynaptic|presynaptic", re.I)
SECRETED = re.compile(r"\bsecreted\b|extracellular", re.I)
GENERIC = re.compile(r"^\s*membrane\s*$", re.I)


def surfaceome() -> pd.DataFrame:
    """One row per UniProt accession: SURFY label and CSPA category."""
    from .xlsx import read_sheet
    s = read_sheet(M16 / "surfy_table_S3.xlsx", "SurfaceomeMasterTable", header=1)
    s = s[["UniProt accession", "UniProt gene", "Surfaceome Label", "Surfaceome Label Source"]].rename(columns={
        "UniProt accession": "accession", "UniProt gene": "surfy_gene", "Surfaceome Label": "surfy_label",
        "Surfaceome Label Source": "surfy_source"})
    c = read_sheet(M16 / "cspa_S2_File.xlsx", "Table A", header=0)
    c = c[c["organism"].eq("Human")][["ID_link", "CSPA category"]].rename(columns={"ID_link": "accession",
                                                                                 "CSPA category": "cspa_category"})
    out = s.merge(c, on="accession", how="outer").dropna(subset=["accession"]).drop_duplicates("accession")
    out = _text(out)
    out.to_parquet(PROCESSED / "surfaceome.parquet", index=False)
    return out


def _text(df: pd.DataFrame) -> pd.DataFrame:
    """Spreadsheet text columns mix numbers and strings (gene names such as '1-Mar'); store them as text."""
    for col in df.columns:
        if df[col].dtype == object:
            df[col] = df[col].map(lambda v: v if v is None or (isinstance(v, float) and np.isnan(v)) else str(v))
    return df


def surface_call(locations, accession: str | None, table: pd.DataFrame | None) -> tuple[str, str]:
    """'surface', 'secreted', 'internal' or 'unknown', with the evidence that decided it.

    Order: CSPA high confidence (observed at the surface) > a specific UniProt location > SURFY (in-silico), which only
    decides when UniProt says nothing more specific than "Membrane". A curated organelle location beats a prediction.
    """
    locs = [x for x in (locations or []) if isinstance(x, str) and x.strip()]
    row = table.loc[accession] if table is not None and accession in getattr(table, "index", []) else None
    surfy = None if row is None or pd.isna(row.get("surfy_label")) else str(row["surfy_label"])
    cspa = None if row is None or pd.isna(row.get("cspa_category")) else str(row["cspa_category"])
    src = f"UniProt: {'; '.join(locs) or 'none'} | SURFY: {surfy or 'not classified'} | CSPA: {cspa or 'not observed'}"
    if cspa and cspa.startswith("1"):
        return "surface", src
    if any(CELL_SURFACE.search(x) for x in locs):
        return "surface", src
    if any(SECRETED.search(x) for x in locs):
        return "secreted", src
    generic_only = not locs or any(GENERIC.match(x) for x in locs)
    if not generic_only:
        return "internal", src
    if surfy == "surface":
        return "surface", src
    if surfy == "nonsurface":
        return "internal", src
    return "unknown", src


SURFACE_VALUE = {"surface": 1.0, "secreted": 1.0, "unknown": 0.5, "internal": 0.0}


def load_surfaceome() -> pd.DataFrame | None:
    p = PROCESSED / "surfaceome.parquet"
    return pd.read_parquet(p).set_index("accession") if p.exists() else None


# ---------------------------------------------------------------- normal-tissue proteome
def _num(x) -> float:
    """TS cells read '0.02;NA_one_rep_in_raw' when the protein was seen once in that tissue: not reliable, NaN."""
    if isinstance(x, str):
        return np.nan if ";" in x or x.strip() in ("", "NA") else float(x)
    return np.nan if x is None else float(x)


def _argmax(df: pd.DataFrame) -> pd.Series:
    """Column of the row maximum, None where the row has no value."""
    return df.fillna(-np.inf).idxmax(axis=1).where(df.notna().any(axis=1))


def normal_proteome() -> pd.DataFrame:
    from .omics import to_hgnc
    from .xlsx import read_sheet
    path = M16 / "gtex_proteome_table_S2.xlsx"
    ts = read_sheet(path, "G protein TS score", header=1)
    med = read_sheet(path, "E protein tissue median", header=1).rename(columns={"gene.id": "ensembl_id"})
    vit_ts = ts[VITAL_PROTEOME].map(_num)
    vit_med = med.set_index("ensembl_id")[VITAL_PROTEOME].map(_num).reindex(ts["ensembl_id"]).reset_index(drop=True)
    out = pd.DataFrame({
        "gene": ts["hgnc_symbol"], "ensembl_id": ts["ensembl_id"],
        "vital_protein_detected": vit_med.notna().sum(axis=1).astype(int),
        "vital_protein_max_rel": vit_med.max(axis=1),
        "vital_protein_max_tissue": _argmax(vit_med),
        "vital_protein_ts_max": vit_ts.max(axis=1),
        "vital_protein_ts_tissue": _argmax(vit_ts),
    })
    out["vital_protein_enriched"] = out["vital_protein_ts_max"] >= TS_ENRICHED
    out = to_hgnc(_text(out.dropna(subset=["gene"])), "gene", rank_by="vital_protein_detected")
    out.to_parquet(PROCESSED / "normal_proteome.parquet", index=False)
    return out


# ---------------------------------------------------------------- PRISM drug response
def prism() -> pd.DataFrame:
    """Per compound: ccRCC lines vs all other lines (log2 fold change of viability, lower = more killing). Per target
    gene: the compound that is most selective for ccRCC lines. Selectivity = mean LFC other - mean LFC ccRCC."""
    from .layers import bh, ccrcc_lines
    from .omics import to_hgnc
    lfc = pd.read_csv(M16 / "prism_20q2_primary_lfc.csv", index_col=0)
    info = pd.read_csv(M16 / "prism_20q2_primary_treatment_info.csv")
    rcc = [a for a in ccrcc_lines()["ach"] if a in lfc.index]
    other = lfc.index.difference(rcc)
    R, O = lfc.loc[rcc], lfc.loc[other]
    nr = R.notna().sum()
    with np.errstate(invalid="ignore"):
        res = stats.ttest_ind(R, O, axis=0, equal_var=False, nan_policy="omit")
    comp = pd.DataFrame({"column_name": lfc.columns, "n_ccrcc_lines": nr.to_numpy(), "n_other_lines": O.notna().sum().to_numpy(),
                         "mean_lfc_ccrcc": R.mean().to_numpy(), "mean_lfc_other": O.mean().to_numpy(),
                         "pvalue": np.asarray(res.pvalue, float)})
    comp = comp.merge(_text(info[["column_name", "broad_id", "name", "moa", "target", "phase", "dose"]].copy()),
                      on="column_name", how="left")
    comp = comp[comp["n_ccrcc_lines"] >= MIN_CCRCC_LINES].copy()
    comp["selectivity"] = comp["mean_lfc_other"] - comp["mean_lfc_ccrcc"]
    comp["padj"] = bh(comp["pvalue"])
    comp.to_parquet(PROCESSED / "prism_compounds.parquet", index=False)
    rows = []
    for r in comp.dropna(subset=["target"]).itertuples():
        targets = [g.strip() for g in str(r.target).split(",") if g.strip()]
        for g in targets:
            # compound_targets: a promiscuous compound (baicalein lists 7 targets) says little about any one of them
            rows.append({"gene": g, "compound": r.name, "moa": r.moa, "phase": r.phase, "compound_targets": len(targets),
                         "selectivity": r.selectivity, "mean_lfc_ccrcc": r.mean_lfc_ccrcc, "mean_lfc_other": r.mean_lfc_other,
                         "pvalue": r.pvalue, "padj": r.padj, "n_ccrcc_lines": int(r.n_ccrcc_lines)})
    per = pd.DataFrame(rows)
    n = per.groupby("gene").size().rename("n_compounds")
    # prefer the most selective compound among those with at most 3 annotated targets; fall back to any compound
    per["specific"] = per["compound_targets"] <= 3
    best = (per.sort_values(["specific", "selectivity"], ascending=False).drop_duplicates("gene").set_index("gene")
            .join(n).reset_index())
    best = to_hgnc(best, "gene", rank_by="selectivity")
    best.to_parquet(PROCESSED / "layer_prism.parquet", index=False)
    return best


# ---------------------------------------------------------------- CheckMate 025
def checkmate() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Clinical table (index RNA id) and expression (genes x RNA ids) for the randomised CheckMate 025 patients.

    PFS_CNSR is read as the event indicator (1 = progression or death): 276 of 311 patients with RNA have 1, as expected
    for advanced RCC, where most patients progress during follow-up."""
    cache_c, cache_x = PROCESSED / "checkmate_clinical.parquet", PROCESSED / "checkmate_expression.parquet"
    if not (cache_c.exists() and cache_x.exists()):
        from .xlsx import read_sheet
        path = M16 / "braun2020_supplementary_tables.xlsx"
        c = read_sheet(path, "S1_Clinical_and_Immune_Data", header=1)
        x = read_sheet(path, "S4A_RNA_Expression", header=1)
        x = x.dropna(subset=["gene_name"]).assign(gene_name=lambda d: d["gene_name"].astype(str))
        x = x.groupby("gene_name").mean(numeric_only=True)
        _text(c).to_parquet(cache_c, index=False)
        x.to_parquet(cache_x)
    c, x = pd.read_parquet(cache_c), pd.read_parquet(cache_x)
    c = c[c["Cohort"].eq("CM-025") & c["RNA_ID"].isin(x.columns)].drop_duplicates("RNA_ID").set_index("RNA_ID")
    c["arm"] = c["Arm"].map({"NIVOLUMAB": "nivolumab_cm025", "EVEROLIMUS": "everolimus_cm025"})
    c["pfs_time"] = pd.to_numeric(c["PFS"], errors="coerce")
    c["pfs_event"] = pd.to_numeric(c["PFS_CNSR"], errors="coerce").fillna(0).astype(int)
    c["responder"] = c["ORR"].map({"CR": 1.0, "PR": 1.0, "CRPR": 1.0, "SD": 0.0, "PD": 0.0})
    return c.dropna(subset=["arm"]), x[c.index]


def run() -> dict:
    sur = surfaceome()
    prot = normal_proteome()
    pr = prism()
    return {"surfaceome": {"accessions": int(len(sur)), "surfy_surface": int((sur.surfy_label == "surface").sum()),
                           "cspa_high_confidence": int(sur.cspa_category.astype(str).str.startswith("1").sum())},
            "normal_proteome": {"genes": int(len(prot)), "vital_enriched": int(prot.vital_protein_enriched.sum())},
            "prism": {"target_genes": int(len(pr)), "selective_padj_0.1": int((pr.padj < 0.1).sum())}}


if __name__ == "__main__":
    print(json.dumps(run(), indent=1))
