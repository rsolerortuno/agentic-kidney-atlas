"""M1: download CPTAC ccRCC (LinkedOmics) bulk files and the HGNC complete set."""
from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone

import pandas as pd
import requests

from .paths import RAW

LINKEDOMICS_BASE = "https://linkedomics.org/data_download/CPTAC-CCRCC/"
CPTAC_FILES = [
    "HS_CPTAC_CCRCC_CLI.tsi",
    "HS_CPTAC_CCRCC_CNV_gene_Tumor.cct",
    "HS_CPTAC_CCRCC_proteome_Tumor.cct",
    "HS_CPTAC_CCRCC_proteome_Normal.cct",
    "HS_CPTAC_CCRCC_phosphoproteome_gene_Tumor.cct",
    "HS_CPTAC_CCRCC_phosphoproteome_gene_Normal.cct",
    "HS_CPTAC_CCRCC_RNAseq_fpkm_log2_Tumor.cct",
    "HS_CPTAC_CCRCC_RNAseq_fpkm_log2_Normal.cct",
]
HGNC_URL = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
MANIFEST = RAW / "manifest.json"


def _sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _fetch(url: str, dest, force: bool = False) -> None:
    if dest.exists() and dest.stat().st_size > 0 and not force:
        return
    tmp = dest.with_suffix(dest.suffix + ".part")
    with requests.get(url, stream=True, timeout=120) as r:
        r.raise_for_status()
        with open(tmp, "wb") as f:
            for chunk in r.iter_content(1 << 20):
                f.write(chunk)
    tmp.replace(dest)


def _dimensions(path) -> dict:
    df = pd.read_csv(path, sep="\t", index_col=0, low_memory=False)
    return {"rows": int(df.shape[0]), "columns": int(df.shape[1])}


def download_all(force: bool = False) -> dict:
    entries = []
    targets = [(LINKEDOMICS_BASE + name, RAW / name) for name in CPTAC_FILES]
    targets.append((HGNC_URL, RAW / "hgnc_complete_set.txt"))
    previous = json.loads(MANIFEST.read_text()) if MANIFEST.exists() else {}
    prev_by_file = {e["file"]: e for e in previous.get("files", [])}
    for url, dest in targets:
        existed = dest.exists()
        _fetch(url, dest, force=force)
        prev = prev_by_file.get(dest.name)
        entries.append({
            "file": dest.name,
            "url": url,
            "downloaded_utc": prev["downloaded_utc"] if (existed and prev and not force)
            else datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "sha256": _sha256(dest),
            "bytes": dest.stat().st_size,
            **_dimensions(dest),
        })
    manifest = {"source": "LinkedOmics CPTAC-CCRCC discovery cohort + HGNC", "files": entries}
    MANIFEST.write_text(json.dumps(manifest, indent=2))
    return manifest


if __name__ == "__main__":
    for e in download_all()["files"]:
        print(f"{e['file']:55s} {e['rows']:>7} x {e['columns']:<5} {e['sha256'][:12]}")
