"""M13 downloads: the extra evidence layers added after the first evaluation.

Every file goes to data/raw/extra/ and is recorded in data/raw/extra/manifest.json (URL, UTC date, sha256, bytes).
Sources, all public and redistributable only through the original providers (nothing here is committed):

- Somatic mutations: cBioPortal REST API, CPTAC renal (rcc_cptac_gdc, same C3L/C3N patient ids as LinkedOmics)
  and TCGA-KIRC PanCancer Atlas (kirc_tcga_pan_can_atlas_2018).
- GENCODE v22 GTF (via the cptac Zenodo record) for coding-sequence length, used to correct mutation recurrence
  for gene length.
- DepMap Public 26Q1 Chronos gene effect (Figshare) + Cellosaurus for renal cell carcinoma line annotation.
- TCGA-KIRC RNA-seq (UCSC Xena, log2(RSEM+1), tumor + adjacent normal) and overall survival (Xena curated table).
- CPTAC ccRCC phosphosite abundance (BCM pipeline) and ESTIMATE tumour purity (WashU), both from the cptac Zenodo
  record 8394329.
- OmniPath enzyme-substrate relations (kinase -> phosphosite).
- Open Targets Platform 26.09 bulk tables: target tractability (all targets), clinical targets and clinical
  precedence (trial start dates, used for the temporal holdout).
"""
from __future__ import annotations

import gzip
import hashlib
import json
import shutil
import time
from datetime import datetime, timezone

import pandas as pd
import requests

from .paths import RAW

EXTRA = RAW / "extra"
MANIFEST = EXTRA / "manifest.json"
ZENODO = "https://zenodo.org/api/records/8394329/files/{}/content"
OT = "https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/26.09/output/"
CBIO = "https://www.cbioportal.org/api"

FILES = {
    "gencode.v22.annotation.gtf.gz": ZENODO.format("washu-ccrcc-mapping-gencode.v22.annotation.gtf.gz"),
    "depmap_26Q1_gene_effect.csv": "https://ndownloader.figshare.com/files/67214582",
    "tcga_kirc_hiseqv2.gz": "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.KIRC.sampleMap%2FHiSeqV2.gz",
    "tcga_kirc_survival.txt": "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/survival%2FKIRC_survival.txt",
    "cptac_ccrcc_phosphosite_tumor.txt": ZENODO.format(
        "bcm-ccrcc-phosphoproteomics-CCRCC_phospho_site_abundance_log2_reference_intensity_normalized_Tumor.txt"),
    "cptac_ccrcc_phosphosite_normal.txt.gz": ZENODO.format(
        "bcm-ccrcc-phosphoproteomics-CCRCC_phospho_site_abundance_log2_reference_intensity_normalized_Normal.txt.gz"),
    "cptac_pancan_purity_estimate.tsv.gz": ZENODO.format(
        "washu-all_cancers-tumor_purity-CPTAC_pancan_RNA_tumor_purity_ESTIMATE_WashU.tsv.gz"),
    "omnipath_enzsub.tsv": "https://omnipathdb.org/enzsub?format=tsv&genesymbols=yes&types=phosphorylation",
    "ot_clinical_target.parquet": OT + "clinical_target/00000000.parquet",
}
OT_MULTIPART = {"ot_target_tractability": "target_tractability", "ot_clinical_precedence": "evidence_clinical_precedence"}
CBIO_STUDIES = {"cptac": ("rcc_cptac_gdc_mutations", "rcc_cptac_gdc_all"),
                "tcga": ("kirc_tcga_pan_can_atlas_2018_mutations", "kirc_tcga_pan_can_atlas_2018_sequenced")}
CELLOSAURUS = ("https://api.cellosaurus.org/search/cell-line?q=di:%22renal%20cell%20carcinoma%22"
               "&fields=id,ac,dr,di&rows=2000&format=txt")


def _sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _get(url: str, dest, **kw) -> None:
    if dest.exists() and dest.stat().st_size > 0:
        return
    tmp = dest.with_suffix(dest.suffix + ".part")
    for attempt in range(5):
        try:
            with requests.get(url, stream=True, timeout=300, **kw) as r:
                r.raise_for_status()
                with open(tmp, "wb") as f:
                    for chunk in r.iter_content(1 << 20):
                        f.write(chunk)
            tmp.replace(dest)
            return
        except (requests.ConnectionError, requests.Timeout, requests.HTTPError):
            if attempt == 4:
                raise
            time.sleep(2 ** attempt)


def _ot_parts(folder: str) -> list[str]:
    import re
    html = requests.get(OT + folder + "/", timeout=60).text
    return sorted(set(re.findall(r'href="(part-[^"]+\.parquet)"', html)))


def _cbio_mutations(profile: str, sample_list: str, dest) -> None:
    if dest.exists() and dest.stat().st_size > 0:
        return
    r = requests.post(f"{CBIO}/molecular-profiles/{profile}/mutations/fetch?projection=DETAILED",
                      json={"sampleListId": sample_list}, timeout=600)
    r.raise_for_status()
    rows = [{"sample": m["sampleId"], "patient": m["patientId"], "gene": (m.get("gene") or {}).get("hugoGeneSymbol"),
             "entrez": m.get("entrezGeneId"), "mutation_type": m.get("mutationType"),
             "protein_change": m.get("proteinChange"), "variant_type": m.get("variantType")} for m in r.json()]
    pd.DataFrame(rows).to_parquet(dest, index=False)
    samples = requests.get(f"{CBIO}/sample-lists/{sample_list}", timeout=120).json()["sampleIds"]
    dest.with_suffix(".samples.json").write_text(json.dumps(samples))


def download_extra() -> dict:
    EXTRA.mkdir(parents=True, exist_ok=True)
    previous = json.loads(MANIFEST.read_text()) if MANIFEST.exists() else {}
    prev = {e["file"]: e for e in previous.get("files", [])}
    targets: list[tuple[str, str]] = list(FILES.items())
    for name, folder in OT_MULTIPART.items():
        for i, part in enumerate(_ot_parts(folder)):
            targets.append((f"{name}_{i}.parquet", OT + folder + "/" + part))
    for name, url in targets:
        if (EXTRA / (name + ".gz")).exists():
            continue
        _get(url, EXTRA / name)
    for key, (profile, slist) in CBIO_STUDIES.items():
        _cbio_mutations(profile, slist, EXTRA / f"mutations_{key}.parquet")
    _get(CELLOSAURUS, EXTRA / "cellosaurus_rcc.txt")
    # the BCM tumour phosphosite table is published uncompressed (~0.3 GB); keep a gzip copy only
    raw = EXTRA / "cptac_ccrcc_phosphosite_tumor.txt"
    gz = EXTRA / "cptac_ccrcc_phosphosite_tumor.txt.gz"
    if raw.exists() and not gz.exists():
        with open(raw, "rb") as fi, gzip.open(gz, "wb", compresslevel=5) as fo:
            shutil.copyfileobj(fi, fo)
        raw.unlink()
    entries = []
    for f in sorted(EXTRA.iterdir()):
        if f.name == "manifest.json" or f.suffix == ".part" or f.stat().st_size == 0:
            continue
        url = dict(targets).get(f.name) or ("cBioPortal API" if f.name.startswith("mutations_") else
                                            CELLOSAURUS if f.name.startswith("cellosaurus") else
                                            dict(targets).get(f.name.removesuffix(".gz"), "derived"))
        p = prev.get(f.name)
        same = p and p.get("bytes") == f.stat().st_size
        entries.append({"file": f.name, "url": url,
                        "downloaded_utc": p["downloaded_utc"] if same else
                        datetime.now(timezone.utc).isoformat(timespec="seconds"),
                        "sha256": p["sha256"] if same else _sha256(f), "bytes": f.stat().st_size})
    manifest = {"source": "M13 extra layers (cBioPortal, DepMap, UCSC Xena, cptac Zenodo, OmniPath, Open Targets)",
                "files": entries}
    MANIFEST.write_text(json.dumps(manifest, indent=1))
    return {"files": len(entries), "bytes": sum(e["bytes"] for e in entries)}


if __name__ == "__main__":
    print(json.dumps(download_extra(), indent=1))
