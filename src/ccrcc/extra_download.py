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

Added for the treatment, subtype, immune, cell-replication and tissue-viewer analyses (same manifest):
- Treatment response: cBioPortal datahub (GitHub LFS) for IMmotion150 (atezolizumab +/- bevacizumab vs sunitinib,
  McDermott et al. 2018) and the nivolumab biomarker trial of Choueiri et al. 2016, both iAtlas-harmonized.
- Immune deconvolution of the CPTAC tumours: CIBERSORT (absolute) and xCell from the same cptac Zenodo record.
- Human Protein Atlas normal-tissue immunohistochemistry (protein in normal organs, for ADC / radioligand safety).
- Second single-cell dataset: Zhang et al. 2021 PNAS ccRCC and normal kidney (CZ CELLxGENE).
- Four single-section Visium datasets of Li et al. 2022 that carry the H&E image (CZ CELLxGENE), for the tissue viewer.
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
    "cptac_ccrcc_cibersort_abs.txt.gz": ZENODO.format("washu-ccrcc-cibersort-CIBERSORT.Output_Abs_ccRCC.txt.gz"),
    "cptac_ccrcc_xcell.txt.gz": ZENODO.format("washu-ccrcc-xcell-ccRCC_xCell.txt.gz"),
    "hpa_normal_ihc.tsv.zip": "https://www.proteinatlas.org/download/tsv/normal_ihc_data.tsv.zip",
}
DATAHUB = "https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/{}/{}"
TREATMENT_STUDIES = {"immotion150": ("rcc_iatlas_immotion150_2018", "data_mrna_seq_tpm.txt"),
                     "choueiri2016": ("ccrcc_iatlas_choueiri_2016", None)}
CXG = "https://datasets.cellxgene.cziscience.com/{}.h5ad"
SINGLECELL_EXTRA = {   # file -> CZ CELLxGENE dataset asset id
    "zhang2021_ccrcc.h5ad": "f77a9049-7cbe-477d-8d3f-b417bba03147",
    "zhang2021_normal.h5ad": "5f0b9675-aed8-4870-8ad8-7a3f38037486",
    "li2022_visium_section_6800STDY12499406.h5ad": "3cb22696-5867-4001-ba23-3a24da12f4ab",
}
# single-section Visium datasets (with H&E) whose asset id is looked up from the collection at download time
VISIUM_SECTIONS = {"6800STDY12499408": "252438d3-8143-48a9-9ac2-ef5991fb21b8",
                   "6800STDY12499410": "a6046b15-a095-43b0-9fb5-b36899d87fdb",
                   "6800STDY12499508": "104cfa2a-8f7e-4dcc-a0a3-8fd46711c904"}
LI2022 = "https://api.cellxgene.cziscience.com/curation/v1/collections/f7cecffa-00b4-4560-a29a-8ad626b8ee08/datasets/{}"
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


def _treatment_studies() -> None:
    for key, (study, expr) in TREATMENT_STUDIES.items():
        if expr is None:   # pick the non-z-score expression file from the datahub listing
            listing = requests.get(f"https://api.github.com/repos/cBioPortal/datahub/contents/public/{study}",
                                   timeout=60).json()
            expr = next(x["name"] for x in listing if x["name"].startswith("data_mrna") and "zscores" not in x["name"])
        for name, dest in ((expr, f"{key}_expression.txt"), ("data_clinical_patient.txt", f"{key}_clinical_patient.txt"),
                           ("data_clinical_sample.txt", f"{key}_clinical_sample.txt")):
            _get(DATAHUB.format(study, name), EXTRA / dest)
            if (EXTRA / dest).stat().st_size < 1000 and (EXTRA / dest).read_text(errors="ignore").startswith("version"):
                raise RuntimeError(f"{dest}: got a git-lfs pointer instead of the file")


def _singlecell_extra() -> None:
    sc = RAW / "singlecell"
    for name, asset in SINGLECELL_EXTRA.items():
        _get(CXG.format(asset), sc / name)
    for section, ds in VISIUM_SECTIONS.items():
        dest = sc / f"li2022_visium_section_{section}.h5ad"
        if dest.exists() and dest.stat().st_size > 0:
            continue
        assets = requests.get(LI2022.format(ds), timeout=60).json()["assets"]
        _get(next(a["url"] for a in assets if a["filetype"] == "H5AD"), dest)


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
    _treatment_studies()
    _singlecell_extra()
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
                                            "cBioPortal datahub" if f.name.split("_")[0] in TREATMENT_STUDIES else
                                            dict(targets).get(f.name.removesuffix(".gz"), "derived"))
        p = prev.get(f.name)
        same = p and p.get("bytes") == f.stat().st_size
        entries.append({"file": f.name, "url": url,
                        "downloaded_utc": p["downloaded_utc"] if same else
                        datetime.now(timezone.utc).isoformat(timespec="seconds"),
                        "sha256": p["sha256"] if same else _sha256(f), "bytes": f.stat().st_size})
    for f in sorted((RAW / "singlecell").glob("*.h5ad")):
        if not (f.name in SINGLECELL_EXTRA or f.name.startswith("li2022_visium_section_")):
            continue
        p = prev.get("singlecell/" + f.name)
        same = p and p.get("bytes") == f.stat().st_size
        entries.append({"file": "singlecell/" + f.name, "url": "CZ CELLxGENE Discover",
                        "downloaded_utc": p["downloaded_utc"] if same else
                        datetime.now(timezone.utc).isoformat(timespec="seconds"),
                        "sha256": p["sha256"] if same else _sha256(f), "bytes": f.stat().st_size})
    manifest = {"source": "M13 extra layers (cBioPortal, DepMap, UCSC Xena, cptac Zenodo, OmniPath, Open Targets)",
                "files": entries}
    MANIFEST.write_text(json.dumps(manifest, indent=1))
    return {"files": len(entries), "bytes": sum(e["bytes"] for e in entries)}


if __name__ == "__main__":
    print(json.dumps(download_extra(), indent=1))
