"""M13 evaluation against outcomes the scores do not see (deterministic, no LLM).

Three gold sets, each matched to the question a score answers:

1. RCC clinical targets (Open Targets 26.09 clinical precedence): targets of drugs tested in renal cell carcinoma
   (papillary / chromophobe / HLRCC excluded) in phase 2 or later, with >= 2 trials or an approval. Drugs whose
   mechanism maps to more than 5 targets (e.g. metformin -> 50 complex I subunits) are removed so that one drug
   does not contribute a block of genes.
2. Retrospective post-2016 evaluation: the subset of (1) whose first RCC trial started on or after 2016-01-01. The data-only
   profiles use measurements, not trial history or literature, so this asks whether they rank genes that entered
   RCC clinical development later.
3. Prognostic biomarkers: the ClearCode34 ccA/ccB subtype classifier (Brooks et al., Eur Urol 2014, Table 1; 34 genes,
   prognosis validated in TCGA). A first choice — verbatim "biomarker_of ccRCC" relations from the disease-level PubMed
   query — turned out empty: every such relation came from the gene-specific queries, which were chosen from the
   omics ranking and would make the test circular. Caveat: ClearCode34 was validated in TCGA, the cohort the
   biomarker profile's prognosis component also uses, so the biomarker AUROC is partly optimistic; the
   "without prognosis" row removes that component.

Metrics: AUROC over the whole gene universe (gold vs all other genes) with a bootstrap 95 % CI over gold genes,
recall@100 / @500 and fold enrichment. The convergence score (which includes literature and drug history) is
reported next to the profiles to show how much circularity inflates it.
"""
from __future__ import annotations

import json
import re

import numpy as np
import pandas as pd

from .literature_entities import is_ccrcc, resolve_gene
from .paths import LITERATURE, RAW

EXTRA = RAW / "extra"
RCC = re.compile(r"renal cell|kidney cancer|renal cancer|renal carcinoma|clear cell renal|kidney carcinoma|\bRCC\b", re.I)
NOT_CC = re.compile(r"papillary|chromophobe|leiomyomatosis|translocation", re.I)
ADVANCED = {"PHASE_2", "PHASE_2_3", "PHASE_3", "PHASE_4", "APPROVAL"}
CUTOFF = "2016-01-01"
MAX_TARGETS_PER_DRUG = 5


def _symbols() -> dict:
    hg = pd.read_csv(RAW / "hgnc_complete_set.txt", sep="\t", usecols=["symbol", "ensembl_gene_id"], dtype=str).dropna()
    return dict(zip(hg["ensembl_gene_id"], hg["symbol"]))


def clinical_targets() -> pd.DataFrame:
    cp = pd.read_parquet(EXTRA / "ot_clinical_precedence_0.parquet",
                         columns=["targetId", "drugId", "clinicalReportId", "clinicalStage", "studyStartDate",
                                  "diseaseFromSource"])
    return clinical_targets_from_records(cp, _symbols())


def clinical_targets_from_records(cp: pd.DataFrame, symbols: dict) -> pd.DataFrame:
    """Retrospective grouping, not a prospective test: features are not frozen at the cutoff.

    Eligibility uses advanced trials; first observed RCC trial uses every phase and
    every drug, including drugs excluded from the eligibility set for promiscuity.
    Missing dates prevent assignment to the post-cutoff group.
    """
    n_targets = cp.groupby("drugId")["targetId"].nunique()
    src = cp["diseaseFromSource"].fillna("")
    r = cp[src.str.contains(RCC) & ~src.str.contains(NOT_CC)].copy()
    r["gene"] = r["targetId"].map(symbols)
    r["start"] = pd.to_datetime(r["studyStartDate"], errors="coerce")
    dates = r.dropna(subset=["gene"]).groupby("gene").agg(
        first_trial=("start", "min"), date_incomplete=("start", lambda s: s.isna().any()))
    eligible = r[r["clinicalStage"].isin(ADVANCED) & (r["drugId"].map(n_targets) <= MAX_TARGETS_PER_DRUG)]
    g = eligible.dropna(subset=["gene"]).groupby("gene").agg(
        trials=("clinicalReportId", "nunique"), drugs=("drugId", "nunique"),
        approved=("clinicalStage", lambda s: "APPROVAL" in set(s)))
    g = g[(g["trials"] >= 2) | g["approved"]]
    g = g.join(dates)
    g["holdout"] = (g["first_trial"] >= pd.Timestamp(CUTOFF)) & ~g["date_incomplete"]
    return g.reset_index()


CLEARCODE34 = {  # Brooks et al. 2014, Table 1 (symbols as published; mapped to current HGNC below)
    "ccA": ["MAPT", "STK32B", "FZD1", "RGS5", "GIPC2", "PDGFD", "EPAS1", "MAOB", "CDH5", "TCEA3", "LEPROTL1", "BNIP3L",
            "EHBP1", "VCAM1", "PHYH", "PRKAA2", "SLC4A4", "ESD", "TLR3", "NRP1", "C11orf1", "ST13", "ARNT", "C13orf1"],
    "ccB": ["SERPINA3", "SLC4A3", "MOXD1", "KCNN4", "ROR2", "FLJ23867", "FOXM1", "UNG2", "GALNT10", "GALNT4"]}


def clearcode34() -> list[str]:
    out = []
    for g in CLEARCODE34["ccA"] + CLEARCODE34["ccB"]:
        r = resolve_gene(g) or {"UNG2": "UNG"}.get(g)
        if r and r not in out:
            out.append(r)
    return out


def literature_biomarkers() -> list[str]:
    q = json.loads((LITERATURE / "query_log.json").read_text())
    disease = set(q["disease_pmids"])
    out = set()
    for f in sorted((LITERATURE / "extracted").glob("*.jsonl")):
        for line in f.read_text(encoding="utf-8").splitlines():
            if not line.strip():
                continue
            r = json.loads(line)
            if r.get("relation") == "biomarker_of" and r.get("pmid") in disease and is_ccrcc(r.get("object", "")):
                g = resolve_gene(r.get("subject", ""))
                if g:
                    out.add(g)
    return sorted(out)


def auroc(score: np.ndarray, gold: np.ndarray) -> float:
    pos, neg = score[gold], score[~gold]
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    ranks = pd.Series(np.concatenate([pos, neg])).rank().to_numpy()
    return float((ranks[:len(pos)].sum() - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg)))


def metrics(t: pd.DataFrame, col: str, gold: set[str], exclude: set[str] | None = None, seed: int = 3) -> dict:
    d = t[~t["gene"].isin(exclude or set())]
    s = d[col].to_numpy(float)
    g = d["gene"].isin(gold).to_numpy()
    # Ties are broken at random, not by input order: the input is sorted by convergence (literature + drugs), which
    # would rank tied known targets first. A separate generator keeps the bootstrap draws below unchanged.
    order = np.lexsort((np.random.default_rng(seed + 1).random(len(s)), -s))
    rank = np.empty(len(s), int)
    rank[order] = np.arange(1, len(s) + 1)
    rng = np.random.default_rng(seed)
    pos_idx = np.nonzero(g)[0]
    boots = []
    for _ in range(500):
        sample = rng.choice(pos_idx, size=len(pos_idx), replace=True)
        boots.append(auroc(np.concatenate([s[sample], s[~g]]),
                           np.concatenate([np.ones(len(sample), bool), np.zeros((~g).sum(), bool)])))
    if not len(pos_idx):
        boots = [float("nan")]
    n = int(g.sum())
    rec100 = float((rank[g] <= 100).mean()) if n else float("nan")
    return {"n_gold": n, "universe": int(len(s)), "auroc": round(auroc(s, g), 3),
            "auroc_ci95": [round(float(np.nanpercentile(boots, 2.5)), 3), round(float(np.nanpercentile(boots, 97.5)), 3)],
            "recall_at_100": round(rec100, 3), "recall_at_500": round(float((rank[g] <= 500).mean()), 3) if n else None,
            "median_rank": int(np.median(rank[g])) if n else None,
            "fold_enrichment_at_100": round(rec100 / (100 / len(s)), 1) if n else None,
            "top_hits": d.iloc[order[:100]].loc[lambda x: x["gene"].isin(gold), "gene"].tolist()[:15]}


def run(t: pd.DataFrame) -> dict:
    """t = integrate.compute() output (scores and components for every gene)."""
    t = t.copy()
    from .config import load
    sw = load()["convergence"]["stream_weights"]
    t["score_data_only"] = t["score"] - t["c_literature"] * sw["literature"] - t["c_tractability"] * sw["tractability"]
    clin = clinical_targets()
    allc, hold = set(clin["gene"]), set(clin.loc[clin["holdout"], "gene"])
    pre = allc - hold
    bio = set(clearcode34())
    w = load()["profiles"]["biomarker"]
    t["score_biomarker_no_prognosis"] = t["score_biomarker"] - w.get("prognosis", 0) * t["c_prognosis"]
    scores = {"target profile": "score_target", "biomarker profile": "score_biomarker",
              "convergence v3 (incl. literature + drugs)": "score", "convergence v3 without literature/tractability":
              "score_data_only"}
    out = {"gold_sets": {"rcc_clinical_targets": sorted(allc), "temporal_holdout_targets": sorted(hold),
                         "clearcode34": sorted(bio), "literature_biomarkers_disease_query": literature_biomarkers(),
                         "cutoff": CUTOFF},
           "clinical_targets": {}, "temporal_holdout": {}, "clearcode34": {}}
    for label, col in scores.items():
        out["clinical_targets"][label] = metrics(t, col, allc)
        # holdout: earlier clinical targets are removed from the background so they neither help nor hurt
        out["temporal_holdout"][label] = metrics(t, col, hold, exclude=pre)
        out["clearcode34"][label] = metrics(t, col, bio)
    out["clearcode34"]["biomarker profile without prognosis"] = metrics(t, "score_biomarker_no_prognosis", bio)
    # ablations: how much of the target-profile result is "druggable proteins are what trials target"?
    wt = load()["profiles"]["target"]
    t["score_target_no_tractability"] = t["score_target"] - wt.get("tractability_predicted", 0) * t["c_tractability_predicted"]
    # binary with random tie order: recall@k then shows what tractability alone gives when nothing orders the tractable set
    noise = np.random.default_rng(0).random(len(t)) * 1e-9          # random order within ties (not the input order)
    t["tractability_only"] = t["c_tractability_predicted"] + noise
    out["n_predicted_tractable"] = int((t["c_tractability_predicted"] > 0).sum())
    for key, gold, excl in (("clinical_targets", allc, None), ("temporal_holdout", hold, pre)):
        out[key]["target profile without tractability"] = metrics(t, "score_target_no_tractability", gold, exclude=excl)
        out[key]["predicted tractability alone"] = metrics(t, "tractability_only", gold, exclude=excl)
    out["clinical_table"] = clin.assign(first_trial=clin["first_trial"].dt.date.astype(str)).to_dict("records")
    return out
