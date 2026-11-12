"""External knowledge clients (Open Targets, STRING, UniProt) with a disk cache.

Every response is cached as JSON under data/cache/<source>/ so reruns and tests are reproducible offline.
"""
from __future__ import annotations

import hashlib
import json
import time

import requests

from .config import load
from .paths import CACHE

OT_URL = "https://api.platform.opentargets.org/api/v4/graphql"
STRING_URL = "https://string-db.org/api/json/interaction_partners"
UNIPROT_URL = "https://rest.uniprot.org/uniprotkb/search"


# Where a new answer is written. The pipeline writes to data/cache, which is frozen in reproducibility/inputs.json.
# The agents' tools (mcp_servers/external.py) set this to data/cache_agent: a lookup made during an agent run for a
# gene the pipeline never queried used to add a file to the frozen tree, after which the input check and `run-all`
# failed ("unlisted input in frozen tree").
WRITE_ROOT = CACHE


def _cached(source: str, key: str, fetch):
    name = hashlib.sha1(key.encode()).hexdigest()[:16] + ".json"
    for root in dict.fromkeys((CACHE, WRITE_ROOT)):       # the frozen cache is always read first
        if (root / source / name).exists():
            return json.loads((root / source / name).read_text(encoding="utf-8"))
    path = WRITE_ROOT / source / name
    path.parent.mkdir(parents=True, exist_ok=True)
    data = fetch()
    path.write_text(json.dumps({"key": key, "retrieved_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                                "data": data}), encoding="utf-8")
    return json.loads(path.read_text(encoding="utf-8"))


def _request(method: str, url: str, **kw) -> requests.Response:
    attempts = 8   # public APIs (UniProt, Open Targets) return transient 5xx; back off up to ~1 min
    for attempt in range(attempts):
        try:
            r = requests.request(method, url, timeout=60, **kw)
        except (requests.ConnectionError, requests.Timeout):
            if attempt == attempts - 1:
                raise
            time.sleep(min(2 ** attempt, 60))
            continue
        if (r.status_code == 429 or r.status_code >= 500) and attempt < attempts - 1:
            time.sleep(min(2 ** attempt, 60))
            continue
        r.raise_for_status()
        return r
    r.raise_for_status()
    return r


def _ot(query: str, variables: dict) -> dict:
    key = query + json.dumps(variables, sort_keys=True)
    return _cached("open_targets", key, lambda: _request(
        "POST", OT_URL, json={"query": query, "variables": variables}).json())["data"]


# ---------------------------------------------------------------- Open Targets
_ASSOC_Q = """query($id:String!,$index:Int!,$size:Int!){ disease(efoId:$id){ id name
  associatedTargets(page:{index:$index,size:$size}){ count rows{ score target{ id approvedSymbol }
  datatypeScores{ id score } } } } }"""
OT_DATATYPES = ("genetic_association", "somatic_mutation", "known_drug", "literature", "rna_expression",
                "animal_model", "affected_pathway")

_TARGET_Q = """query($id:String!){ target(ensemblId:$id){ id approvedSymbol approvedName
  tractability{ label modality value }
  pathways{ pathwayId pathway topLevelTerm }
  drugAndClinicalCandidates{ count rows{ maxClinicalStage drug{ id name } } } } }"""


def ot_disease_associations(n: int | None = None) -> list[dict]:
    cfg = load()["disease"]
    n = n or cfg["open_targets_top"]
    size, rows = 500, []
    for index in range(0, (n + size - 1) // size):
        res = _ot(_ASSOC_Q, {"id": cfg["open_targets_id"], "index": index, "size": size})
        page = res["data"]["disease"]["associatedTargets"]["rows"]
        for r in page:
            dt = {d["id"]: d["score"] for d in r.get("datatypeScores") or []}
            rows.append({"ensembl": r["target"]["id"], "symbol": r["target"]["approvedSymbol"], "score": r["score"],
                         **{f"dt_{k}": float(dt.get(k, 0.0)) for k in OT_DATATYPES}})
        if len(page) < size:
            break
    return rows[:n]


def ot_target(ensembl_id: str) -> dict | None:
    t = _ot(_TARGET_Q, {"id": ensembl_id})["data"]["target"]
    if t is None:
        return None
    tract = [x for x in (t.get("tractability") or []) if x["value"]]
    drugs = (t.get("drugAndClinicalCandidates") or {}).get("rows") or []
    return {
        "ensembl": t["id"], "symbol": t["approvedSymbol"], "name": t["approvedName"],
        "tractable_modalities": sorted({x["modality"] for x in tract}),
        "tractability_labels": sorted({f"{x['modality']}: {x['label']}" for x in tract}),
        "has_approved_or_clinical_drug": any(d["maxClinicalStage"] in {"APPROVAL", "PHASE_4", "PHASE_3", "PHASE_2", "PHASE_1"}
                                             for d in drugs),
        "drugs": [{"chembl": d["drug"]["id"], "name": d["drug"]["name"], "max_stage": d["maxClinicalStage"]}
                  for d in drugs],
        "pathways": [{"reactome": p["pathwayId"], "name": p["pathway"], "top_level": p["topLevelTerm"]}
                     for p in (t.get("pathways") or [])],
    }


# ---------------------------------------------------------------- STRING
def string_partners(symbol: str, limit: int = 10, min_score: float = 0.7) -> list[dict]:
    params = {"identifiers": symbol, "species": 9606, "limit": limit, "caller_identity": "ccrcc-graphrag"}
    def fetch():
        try:
            return _request("GET", STRING_URL, params=params).json()
        except requests.HTTPError as e:
            if e.response is not None and e.response.status_code in (400, 404):  # identifier unknown to STRING
                return []
            raise

    data = _cached("string", json.dumps(params, sort_keys=True), fetch)["data"]
    if isinstance(data, dict):  # error payload
        return []
    return [{"partner": r["preferredName_B"], "score": r["score"], "experimental": r["escore"],
             "database": r["dscore"], "textmining": r["tscore"]} for r in data if r["score"] >= min_score]


# ---------------------------------------------------------------- UniProt
def uniprot(symbol: str) -> dict | None:
    params = {"query": f"gene_exact:{symbol} AND organism_id:9606 AND reviewed:true",
              "fields": "accession,protein_name,cc_function,cc_subcellular_location", "format": "json", "size": 1}
    data = _cached("uniprot", json.dumps(params, sort_keys=True),
                   lambda: _request("GET", UNIPROT_URL, params=params).json())["data"]
    if not data.get("results"):
        return None
    r = data["results"][0]
    name = r.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value")
    function, locations = "", []
    for c in r.get("comments", []):
        if c.get("commentType") == "FUNCTION":
            function = " ".join(t["value"] for t in c.get("texts", []))
        elif c.get("commentType") == "SUBCELLULAR LOCATION":
            locations += [s["location"]["value"] for s in c.get("subcellularLocations", []) if "location" in s]
    return {"accession": r["primaryAccession"], "name": name, "function": function[:800],
            "subcellular_location": sorted(set(locations))}


# ---------------------------------------------------------------- Europe PMC (literature volume)
EPMC_URL = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
CCRCC_TERMS = ('("renal cell carcinoma" OR "clear cell renal" OR "ccRCC" OR "kidney cancer" OR "renal carcinoma" '
               'OR "renal cancer")')


def europepmc_count(query: str) -> int:
    """Number of Europe PMC records matching a query (cached; used to size the literature around a gene)."""
    params = {"query": query, "format": "json", "pageSize": 1, "resultType": "lite"}
    data = _cached("europepmc", json.dumps(params, sort_keys=True),
                   lambda: _request("GET", EPMC_URL, params=params).json())["data"]
    time.sleep(0.05)
    return int(data["hitCount"])


# ---------------------------------------------------------------- ChEMBL (preclinical chemical matter and patents)
CHEMBL_URL = "https://www.ebi.ac.uk/chembl/api/data/"
PATENT_SOURCES = {37, 38}    # BindingDB patent bioactivities, SureChEMBL patent bioactivities


def _chembl(endpoint: str, params: dict) -> dict:
    params = {**params, "format": "json"}
    return _cached("chembl", endpoint + json.dumps(params, sort_keys=True),
                   lambda: _request("GET", CHEMBL_URL + endpoint, params=params).json())["data"]


def chembl_target_summary(accession: str, min_pchembl: float = 6.0) -> dict:
    """Potent compounds (pChEMBL >= min_pchembl, i.e. <= 1 uM) reported against the single-protein ChEMBL target of a
    UniProt accession. `potent_activities` counts all such records; molecules and patent records (BindingDB / SureChEMBL
    patent sources) are counted in the first and last 1,000 records, which is enough to tell 'none', 'a few' and 'a series'."""
    t = _chembl("target.json", {"target_components__accession": accession, "target_type": "SINGLE PROTEIN",
                                "organism": "Homo sapiens", "limit": 5})
    targets = [x["target_chembl_id"] for x in t.get("targets", [])]
    if not targets:
        return {"chembl_target": None, "potent_activities": 0, "potent_molecules": 0, "patent_activities": 0}
    acts = _chembl("activity.json", {"target_chembl_id": targets[0], "pchembl_value__gte": min_pchembl, "limit": 1000,
                                     "only": "molecule_chembl_id,src_id"})
    rows = acts.get("activities", [])
    total = int(acts["page_meta"]["total_count"])
    if total > 1000:   # patent records were loaded last (highest activity ids): look at the last page as well
        rows = rows + _chembl("activity.json", {"target_chembl_id": targets[0], "pchembl_value__gte": min_pchembl,
                                                "limit": 1000, "offset": total - 1000,
                                                "only": "molecule_chembl_id,src_id"}).get("activities", [])
    return {"chembl_target": targets[0], "potent_activities": int(acts["page_meta"]["total_count"]),
            "potent_molecules": len({a["molecule_chembl_id"] for a in rows}),
            "patent_activities": sum(1 for a in rows if a.get("src_id") in PATENT_SOURCES)}
