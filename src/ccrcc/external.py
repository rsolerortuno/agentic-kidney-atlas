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


def _cached(source: str, key: str, fetch):
    d = CACHE / source
    d.mkdir(parents=True, exist_ok=True)
    path = d / (hashlib.sha1(key.encode()).hexdigest()[:16] + ".json")
    if path.exists():
        return json.loads(path.read_text(encoding="utf-8"))
    data = fetch()
    path.write_text(json.dumps({"key": key, "retrieved_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                                "data": data}), encoding="utf-8")
    return json.loads(path.read_text(encoding="utf-8"))


def _request(method: str, url: str, **kw) -> requests.Response:
    for attempt in range(5):
        r = requests.request(method, url, timeout=60, **kw)
        if r.status_code == 429 or r.status_code >= 500:
            time.sleep(2 ** attempt)
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
  associatedTargets(page:{index:$index,size:$size}){ count rows{ score target{ id approvedSymbol } } } } }"""

_TARGET_Q = """query($id:String!){ target(ensemblId:$id){ id approvedSymbol approvedName
  tractability{ label modality value }
  drugAndClinicalCandidates{ count rows{ maxClinicalStage drug{ id name } } } } }"""


def ot_disease_associations(n: int | None = None) -> list[dict]:
    cfg = load()["disease"]
    n = n or cfg["open_targets_top"]
    size, rows = 500, []
    for index in range(0, (n + size - 1) // size):
        res = _ot(_ASSOC_Q, {"id": cfg["open_targets_id"], "index": index, "size": size})
        page = res["data"]["disease"]["associatedTargets"]["rows"]
        rows += [{"ensembl": r["target"]["id"], "symbol": r["target"]["approvedSymbol"], "score": r["score"]}
                 for r in page]
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
    }


# ---------------------------------------------------------------- STRING
def string_partners(symbol: str, limit: int = 10, min_score: float = 0.7) -> list[dict]:
    params = {"identifiers": symbol, "species": 9606, "limit": limit, "caller_identity": "ccrcc-graphrag"}
    data = _cached("string", json.dumps(params, sort_keys=True),
                   lambda: _request("GET", STRING_URL, params=params).json())["data"]
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
