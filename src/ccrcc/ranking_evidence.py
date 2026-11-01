"""Citable, content-addressed ranking results; never a pool of interchangeable numbers."""
from __future__ import annotations

import hashlib
import json

from .graph import GraphStore, evidence_id


def load(conv: dict, store: GraphStore | None = None) -> list[str]:
    store = store or GraphStore()
    config = {k:conv.get(k) for k in ("config", "profile_weights", "sensitivity_config")}
    config_hash = hashlib.sha256(json.dumps(config, sort_keys=True).encode()).hexdigest()
    ids = []
    with store.batch():
        for profile, rows in {"convergence":conv["ranking"], **conv.get("profiles", {})}.items():
            source = f"ANALYSIS:ranking:{profile}:{config_hash[:16]}"
            store.add_node(source, "Dataset", f"Computed {profile} ranking", config_sha256=config_hash)
            for row in rows:
                node = store.resolve(row["gene"])
                if not node or node["type"] != "Gene":
                    continue
                values = {k:row[k] for k in ("gene", "rank", "score", "rank_interval_90") if k in row}
                digest = hashlib.sha256(json.dumps(values, sort_keys=True).encode()).hexdigest()
                ev = evidence_id("analysis_ranking", profile, config_hash, digest)
                store.add_evidence(ev, "analysis_ranking", source, conv.get("generated"),
                                   **values, profile=profile, result_sha256=digest,
                                   interpretation="Computed ranking and weight sensitivity, not biological validation")
                store.add_claim(node["id"], "ranked_in", source, ev)
                ids.append(ev)
    return ids
