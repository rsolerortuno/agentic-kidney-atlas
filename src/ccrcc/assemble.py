"""Assemble runs/<date>/findings.json from deterministic convergence, verified agent findings and the
orchestrator narrative. Narrative statements are verified the same way as agent findings."""
from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

from .config import load
from .paths import CONFIG, RAW
from .verify import AGENT_FILES, Verifier

MODELS = {"orchestrator": "Claude Opus 5.5 (effort medium) — Claude Code main session",
          "agents": "Claude Sonnet 5.5 — Claude Code subagents (agente-omico, agente-literatura, agente-externo)",
          "relation_extraction": "Claude Sonnet 5.5 — extractor-relaciones subagent"}


def assemble(run_dir: Path) -> dict:
    cfg = load()
    conv = json.loads((run_dir / "convergence.json").read_text(encoding="utf-8"))
    agents = {}
    for name in AGENT_FILES:
        p = run_dir / name.replace(".json", ".verified.json")
        if p.exists():
            agents[name.removeprefix("agent_").removesuffix(".json")] = json.loads(p.read_text(encoding="utf-8"))

    v = Verifier()
    narrative_path = run_dir / "narrative.json"
    narrative = json.loads(narrative_path.read_text(encoding="utf-8")) if narrative_path.exists() else {}
    dropped = []

    def keep(paras, where):
        out = []
        for p in paras:
            ok, log = v.check_statement(p)
            (out.append(ok) if ok else dropped.append({"where": where, **log}))
        return out

    narrative_v = {"summary": keep(narrative.get("summary", []), "summary"),
                   "sections": {k: keep(v_, k) for k, v_ in narrative.get("sections", {}).items()},
                   "limitations": narrative.get("limitations", [])}
    targets = keep(narrative.get("targets", []), "targets")

    manifest = json.loads((RAW / "manifest.json").read_text())
    findings = {
        "question": cfg["question"],
        "disease": cfg["disease"],
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "models": MODELS,
        "config_sha256": hashlib.sha256(CONFIG.read_bytes()).hexdigest(),
        "datasets": [{k: f[k] for k in ("file", "url", "downloaded_utc", "sha256", "rows", "columns")}
                     for f in manifest["files"]],
        "convergence": {"universe_size": conv["universe_size"], "config": conv["config"],
                        "ranking": conv["ranking"]},
        "agent_findings": agents,
        "narrative": narrative_v,
        "targets": targets,
        "dropped_during_assembly": dropped,
    }
    out = run_dir / "findings.json"
    out.write_text(json.dumps(findings, indent=1, ensure_ascii=False), encoding="utf-8")
    return {"path": str(out), "agents": list(agents), "narrative_dropped": len(dropped)}
