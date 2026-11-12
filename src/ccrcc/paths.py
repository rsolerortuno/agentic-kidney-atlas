"""Canonical project paths."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data"
RAW = DATA / "raw"
PROCESSED = DATA / "processed"
CACHE = DATA / "cache"
AGENT_CACHE = DATA / "cache_agent"   # answers fetched by agent tools; not part of the frozen inputs
LITERATURE = DATA / "literature"
GRAPH_DB = DATA / "graph" / "kg.sqlite"
RUNS = ROOT / "runs"
REPORTS = ROOT / "reports"
CONFIG = ROOT / "config.yaml"

for _d in (RAW, PROCESSED, CACHE, LITERATURE, GRAPH_DB.parent, RUNS, REPORTS):
    _d.mkdir(parents=True, exist_ok=True)
