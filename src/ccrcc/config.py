"""Config loader."""
from functools import lru_cache

import yaml

from .paths import CONFIG


@lru_cache(maxsize=1)
def load() -> dict:
    return yaml.safe_load(CONFIG.read_text(encoding="utf-8"))


def analysis_targets() -> list[str]:
    """Genes that get cell-type, spatial and structure analysis: the configured report candidates plus the
    discovery shortlist (data/processed/discovery_shortlist.json, written by `ccrcc discovery`)."""
    import json

    from .paths import PROCESSED
    genes = list(load()["cell_context"]["targets"])
    short = PROCESSED / "discovery_shortlist.json"
    if short.exists():
        genes += [g for g in json.loads(short.read_text()) if g not in genes]
    return genes
