"""Config loader."""
from functools import lru_cache

import yaml

from .paths import CONFIG


@lru_cache(maxsize=1)
def load() -> dict:
    return yaml.safe_load(CONFIG.read_text(encoding="utf-8"))
