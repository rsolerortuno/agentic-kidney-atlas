"""BM25 index over literature chunks (title/abstract sentences)."""
from __future__ import annotations

import re
from functools import lru_cache

import pandas as pd
from rank_bm25 import BM25Okapi

from .paths import LITERATURE

_TOKEN = re.compile(r"[A-Za-z0-9][A-Za-z0-9\-]+")


def tokenize(text: str) -> list[str]:
    return [t.lower() for t in _TOKEN.findall(text)]


@lru_cache(maxsize=1)
def _index() -> tuple[BM25Okapi, pd.DataFrame]:
    chunks = pd.read_parquet(LITERATURE / "chunks.parquet")
    return BM25Okapi([tokenize(t) for t in chunks["text"]]), chunks


def search_chunks(query: str, k: int = 10) -> list[dict]:
    bm25, chunks = _index()
    scores = bm25.get_scores(tokenize(query))
    top = scores.argsort()[::-1][:k]
    return [{"chunk_id": chunks.iloc[i]["chunk_id"], "pmid": int(chunks.iloc[i]["pmid"]),
             "text": chunks.iloc[i]["text"], "score": round(float(scores[i]), 3)} for i in top if scores[i] > 0]
