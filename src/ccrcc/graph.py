"""Embedded knowledge graph: SQLite storage + NetworkX view behind a small `GraphStore` interface.

Every relation is reified: (subject) -[predicate]-> (object) is stored as a claim that MUST reference an
Evidence record (type, source, statistics, date). Swapping the backend for Neo4j only requires another
implementation of this interface; agents (via kg-mcp) and the report only use these methods.
"""
from __future__ import annotations

import hashlib
import json
import sqlite3
from contextlib import contextmanager
from pathlib import Path

import networkx as nx

from .paths import GRAPH_DB

SCHEMA = """
CREATE TABLE IF NOT EXISTS nodes (
    id TEXT PRIMARY KEY, type TEXT NOT NULL, name TEXT NOT NULL, props TEXT NOT NULL DEFAULT '{}');
CREATE TABLE IF NOT EXISTS evidence (
    id TEXT PRIMARY KEY, type TEXT NOT NULL, source TEXT NOT NULL, date TEXT,
    props TEXT NOT NULL DEFAULT '{}');
CREATE TABLE IF NOT EXISTS claims (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    subject TEXT NOT NULL REFERENCES nodes(id),
    predicate TEXT NOT NULL,
    object TEXT NOT NULL REFERENCES nodes(id),
    evidence_id TEXT NOT NULL REFERENCES evidence(id),
    UNIQUE(subject, predicate, object, evidence_id));
CREATE INDEX IF NOT EXISTS ix_claims_subject ON claims(subject);
CREATE INDEX IF NOT EXISTS ix_claims_object ON claims(object);
CREATE INDEX IF NOT EXISTS ix_claims_evidence ON claims(evidence_id);
CREATE INDEX IF NOT EXISTS ix_nodes_name ON nodes(name);
CREATE INDEX IF NOT EXISTS ix_nodes_type ON nodes(type);
"""

NODE_TYPES = {"Gene", "Protein", "Pathway", "Drug", "Disease", "Dataset", "Publication", "Chunk", "Concept",
              "CellType", "Pocket"}


def evidence_id(ev_type: str, *parts) -> str:
    h = hashlib.sha1("|".join(map(str, parts)).encode()).hexdigest()[:12]
    return f"EV:{ev_type}:{h}"


class GraphStore:
    def __init__(self, path: Path = GRAPH_DB):
        self.path = Path(path)
        self.conn = sqlite3.connect(self.path)
        self.conn.row_factory = sqlite3.Row
        self.conn.execute("PRAGMA foreign_keys = ON")
        self.conn.executescript(SCHEMA)

    # ------------------------------------------------------------ writing
    @contextmanager
    def batch(self):
        try:
            yield self
            self.conn.commit()
        except Exception:
            self.conn.rollback()
            raise

    def reset(self) -> None:
        self.conn.executescript("DROP TABLE IF EXISTS claims; DROP TABLE IF EXISTS evidence; DROP TABLE IF EXISTS nodes;")
        self.conn.executescript(SCHEMA)
        self.conn.commit()

    def add_node(self, node_id: str, node_type: str, name: str, **props) -> str:
        if node_type not in NODE_TYPES:
            raise ValueError(f"unknown node type {node_type}")
        self.conn.execute(
            "INSERT INTO nodes(id, type, name, props) VALUES (?,?,?,?) "
            "ON CONFLICT(id) DO UPDATE SET props = json_patch(nodes.props, excluded.props)",
            (node_id, node_type, name, json.dumps(props, default=str)))
        return node_id

    def add_evidence(self, ev_id: str, ev_type: str, source: str, date: str | None = None, **props) -> str:
        # Ids hash (type, key parts), not values: a reloaded layer must replace the old values, not keep them silently.
        self.conn.execute("INSERT INTO evidence(id, type, source, date, props) VALUES (?,?,?,?,?) "
                          "ON CONFLICT(id) DO UPDATE SET type = excluded.type, source = excluded.source, "
                          "date = excluded.date, props = excluded.props",
                          (ev_id, ev_type, source, date, json.dumps(props, default=str)))
        return ev_id

    def add_claim(self, subject: str, predicate: str, obj: str, ev_id: str) -> None:
        self.conn.execute("INSERT OR IGNORE INTO claims(subject, predicate, object, evidence_id) VALUES (?,?,?,?)",
                          (subject, predicate, obj, ev_id))

    # ------------------------------------------------------------ reading
    @staticmethod
    def _node(row) -> dict:
        return {"id": row["id"], "type": row["type"], "name": row["name"], **json.loads(row["props"])}

    @staticmethod
    def _evidence(row) -> dict:
        return {"evidence_id": row["id"], "type": row["type"], "source": row["source"], "date": row["date"],
                **json.loads(row["props"])}

    def get_node(self, node_id: str) -> dict | None:
        r = self.conn.execute("SELECT * FROM nodes WHERE id = ?", (node_id,)).fetchone()
        return self._node(r) if r else None

    def resolve(self, entity: str) -> dict | None:
        """Resolve an id or a name (case-insensitive; genes preferred) to a node."""
        n = self.get_node(entity)
        if n:
            return n
        r = self.conn.execute(
            "SELECT * FROM nodes WHERE lower(name) = lower(?) "
            "ORDER BY CASE type WHEN 'Gene' THEN 0 WHEN 'Disease' THEN 1 WHEN 'Drug' THEN 2 ELSE 3 END LIMIT 1",
            (entity,)).fetchone()
        return self._node(r) if r else None

    def find_nodes(self, query: str, node_type: str | None = None, limit: int = 20) -> list[dict]:
        sql = "SELECT * FROM nodes WHERE name LIKE ?" + (" AND type = ?" if node_type else "") + " LIMIT ?"
        args = [f"%{query}%"] + ([node_type] if node_type else []) + [limit]
        return [self._node(r) for r in self.conn.execute(sql, args)]

    def claims_for(self, node_id: str, predicate: str | None = None, limit: int = 200) -> list[dict]:
        sql = ("SELECT c.subject, c.predicate, c.object, c.evidence_id, s.name AS sname, o.name AS oname, "
               "s.type AS stype, o.type AS otype FROM claims c JOIN nodes s ON s.id = c.subject "
               "JOIN nodes o ON o.id = c.object WHERE (c.subject = ? OR c.object = ?)"
               + (" AND c.predicate = ?" if predicate else "") + " ORDER BY c.id LIMIT ?")
        args = [node_id, node_id] + ([predicate] if predicate else []) + [limit]
        return [dict(r) for r in self.conn.execute(sql, args)]

    def neighbors(self, node_id: str, limit: int = 100) -> list[dict]:
        out = {}
        for c in self.claims_for(node_id, limit=limit * 5):
            other, oname, otype = ((c["object"], c["oname"], c["otype"]) if c["subject"] == node_id
                                   else (c["subject"], c["sname"], c["stype"]))
            key = (other, c["predicate"])
            out.setdefault(key, {"id": other, "name": oname, "type": otype, "predicate": c["predicate"],
                                 "direction": "out" if c["subject"] == node_id else "in", "evidence_ids": []})
            out[key]["evidence_ids"].append(c["evidence_id"])
        return list(out.values())[:limit]

    def get_evidence(self, ids: list[str]) -> list[dict]:
        if not ids:
            return []
        q = ",".join("?" * len(ids))
        return [self._evidence(r) for r in self.conn.execute(f"SELECT * FROM evidence WHERE id IN ({q})", ids)]

    def evidence_for(self, node_id: str, limit: int = 100) -> list[dict]:
        claims = self.claims_for(node_id, limit=limit)
        ev = {e["evidence_id"]: e for e in self.get_evidence([c["evidence_id"] for c in claims])}
        return [{"subject": c["sname"], "predicate": c["predicate"], "object": c["oname"],
                 **ev.get(c["evidence_id"], {"evidence_id": c["evidence_id"]})} for c in claims]

    def check_evidence_ids(self, ids: list[str]) -> dict:
        found = {e["evidence_id"] for e in self.get_evidence(list(set(ids)))}
        return {"valid": sorted(found), "invalid": sorted(set(ids) - found)}

    def stats(self) -> dict:
        c = self.conn
        return {
            "nodes_by_type": dict(c.execute("SELECT type, count(*) FROM nodes GROUP BY type").fetchall()),
            "evidence_by_type": dict(c.execute("SELECT type, count(*) FROM evidence GROUP BY type").fetchall()),
            "claims_by_predicate": dict(c.execute("SELECT predicate, count(*) FROM claims GROUP BY predicate").fetchall()),
            "claims_total": c.execute("SELECT count(*) FROM claims").fetchone()[0],
            "claims_without_resolvable_evidence": c.execute(
                "SELECT count(*) FROM claims c LEFT JOIN evidence e ON e.id = c.evidence_id "
                "WHERE e.id IS NULL OR e.source IS NULL OR e.source = ''").fetchone()[0],
        }

    def to_networkx(self, node_ids: set[str] | None = None) -> nx.MultiDiGraph:
        g = nx.MultiDiGraph()
        for r in self.conn.execute("SELECT * FROM nodes"):
            if node_ids is None or r["id"] in node_ids:
                g.add_node(r["id"], **self._node(r))
        for r in self.conn.execute("SELECT * FROM claims"):
            if node_ids is None or (r["subject"] in node_ids and r["object"] in node_ids):
                g.add_edge(r["subject"], r["object"], key=r["evidence_id"], predicate=r["predicate"],
                           evidence_id=r["evidence_id"])
        return g
