"""Validate reports/ccrcc_report.html:
- every claim element (class "claim") carries >= 1 evidence id;
- every evidence id exists in the knowledge graph and has an anchor in the evidence appendix;
- every in-page citation link resolves to an anchor;
- the file is self-contained (no external <script src> / <link href> needed to open it).
Exit code 1 on any violation.
"""
import re
import sys
from html.parser import HTMLParser

from ccrcc.graph import GraphStore
from ccrcc.paths import REPORTS
from ccrcc.render import anchor


class Scan(HTMLParser):
    def __init__(self):
        super().__init__()
        self.claims, self.anchors, self.links, self.external = [], set(), [], []

    def handle_starttag(self, tag, attrs):
        a = dict(attrs)
        if "claim" in (a.get("class") or "").split():
            self.claims.append(a.get("data-evidence", "").split())
        if a.get("id"):
            self.anchors.add(a["id"])
        if tag == "a" and (a.get("href") or "").startswith("#"):
            self.links.append(a["href"][1:])
        if tag == "script" and a.get("src"):
            self.external.append(a["src"])
        if tag == "link" and a.get("rel") == "stylesheet" and a.get("href"):
            self.external.append(a["href"])


def main(path=REPORTS / "ccrcc_report.html") -> int:
    s = Scan()
    s.feed(path.read_text(encoding="utf-8"))
    errors = []
    empty = sum(1 for c in s.claims if not c)
    if empty:
        errors.append(f"{empty} claim(s) without evidence ids")
    ids = sorted({i for c in s.claims for i in c})
    chk = GraphStore().check_evidence_ids(ids)
    if chk["invalid"]:
        errors.append(f"{len(chk['invalid'])} evidence id(s) not in graph: {chk['invalid'][:5]}")
    missing_anchor = [i for i in ids if anchor(i) not in s.anchors]
    if missing_anchor:
        errors.append(f"{len(missing_anchor)} evidence id(s) without appendix anchor: {missing_anchor[:5]}")
    broken = [l for l in s.links if l not in s.anchors]
    if broken:
        errors.append(f"{len(broken)} broken in-page link(s): {broken[:5]}")
    if s.external:
        errors.append(f"external resources required: {s.external}")
    print(f"claims: {len(s.claims)} | distinct evidence ids: {len(ids)} | citation links: {len(s.links)} | "
          f"external resources: {len(s.external)}")
    for e in errors:
        print("ERROR:", e)
    print("report valid" if not errors else "report INVALID")
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())
