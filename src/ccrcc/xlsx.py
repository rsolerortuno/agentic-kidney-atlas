"""Minimal .xlsx reader (shared strings + inline values), so supplementary tables need no extra dependency.

Only what published supplementary tables use: one header row per sheet, text and numbers. Formulas are read as their
cached values; dates stay numbers.
"""
from __future__ import annotations

import re
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

import pandas as pd

NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
      "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships"}
REL_NS = "{http://schemas.openxmlformats.org/package/2006/relationships}"


def _col(ref: str) -> int:
    n = 0
    for ch in re.match(r"[A-Z]+", ref).group(0):
        n = n * 26 + ord(ch) - 64
    return n - 1


def sheet_names(path: Path) -> list[str]:
    with zipfile.ZipFile(path) as z:
        wb = ET.fromstring(z.read("xl/workbook.xml"))
    return [s.get("name") for s in wb.find("m:sheets", NS)]


def read_sheet(path: Path, sheet: str | int = 0, header: int = 0) -> pd.DataFrame:
    with zipfile.ZipFile(path) as z:
        wb = ET.fromstring(z.read("xl/workbook.xml"))
        sheets = wb.find("m:sheets", NS)
        el = sheets[sheet] if isinstance(sheet, int) else next(s for s in sheets if s.get("name") == sheet)
        rid = el.get(f"{{{NS['r']}}}id")
        rels = ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
        target = next(r.get("Target") for r in rels.iter(f"{REL_NS}Relationship") if r.get("Id") == rid)
        target = target.lstrip("/")
        target = target if target.startswith("xl/") else "xl/" + target
        shared = []
        if "xl/sharedStrings.xml" in z.namelist():
            for si in ET.fromstring(z.read("xl/sharedStrings.xml")).iterfind("m:si", NS):
                shared.append("".join(t.text or "" for t in si.iter(f"{{{NS['m']}}}t")))
        rows = []
        # streamed: supplementary expression matrices are hundreds of MB of XML
        with z.open(target) as fh:
            for _, row in ET.iterparse(fh):
                if row.tag != f"{{{NS['m']}}}row":
                    continue
                vals = {}
                for c in row.iterfind("m:c", NS):
                    t, v = c.get("t"), c.find("m:v", NS)
                    if t == "s" and v is not None:
                        val = shared[int(v.text)]
                    elif t == "inlineStr":
                        val = "".join(x.text or "" for x in c.iter(f"{{{NS['m']}}}t"))
                    elif v is not None:
                        try:
                            val = float(v.text)
                        except ValueError:
                            val = v.text
                    else:
                        val = None
                    vals[_col(c.get("r"))] = val
                rows.append(vals)
                row.clear()
    width = max((max(r) + 1 for r in rows if r), default=0)
    grid = [[r.get(i) for i in range(width)] for r in rows]
    cols = [str(c) if c is not None else f"col{i}" for i, c in enumerate(grid[header])]
    return pd.DataFrame(grid[header + 1:], columns=cols)
