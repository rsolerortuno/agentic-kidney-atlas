# PROGRESS

Spec: [OBJETIVO.md](OBJETIVO.md). A milestone is closed only when its check passes; evidence is recorded below.

**Current milestone:** M4

| Milestone | Status | Evidence |
|---|---|---|
| M0 Environment | ✅ passed | uv 0.13.6 (winget, user scope); Python 3.12.16 managed by uv; `uv run python -c "import pandas, mcp, networkx"` → ok; `uv run pytest` → 1 passed |
| M1 Data | ✅ passed | `data/raw/manifest.json` (URL, UTC date, sha256, dims). 110 tumors in every tumor matrix; normals: proteome 84, phospho 84, RNA 75. HGNC 45,289 rows. `tests/test_data.py` green |
| M2 Omics analysis | ✅ passed | `data/processed/*.parquet`. Paired t-test (matched cases) + BH. RNA: CA9 log2FC 6.31, NDUFA4L2 6.75, ANGPTL4 5.57, EGLN3 4.41, VEGFA 2.92 (all padj < 1e-32). Protein: CA9 2.27, NDUFA4L2 3.77, ANGPTL4 1.73 (padj < 1e-35). EGLN3 is not in the CPTAC proteome; VEGFA protein has 14 tumor / 4 normal values (not testable). 3p is the top arm loss (median loss frequency 0.80). ID mapping: 18,304 approved, 980 via previous/alias symbol, 684 unmapped. `tests/test_omics.py` green |
| M3 Literature | ✅ passed | PubTator3: 319 publications (top 100 disease hits + 5 per top-50 omics gene, searched with the literature-era symbol, e.g. NDUFA4L2 for HGNC COXFA4L2), 7,270 normalized annotations, 3,487 sentence chunks with a BM25 index (`ccrcc.search`). 100 abstracts in 4 batches processed by 4 `extractor-relaciones` subagents (Sonnet 5.5): 198 relations, 198 kept / 0 rejected by the verbatim validator (82 PMIDs). Validator negative tests (paraphrase, out-of-schema relation, foreign PMID, bad confidence, invalid JSON) pass. `tests/test_literature.py` green |
| M4 Knowledge graph | ⏳ | |
| M5 MCP servers | ⏳ | |
| M6 Agents | ⏳ | |
| M7 HTML report | ⏳ | |
| M8 Evaluation | ⏳ | |
| M9 Reproducibility | ⏳ | |

## Blockers / notes
- The original session folder path exceeded the Windows 260-character limit inside `.venv` (networkx import failures). The project was moved to `C:\Users\rafae\Documents\ccrcc-graphrag`. Enabling long paths is a system setting and was not touched.
- **BLOCKER (user action needed for M6/M8):** the headless Claude Code CLI bundled with the desktop app (`%LOCALAPPDATA%\Claude-Data\claude_pzs8sxrjxfjjc\claude-code\2.1.309\...\claude.exe`) answers `Not logged in · Please run /login`. Run it once in a terminal and execute `/login` (browser OAuth with the Pro account). Credentials cannot be entered on the user's behalf.
