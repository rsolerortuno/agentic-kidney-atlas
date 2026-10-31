"""M5 checks: one test per MCP tool (real data, external calls served from the disk cache) + stdio smoke test."""
import asyncio
import sys

import pytest
from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client

from ccrcc.mcp_servers import external as ext_mcp
from ccrcc.mcp_servers import kg as kg_mcp
from ccrcc.mcp_servers import omics as omics_mcp
from ccrcc.paths import GRAPH_DB, ROOT

pytestmark = pytest.mark.skipif(not GRAPH_DB.exists(), reason="build the graph first")


# ---------------------------------------------------------------- omics-mcp
def test_get_de_results_accepts_legacy_symbol():
    rows = omics_mcp.get_de_results("NDUFA4L2")
    by_layer = {r["layer"]: r for r in rows}
    assert by_layer["protein"]["gene"] == "COXFA4L2"
    assert by_layer["protein"]["significant"] and by_layer["protein"]["log2fc"] > 0
    assert by_layer["protein"]["evidence_id"].startswith("EV:omics_de:")


def test_top_genes():
    rows = omics_mcp.top_genes("protein", 10, "up")
    assert len(rows) == 10 and all(r["log2fc"] > 0 and r["significant"] for r in rows)


def test_rna_protein_correlation():
    r = omics_mcp.rna_protein_correlation("CA9")
    assert r["n_tumors"] >= 20 and -1 <= r["spearman_rho"] <= 1


def test_cnv_frequency():
    r = omics_mcp.cnv_frequency("VHL")
    assert r["arm"] == "3p" and r["loss_freq"] > 0.5 and r["evidence_id"]


# ---------------------------------------------------------------- kg-mcp
def test_neighbors():
    r = kg_mcp.neighbors("CA9", predicate="targets")
    assert r["total"] > 0 and all(n["evidence_ids"] for n in r["neighbors"])


def test_evidence_for():
    r = kg_mcp.evidence_for("CA9", evidence_type="literature_relation")
    assert r["total"] > 0 and all(e["pmid"] and e["sentence"] for e in r["evidence"])


def test_search_literature():
    hits = kg_mcp.search_literature("CA9 immunostain", 5)
    assert hits and hits[0]["evidence_id"] and "pubmed" in hits[0]["url"]


def test_check_and_get_evidence():
    ev = kg_mcp.evidence_for("CA9", limit=3)["evidence"][0]["evidence_id"]
    assert kg_mcp.check_evidence_ids([ev, "EV:fake:000"]) == {"valid": [ev], "invalid": ["EV:fake:000"]}
    assert kg_mcp.get_evidence([ev])[0]["evidence_id"] == ev


def test_convergence_ranking_tool():
    assert isinstance(kg_mcp.convergence_ranking(5), list)


# ---------------------------------------------------------------- external-mcp
@pytest.mark.network
def test_open_targets_association():
    r = ext_mcp.open_targets_association("VHL")
    assert r["rank"] <= 5 and r["evidence_id"]


@pytest.mark.network
def test_open_targets_target():
    r = ext_mcp.open_targets_target("CA9")
    assert "AB" in r["tractable_modalities"]
    assert any(d["name"] == "GIRENTUXIMAB" and d["evidence_id"] for d in r["drugs"])


@pytest.mark.network
def test_string_interactions():
    assert any(p["partner"] in ("HIF1A", "EPAS1") for p in ext_mcp.string_interactions("CA9"))


@pytest.mark.network
def test_uniprot_entry():
    assert ext_mcp.uniprot_entry("CA9")["accession"] == "Q16790"


def test_disease_top_targets():
    top = [r["gene"] for r in ext_mcp.disease_top_targets(10)]
    assert "VHL" in top and "PBRM1" in top


# ---------------------------------------------------------------- stdio protocol smoke test
@pytest.mark.parametrize("module,expected", [
    ("ccrcc.mcp_servers.omics", {"get_de_results", "top_genes", "rna_protein_correlation", "cnv_frequency"}),
    ("ccrcc.mcp_servers.kg", {"neighbors", "evidence_for", "search_literature", "check_evidence_ids"}),
    ("ccrcc.mcp_servers.external", {"open_targets_target", "string_interactions", "uniprot_entry"}),
    ("ccrcc.mcp_servers.cell", {"celltype_expression", "spatial_context", "structure_pockets", "mgd_ranking"}),
])
def test_stdio_server_lists_tools(module, expected):
    async def run():
        params = StdioServerParameters(command=sys.executable, args=["-m", module], cwd=str(ROOT))
        async with stdio_client(params) as (r, w):
            async with ClientSession(r, w) as s:
                await s.initialize()
                return {t.name for t in (await s.list_tools()).tools}
    assert expected <= asyncio.run(run())
