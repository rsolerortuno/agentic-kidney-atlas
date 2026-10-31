"""Call an MCP tool function directly (development / debugging without an MCP client).

  uv run python scripts/mcp_call.py omics get_de_results '{"gene": "CA9"}'
"""
import importlib
import json
import sys

SERVERS = {"omics": "ccrcc.mcp_servers.omics", "kg": "ccrcc.mcp_servers.kg", "external": "ccrcc.mcp_servers.external"}

if __name__ == "__main__":
    server, tool = sys.argv[1], sys.argv[2]
    args = json.loads(sys.argv[3]) if len(sys.argv) > 3 else {}
    fn = getattr(importlib.import_module(SERVERS[server]), tool)
    print(json.dumps(fn(**args), indent=1, ensure_ascii=False, default=str))
