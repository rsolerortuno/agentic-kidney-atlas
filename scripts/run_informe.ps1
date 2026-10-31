# Run the orchestrated multi-agent analysis headless with Claude Code (Pro subscription login, no API key).
# The project MCP servers and the tools the agents need are passed explicitly, so no settings file is required.
$ErrorActionPreference = "Stop"
Set-Location (Split-Path -Parent $PSScriptRoot)
# Reload PATH from the registry so tools installed after this shell started (e.g. uv via winget) are found;
# the MCP servers in .mcp.json and the orchestrator's scripts are launched with `uv`.
$env:Path = [Environment]::GetEnvironmentVariable("Path", "User") + ";" + [Environment]::GetEnvironmentVariable("Path", "Machine") + ";" + $env:Path
if (-not (Get-Command uv -ErrorAction SilentlyContinue)) { throw "uv not found on PATH" }
$claude = (Get-Command claude -ErrorAction SilentlyContinue).Source
if (-not $claude) {
  $claude = Get-ChildItem "$env:LOCALAPPDATA\Claude-Data\*\claude-code\*\*\claude.exe" -ErrorAction SilentlyContinue |
            Sort-Object LastWriteTime -Descending | Select-Object -First 1 -ExpandProperty FullName
}
if (-not $claude) { throw "Claude Code CLI not found" }
# On Windows Claude Code runs shell commands through the PowerShell tool, so allow both shells.
$allowed = "mcp__omics-mcp,mcp__kg-mcp,mcp__external-mcp,mcp__cell-mcp,Bash(uv run python scripts/*),PowerShell(uv run python scripts/*),Read,Write,Edit,Agent,Task"
& $claude -p "/informe" --model claude-opus-5-5 --mcp-config .mcp.json --allowedTools $allowed
