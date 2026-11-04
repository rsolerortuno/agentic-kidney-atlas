# One run of the pipeline architecture benchmark. -Arm multi (/informe with subagents) or single (/informe_single).
# Writes runs/<Run>/ through the normal scripts and Claude Code's JSON result to session_<Run>.json next to this file.
# -Model sets the orchestrator (or single agent) model; specialists keep the model in their own definition.
param([string]$Arm, [string]$Run, [string]$Model = "claude-opus-5-5")
Set-Location (Split-Path -Parent (Split-Path -Parent $PSScriptRoot))
$env:Path = [Environment]::GetEnvironmentVariable("Path", "User") + ";" + [Environment]::GetEnvironmentVariable("Path", "Machine") + ";" + $env:Path
$claude = (Get-Command claude -ErrorAction SilentlyContinue).Source
if (-not $claude) {
  $claude = Get-ChildItem "$env:LOCALAPPDATA\Claude-Data\*\claude-code\*\*\claude.exe" | Sort-Object LastWriteTime -Descending | Select-Object -First 1 -ExpandProperty FullName
}
$tools = "mcp__omics-mcp,mcp__kg-mcp,mcp__external-mcp,mcp__cell-mcp,Bash(uv run python scripts/*),PowerShell(uv run python scripts/*),Read,Write,Edit"
if ($Arm -eq "multi") { $cmd = ".claude/commands/informe.md"; $tools = $tools + ",Agent,Task" } else { $cmd = ".claude/commands/informe_single.md" }
$prompt = @"
Read ${cmd} and follow steps 1 to 5 exactly, with RUN = runs/${Run}: pass ${Run} as the date argument to
scripts/integrate.py, assemble.py and verify.py. Do not read, modify or delete any other folder under runs/.
Skip step 6 (do not render the report). This run is part of an architecture benchmark.
"@
$out = Join-Path $PSScriptRoot "session_${Run}.json"
& $claude -p $prompt --model $Model --mcp-config .mcp.json --allowedTools $tools --output-format json --max-budget-usd 40 > $out 2> (Join-Path $PSScriptRoot "session_${Run}.err")
"exit=$LASTEXITCODE"
