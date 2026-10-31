# Linux / macOS shortcuts (Windows: use the uv commands directly and scripts/run_informe.ps1).
.PHONY: all test informe eval baselines report

all:            ## rebuild every deterministic artifact from scratch and render the report
	uv run ccrcc run-all

test:
	uv run pytest -q

informe:        ## agent step (Claude Code subscription login; no API key)
	./scripts/run_informe.sh

eval:           ## recall, outcome gold sets, temporal holdout, sensitivity, baselines -> evals/results.md
	uv run ccrcc eval

baselines:      ## closed-book LLM runs + agent-vs-ranking + replicate comparison -> evals/baselines/
	uv run ccrcc baselines

report:
	uv run python scripts/render.py
