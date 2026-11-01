# Repository cleanup audit

The cleanup follows the scientific/code revision. It removes only artifacts with
no active consumer; age alone is not a deletion criterion.

Removed:

- `runs/dev_v1_agents/`: seven exploratory agent/verification JSON files. Current
  assembly and baseline readers select dated runs and explicit replicates; no
  active code or documentation links to this exploratory directory.
- `IMPLEMENTATION_STATUS.md`: transient implementation notes; useful conclusions
  are preserved in the audit documents.
- `reports/preview_revision.html`: temporary preview superseded by the canonical
  `reports/ccrcc_report.html`.

Recoverable copies are retained outside the repository in the audit workspace.
`.gitignore` excludes this development directory, temporary report previews and
local implementation notes. Cleanup does not touch the original user checkout.

Retained intentionally:

- All frozen source inputs, reviewed saved model answers and the manifest.
- Replicate findings consumed by `baselines.run`, and their associated run records.
- M12 archived findings used in `evals/history.md`, and historical logs linked by
  `PROGRESS.md`. These are audit evidence, not unused scratch files.
- The eight script entry points: CLI/report validation, agent workflow wrappers,
  cross-platform launchers and the documented direct MCP-call utility.
- Vendored viewer code and license, required for the standalone offline report.

The full final reconstruction completed successfully. All 43 compared tables and
the recalculated evaluation matched the reviewed results. The cleaned repository
retains the canonical report and the reproducibility evidence; temporary work and
recovery copies remain outside it.
