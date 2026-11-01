# Scientific and provenance review of saved findings

This records corrections to the saved analysis, not a claim of complete scientific
validation. The previous text remains available in Git history. The reviewed report has been regenerated. Its saved automatic verification
checks 130 statements with zero unsupported, semantic or stale-evidence failures;
this result does not replace scientific validation. The final clean reproduction completed successfully; see reproduction.md and
cleanup.md for execution evidence and the repository cleanup.

## Omics findings, run 2026-11-01

All twelve original genes are retained: CA9, CD70, NNMT, ENPP3, SCARB1, SLC16A3,
PSAT1, PCK1, ALDH6A1, FABP6, HK2 and CA4. The revised findings contain 139 explicit
metric bindings. Each binding identifies one cited record, one field, its source
and available context. Verification against the saved provenance basis found no
failures in these twelve findings.

Corrections made:

- RNA, protein and TCGA effects are separately identified instead of sharing an
  ambiguous set of numbers. Adjusted p-values use the precision actually stored
  in the cited graph records.
- Purity extrapolation is described as model-based; it cannot by itself establish
  a tumour-cell origin. Cell-attribution evidence must be assessed separately.
- Missing citations for purity, survival, dependency and copy-number measurements
  were added from their actual records.
- FABP6's previously unspecified correlation of 0.309 belongs to protein. RNA's
  correlation is 0.409. Both are now explicitly identified.
- Copy-number correlations support association, not proof that gain or loss
  causes expression change.
- Univariate survival associations are not presented as independently validated
  prognostic biomarkers. Nonsignificant survival results remain visible.
- Measured non-selective dependencies are preserved as evidence, rather than
  omitted from the graph. Their measured effects, selectivity and adjusted
  p-values are cited; absence of selectivity is not equated to absence of data.
- Computed ranks and scores have content-addressed evidence records identifying
  their profile and configuration. They cannot support another metric merely by
  matching its number, and they are not biological validation.
- Clinical development is described using the frozen source snapshot and is not
  treated as proof of efficacy.

The revised text is a reviewed correction of saved agent output. It is not a new
execution of the historical Claude agent models.

## Cell/structure and external evidence review

The nineteen cell-agent findings are retained and now carry 218 exact numeric
bindings. The twelve external-agent findings are retained, with the MTOR
comparison separated from the KDM5C paragraph into its own gene-tagged finding:
thirteen findings with 136 bindings. Omics, cell and external blocks currently
pass together: 44 findings, 493 numeric bindings, no failures detected including
the saved evidence-basis checks. This is not a complete scientific validation.

Additional corrections:

- `pct_expressing` and discovery `source_pct` already store percentages on the
  0ÃƒÂ¢Ã¢â€šÂ¬Ã¢â‚¬Å“100 scale. `native_percent` preserves that scale; fraction-to-percent
  conversion is rejected for those fields.
- The girentuximab statement had cited an ethoxzolamide record for CA9, and the
  cusatuzumab statement had cited a vorsetuzumab record for CD70. These now cite
  the actual compounds. Numeric/gene checks alone could not detect that error.
- The external block lists each frozen drug record with its own evidence ID and
  distinguishes target-level development annotations from ccRCC indications.
- Cell-attribution statements are scoped to the sampled cohorts. CD70's differing
  assignments across cohorts remain explicit. Purity extrapolation is not treated
  as proof of cell of origin.
- New graph cavity records use `structural_context`, an `exposed_loop` component
  and a `structural_context_score`. They do not label an exposed loop as a degron.
  Low scores cannot exclude molecular-glue feasibility. Legacy stored raw score
  keys remain readable, but the active MCP tool is `structural_context_ranking`.
- EPAS1's selective-dependency call depends on the configured adjusted-p threshold
  of 0.1; its adjusted p-value of 0.061 is above 0.05. That limitation is explicit.
- Low mutation counts remain measured evidence instead of disappearing from the
  graph. A numerically stored zero p-value is not described as a literal zero
  probability.
- The limited DE direction check no longer misreads a sentence about higher
  expression and survival as a tumour-versus-normal expression-change claim.

## Literature, applied findings and integrated narrative

All five blocks have now been reviewed: omics 12, literature 11, external 13,
cell/structure 19 and applied 22 findings (77 total). None was discarded by the
verifier. The integrated narrative retains 14 sections and 11 target entries.

The literature review corrected two substantive readings: PMID41971124 reports
PSAT1 downregulation, so pathway enrichment in that article is not evidence of an
opposite expression direction; its intervention experiments concern HPGD.
PMID37523719 does include clear-cell renal cancer cell-line experiments, contrary
to the earlier summary. Reviewed full-abstract evidence records preserve source
hashes and exact supporting excerpts. Reported inhibitor measurements remain
source-reported experimental results, not predictions made here.

Treatment-effect point estimates with nonsignificant tests are described as such.
No subtype difference is not evidence of a uniform pan-subtype biomarker. The
candidate evidence cards distinguish missing, non-favorable and contradictory
observations, and state the intended modality and its evidential limits.

## Reproduction and cleanup audit

See reproduction.md for the frozen input inventory, final exit-zero clean run,
exact comparison of 43 tables and identical recalculated evaluation. The run-order
bug found in the first attempt was corrected and verified from clean inputs.

The obsolete `runs/dev_v1_agents` outputs and temporary report preview have been
removed after checking their consumers. The transient implementation log was
moved outside the repository after preserving its useful facts in these audit
documents. See cleanup.md for retained evidence and recoverable backups. Source
data and frozen model inputs were preserved.
