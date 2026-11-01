# Numeric claims and the limits of automatic verification

Passing the checks means **no failures detected by the automatic checks**. It is
not complete scientific validation. The checks establish evidence linkage,
specific numeric comparisons and a limited up/down consistency rule. They do not
prove causal interpretations, treatment efficacy or the truth of every sentence.

Do not support a number by pointing to a ranking row or an unrelated value that
happens to match. A finding's `evidence_ids` must include the actual source record.
Gene-specific findings must use evidence for that gene.

## Simple statements

Legacy prose accepts only a narrow metric grammar: `log2FC`, `padj`/`q`, `pvalue`/`p`,
`logHR`, a count followed by `pairs`, or a percentage followed by
`of its tumour transcripts`. The cited records must supply exactly one value for
the recognized field. Multiple records with that field are ambiguous, even if
their numbers agree. Use the structured form below when citing multiple layers.

Sign and printed precision are preserved. Percent scaling is allowed only for
registered fraction fields; fold changes and p-values cannot become percentages.
Small integer measurements are checked too. Years, named identifier labels and
explicit significance thresholds have limited lexical exemptions; those
exemptions are not verification of the associated claim.

## Structured metrics, including nested evidence

Every statement can carry a `numeric_claims` list:

```json
{
  "gene": "CA9",
  "statement": "Observed association [gene=CA9; type=treatment; source=trial1; cohort=validation; arms.atezo.log_hr=-0.31].",
  "evidence_ids": ["EV:example"],
  "numeric_claims": [
    {"evidence_id": "EV:example", "field": ["arms", "atezo", "log_hr"],
     "display": "-0.31", "unit": "native"}
  ]
}
```

This is an illustrative schema, not a biological result. Use
`ccrcc.semantic.numeric_claim_text(claim, evidence_record)` to construct the exact
annotation from a real cited record. It includes the gene, evidence type, source
and available layer/cohort/cell-type/arm/contrast metadata, then the complete
metric path. It must occur exactly once in the statement. Do not edit its context
labels or copy an annotation from another record. The original record and its
provenance remain available through the cited evidence ID.

`field` is a nonempty list of dictionary keys or nonnegative list indices.
The endpoint must be a finite numeric value, not a boolean or a number embedded
in arbitrary text. `display` is one signed numeric literal, with optional
scientific notation. `unit` is `native`, or `percent` for a registered fraction
field; omit `%` from `display`. Use `native_percent` for registered fields already
stored on the 0–100 scale, such as single-cell `pct_expressing`; this adds the
percent sign without multiplying again. Other transformations are rejected.

All numeric text outside these annotations is still checked by the legacy
grammar. Free-form interpretation surrounding an annotation remains subject to
scientific review; the annotation does not validate it.

## Review failures rather than hiding them

Read `verification_agents.json` and `verification_final.json`. Correct the
original agent/narrative inputs and re-run verification and assembly. A semantic
failure dropped during assembly or agent filtering remains a final failure until
the inputs are corrected and checked again. Do not shorten the report merely to
obtain a zero-failure counter, refresh stale provenance to silence a failure, or
describe archived runs as having passed a newer verifier without rechecking them.
