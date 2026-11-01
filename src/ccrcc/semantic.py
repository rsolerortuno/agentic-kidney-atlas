"""Semantic verification of statements against the evidence they cite (deterministic, no LLM).

Existence and gene linkage (verify.py) prove that a statement points at real evidence about the right gene. They do
not prove that the statement says what the evidence says. This module adds two checks:

1. Numbers. Match a specific metric in one unambiguous cited record, with its sign,
   scale and printed precision. Structured numeric_claims support nested metrics
   with explicit record context. Ranking/global number pools never supply support.
2. Direction. If a gene statement uses only "up" words (or only "down" words) and every significant omics_de
   evidence it cites for that gene points the other way, the statement contradicts its evidence.
"""
from __future__ import annotations

import math
import re


NUM = re.compile(r"(?<![\w.])([-−])?(\d+(?:\.\d+)?(?:[eE][-−+]?\d+)?)(\s?%)?(?!\w)")
# Words after which a number is a label, not a measurement. Whole words only: single letters ("p", "q") would exempt
# p- and q-values and any word ending in that letter ("from 0.5"). Numbers glued to letters ("M13", "CD70", "3p") are
# never extracted by NUM in the first place.
EXEMPT_BEFORE = re.compile(r"(\b(phase|rank|ranked|top|grade|stage|tier|version|step|chromosome|chr|pmid:?|uniprot|chembl)"
                           r"|#|[A-Za-z]-)\s*$", re.I)
THRESHOLD_BEFORE = re.compile(r"(\bat|<|≤|<=|>|≥|>=|threshold|cut-?off|alpha|α|fdr|bh)\s*$", re.I)
THRESHOLDS = {0.001, 0.01, 0.05, 0.1, 0.25}   # significance levels written as such ("not significant at 0.05")
DRUG_AFTER = re.compile(r"^\s+[A-Z][A-Z0-9-]{3,}")   # "IODINE I 124 GIRENTUXIMAB": part of a drug name
UP = re.compile(r"\b(up|upregulated|up-regulated|increased|elevated|higher|overexpress\w*|induced|gained)\b", re.I)
DOWN = re.compile(r"\b(down|downregulated|down-regulated|decreased|reduced|lower|depleted|suppressed)\b", re.I)


def numbers_in(text: str, *, with_spans: bool = False) -> list[dict]:
    out = []
    for m in NUM.finditer(text.replace("−", "-")):
        raw = m.group(2)
        before = text[max(0, m.start() - 12):m.start()]
        val = float(raw)
        decimals = len(raw.split(".")[1].split("e")[0].split("E")[0]) if "." in raw else 0
        sci = "e" in raw.lower()
        is_int = "." not in raw and not sci
        if is_int and 1900 <= val <= 2100 and not metric_field(before, text[m.end():m.end() + 60]):
            continue
        if EXEMPT_BEFORE.search(before) or DRUG_AFTER.match(text[m.end():m.end() + 20]):
            continue
        if val in THRESHOLDS and THRESHOLD_BEFORE.search(before):
            continue
        out.append({"text": m.group(0).strip(), "value": -val if m.group(1) else val, "decimals": decimals,
                    "sci": sci, "percent": bool(m.group(3))})
        if with_spans:
            out[-1].update(start=m.start(), end=m.end())
    return out


def _matches(n: dict, candidates) -> bool:
    v = n["value"]
    if not math.isfinite(v):
        return False
    for c in candidates:
        if not math.isfinite(c):
            continue
        target = c * 100 if n["percent"] else c
        if n["sci"]:
            # Scientific notation: tolerance follows the precision actually printed.
            mantissa = n["text"].lower().split("e")[0].lstrip("+-")
            digits = len(mantissa.split(".")[1]) if "." in mantissa else 0
            exponent = int(n["text"].lower().split("e")[1].strip().rstrip("%"))
            tol = 0.5 * 10 ** (exponent - digits)
        else:
            tol = 0.5 * 10 ** (-n["decimals"])
        if abs(v - target) <= tol + abs(target) * 1e-12:
            return True

    return False


# Only these source fields store fractions that can be displayed as percentages.
# A log fold-change, p-value or score must never silently become a percentage.
FRACTION_FIELDS = {"share_of_tumor_transcripts", "freq", "loss_freq", "gain_freq",
                   "p_top10", "p_top25", "p_top100", "cpm_share", "source_share"}
PERCENT_FIELDS = {"pct_expressing", "source_pct",
                  "tumour_cells_expressing_pct_weakest_quarter_of_donors",
                  "tumour_cells_expressing_pct_median_donor"}  # already stored on the 0–100 scale
CONTEXT_FIELDS = ("gene", "type", "source", "layer", "cohort", "cell_type", "celltype", "cell_group", "tissue", "arm", "contrast", "profile", "track")


def numeric_claim_text(claim: dict, evidence: dict) -> str:
    """Canonical, auditable annotation for a number, including its source context.

    Example claim: {evidence_id, field: ['arms', 'atezo', 'pvalue'],
                    display: '0.032', unit: 'native'}.
    The exact returned annotation must occur once in the statement. No free-form
    label can relabel a protein result as RNA, or a p-value as an effect size.
    """
    path = claim.get("field")
    if not isinstance(path, list) or not path or any(not isinstance(k, (str, int)) or isinstance(k, bool) for k in path):
        raise ValueError("field must be a nonempty path of keys/list indices")
    value = evidence
    try:
        for key in path:
            if isinstance(value, list) and (not isinstance(key, int) or key < 0):
                raise ValueError("invalid list index")
            value = value[key]
    except (KeyError, IndexError, TypeError) as exc:
        raise ValueError("metric path not found") from exc
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError("metric is not a finite number")
    display = claim.get("display", "")
    if not isinstance(display, str) or not NUM.fullmatch(display):
        raise ValueError("display must be one numeric literal")
    unit = claim.get("unit")
    if unit not in ("native", "percent", "native_percent") or "%" in display:
        raise ValueError("unit must be native, percent or native_percent; display must omit %")
    if unit == "percent" and path[-1] not in FRACTION_FIELDS:
        raise ValueError("metric is not registered as a fraction")
    if unit == "native_percent" and path[-1] not in PERCENT_FIELDS:
        raise ValueError("metric is not registered as a native percentage")
    # Prefix makes all magnitudes, including small integers and years, measurable.
    tokens = numbers_in("log2FC " + display + ("%" if unit == "percent" else ""))
    if len(tokens) != 1 or not _matches(tokens[0], [value]):
        raise ValueError("display does not match metric sign/scale/precision")
    context = "; ".join(f"{k}={evidence[k]}" for k in CONTEXT_FIELDS if evidence.get(k) is not None)
    return f"[{context}; {'.'.join(map(str, path))}={display}{'%' if unit in ('percent', 'native_percent') else ''}]"


def check(text: str, gene: str | None, evidence: list[dict], extra_values: set[float] | None = None,
          numeric_claims: list[dict] | None = None) -> dict:
    """Return unmatched numbers and direction contradictions for one statement."""
    # Never accept a number merely because it occurs elsewhere in the evidence,
    # a ranking row, a date, an identifier or a global summary.
    unmatched, binding_errors = [], []
    remaining = text
    indexed = {e.get("evidence_id"): e for e in evidence if e.get("evidence_id")}
    if numeric_claims is not None and not isinstance(numeric_claims, list):
        binding_errors.append("numeric_claims must be a list")
        numeric_claims = []
    for i, claim in enumerate(numeric_claims or []):
        try:
            if not isinstance(claim, dict):
                raise ValueError("claim must be an object")
            ev = indexed.get(claim.get("evidence_id"))
            if ev is None or (gene and ev.get("gene", "").upper() != gene.upper()):
                raise ValueError("claim must reference cited evidence for the same gene")
            annotation = numeric_claim_text(claim, ev)
            if remaining.count(annotation) != 1:
                raise ValueError("canonical metric/context annotation must occur exactly once")
            remaining = remaining.replace(annotation, " " * len(annotation), 1)
        except (ValueError, TypeError) as exc:
            binding_errors.append(f"numeric_claims[{i}]: {exc}")
    for n in numbers_in(remaining, with_spans=True):
        start, end = n["start"], n["end"]
        before, after = remaining[max(0, start - 80):start], remaining[end:end + 60]
        field = metric_field(before, after)
        candidates = []
        context_mismatch = False
        layers = re.findall(r"\b(RNA|protein|proteome|phospho|phosphoproteome)\b", before, re.I)
        named_layer = {"proteome":"protein", "phosphoproteome":"phospho"}.get(layers[-1].lower(), layers[-1].lower()) if layers else None
        for e in evidence:
            if gene and (e.get("gene") or "").upper() != gene.upper():
                continue
            if field and isinstance(e.get(field), (int, float)) and not isinstance(e.get(field), bool):
                candidates.append(e[field])
                if named_layer and e.get("layer") and e["layer"] != named_layer:
                    context_mismatch = True
        if context_mismatch or len(candidates) != 1 or (n["percent"] and field not in FRACTION_FIELDS) or not _matches(n, candidates):
            unmatched.append(n["text"])
    contradiction = None
    if gene:
        de = [e for e in evidence if e.get("type") == "omics_de" and e.get("significant")
              and (e.get("gene") or "").upper() == gene.upper()]
        signs = {1 if e["log2fc"] > 0 else -1 for e in de if e.get("log2fc") is not None}
        # A survival statement such as "higher expression predicts better survival"
        # is not a claim that expression increased in tumour versus normal tissue.
        # Keep this a narrow DE wording check, not a bag of direction words.
        for sentence in re.split(r"(?<=[.!?])\s+", text):
            if re.search(r"\b(survival|prognos\w*|hazard|dependency|rank\w*|score\w*)\b", sentence, re.I):
                continue
            up, down = bool(UP.search(sentence)), bool(DOWN.search(sentence))
            if len(signs) == 1 and up != down:
                said = 1 if up else -1
                if said != next(iter(signs)):
                    contradiction = f"text says {'up' if up else 'down'}, cited omics evidence says " \
                                    f"{'up' if next(iter(signs)) > 0 else 'down'}"
                    break
    return {"unmatched_numbers": unmatched, "direction_contradiction": contradiction,
            "numeric_binding_errors": binding_errors}


# Narrow grammar for legacy prose. Complex/nested comparisons must be emitted as
# structured claims rather than silently checked against a bag of numbers.
def metric_field(before: str, after: str) -> str | None:
    patterns = [
        (r"log2\s*FC\s*[=:]?\s*$", "log2fc"),
        (r"(?:padj|q(?:-value)?)\s*[=:]?\s*$", "padj"),
        (r"(?:pvalue|p-value|\bp)\s*[=:]?\s*$", "pvalue"),
        (r"(?:Cox\s+)?log\s*HR\s*[=:]?\s*$", "log_hr_per_sd"),
    ]
    for pattern, field in patterns:
        if re.search(pattern, before, re.I):
            return field
    if re.match(r"\s+pairs\b", after, re.I):
        return "n_pairs"
    if re.match(r"\s+of (?:its |the )?tumou?r transcripts", after, re.I):
        return "share_of_tumor_transcripts"
    return None
