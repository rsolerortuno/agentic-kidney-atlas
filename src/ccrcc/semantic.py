"""Semantic verification of statements against the evidence they cite (deterministic, no LLM).

Existence and gene linkage (verify.py) prove that a statement points at real evidence about the right gene. They do
not prove that the statement says what the evidence says. This module adds two checks:

1. Numbers. Match a specific metric in one unambiguous cited record, with its sign,
   scale and printed precision. Structured numeric_claims support nested metrics
   with explicit record context. Ranking/global number pools never supply support.
2. Direction. If a gene statement uses only "up" words (or only "down" words) and every significant omics_de
   evidence it cites for that gene points the other way, the statement contradicts its evidence. Statements without
   a gene (summary, sections) are checked sentence by sentence for sentences that name exactly one cited gene.
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
UP = re.compile(r"\b(up|upregulated|up-regulated|increased|elevated|higher|overexpress\w*|induced|gained|raised|rises?|rising)\b", re.I)
DOWN = re.compile(r"\b(down|downregulated|down-regulated|decreased|reduced|lower|lowered|depleted|suppressed|underexpress\w*|diminished)\b", re.I)
# A negated sentence ("is not higher") states the opposite of its direction word; it is skipped, not read literally.
NEGATION = re.compile(r"\b(not|no|never|neither|nor|without|n't)\b|n't\b", re.I)
# ordering words: "placed lower", "listed higher" describe a position in a list, not expression
DIRECTION_EXEMPT = re.compile(r"\b(survival|prognos\w*|hazard|dependency|rank\w*|score\w*|placed|listed|ordered|"
                              r"positioned|priorit\w*)\b", re.I)


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
            if ev is None or (gene and (ev.get("gene") or "").upper() != gene.upper()):
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
    return {"unmatched_numbers": unmatched, "direction_contradiction": _direction(text, gene, evidence),
            "numeric_binding_errors": binding_errors}


def _de_signs(evidence: list[dict]) -> dict[str, dict[str, set[int]]]:
    """gene -> layer -> signs of the significant omics_de changes cited (layer "*" pools all layers)."""
    signs: dict[str, dict[str, set[int]]] = {}
    for e in evidence:
        if e.get("type") == "omics_de" and e.get("significant") and e.get("log2fc") is not None and e.get("gene"):
            sign = 1 if e["log2fc"] > 0 else -1
            g = signs.setdefault(e["gene"].upper(), {})
            g.setdefault("*", set()).add(sign)
            g.setdefault(str(e.get("layer") or "*").lower(), set()).add(sign)
    return signs


# Clauses are checked one at a time, so "protein down ... but RNA up" or a survival clause next to an expression clause
# no longer hides the whole sentence (verifier sabotage test, 2026-11-03: 16 % of flipped directions were caught).
CLAUSE = re.compile(r"[;,]|\bbut\b|\bwhile\b|\bwhereas\b|\band\b|\bwith\b", re.I)
LAYER_WORD = re.compile(r"\b(RNA|mRNA|transcripts?|protein|proteome|phospho\w*)\b", re.I)
PRINTED_LFC = re.compile(r"log2\s*FC\s*[=:]?\s*([-−]?\d+(?:\.\d+)?)", re.I)
RELATIVE_TO_NORMAL = re.compile(r"\b(higher|lower|up|down|increased|decreased|elevated|reduced)\s+in\s+(the\s+)?"
                                r"(normal|healthy|adjacent)", re.I)
NORMAL_ONLY = re.compile(r"\b(normal|healthy|adjacent)\b", re.I)
TUMOUR = re.compile(r"\btumou?rs?\b", re.I)


def _layer_of(clause: str) -> str | None:
    m = LAYER_WORD.findall(clause)
    if not m:
        return None
    w = m[-1].lower()
    return "rna" if w in ("rna", "mrna", "transcript", "transcripts") else "protein" if w.startswith("prot") else "phospho"


def _direction(text: str, gene: str | None, evidence: list[dict]) -> str | None:
    """Direction words that contradict the cited significant omics_de evidence.

    With a statement gene every sentence is read as being about it. Without one (summary and section paragraphs) a
    sentence is checked only when it names exactly one gene that has cited omics_de evidence, so multi-gene
    comparisons are never attributed to the wrong gene.
    """
    signs = _de_signs(evidence)
    for sentence in re.split(r"(?<=[.!?])\s+", text):
        if gene:
            target = gene.upper()
        else:
            named = [g for g in signs if re.search(rf"(?<![\w-]){re.escape(g)}(?![\w-])", sentence, re.I)]
            if len(named) != 1:
                continue
            target = named[0]
        where = "" if gene else f" ({target})"
        # canonical annotations use ";" between fields: keep each annotation inside its clause
        flat = re.sub(r"\[[^\]]*\]", lambda m: m.group(0).replace(";", " ").replace(",", " "), sentence)
        for clause in CLAUSE.split(flat):
            # A survival statement such as "higher expression predicts better survival"
            # is not a claim that expression increased in tumour versus normal tissue.
            # Keep this a narrow DE wording check, not a bag of direction words.
            if DIRECTION_EXEMPT.search(clause) or NEGATION.search(clause):
                continue
            # in a survival / ranking sentence only clauses that name the tumour, a layer or a log2FC are about
            # tumour-versus-normal expression ("higher expression is associated with worse survival" is not)
            if DIRECTION_EXEMPT.search(sentence) and not (TUMOUR.search(clause) or LAYER_WORD.search(clause)
                                                          or PRINTED_LFC.search(clause)):
                continue
            up, down = bool(UP.search(clause)), bool(DOWN.search(clause))
            if up == down:
                continue
            said = 1 if up else -1
            if RELATIVE_TO_NORMAL.search(clause):
                said = -said                      # "lower in normal kidney" = higher in tumour
            elif NORMAL_ONLY.search(clause) and not TUMOUR.search(clause):
                continue
            word = "up" if said > 0 else "down"
            m = PRINTED_LFC.search(clause)
            if m and float(m.group(1).replace("−", "-")) != 0 and (float(m.group(1).replace("−", "-")) > 0) != (said > 0):
                return f"text says {word}{where} next to a printed log2FC of the opposite sign ({m.group(1)})"
            layer = _layer_of(clause)
            s = signs.get(target, {})
            seen_set = s.get(layer) if layer and layer in s else s.get("*", set())
            if len(seen_set) == 1 and said != next(iter(seen_set)):
                seen = next(iter(seen_set))
                return f"text says {word}{where}, cited omics evidence says {'up' if seen > 0 else 'down'}"
    return None


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
