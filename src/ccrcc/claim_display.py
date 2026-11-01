"""Compact presentation of exact numeric annotations; source details remain accessible."""
from html import escape
from .semantic import numeric_claim_text


def statement_html(statement: dict, evidence: dict) -> str:
    text = statement.get('statement') or statement.get('text') or statement.get('rationale') or ''
    spans = []
    for claim in statement.get('numeric_claims', []):
        record = evidence.get(claim.get('evidence_id'))
        if record is None:
            continue
        try:
            annotation = numeric_claim_text(claim, record)
        except (ValueError, TypeError):
            continue
        if text.count(annotation) != 1:
            continue
        start = text.index(annotation)
        context = ' · '.join(str(record[k]) for k in ('gene','layer','profile') if record.get(k))
        metric = '.'.join(map(str,claim['field']))
        value = claim['display'] + ('%' if claim['unit'] in ('percent','native_percent') else '')
        label = f'{context} · {metric} = {value}' if context else f'{metric} = {value}'
        spans.append((start, start+len(annotation), f'<span class="numeric-claim" tabindex="0" title="{escape(annotation,quote=True)}">{escape(label)}</span>'))
    out, end = [], 0
    for start, stop, html in sorted(set(spans)):
        if start < end:
            continue
        out.extend((escape(text[end:start]),html)); end=stop
    out.append(escape(text[end:]))
    return ''.join(out)
