from ccrcc.claim_display import statement_html
from ccrcc.semantic import numeric_claim_text


def test_compact_annotation_preserves_context_and_escapes_source_and_prose():
    record={'evidence_id':'EV:1','gene':'CA9','source':'<script>bad()</script>', 'layer':'rna','log2fc':2.3}
    claim={'evidence_id':'EV:1','field':['log2fc'],'display':'2.3','unit':'native'}
    text='<img src=x onerror=bad()> '+numeric_claim_text(claim,record)
    html=statement_html({'text':text,'numeric_claims':[claim]}, {'EV:1':record})
    assert '<img' not in html and '<script>' not in html
    assert 'CA9 · rna · log2fc = 2.3' in html
    assert 'source=&lt;script&gt;' in html
    assert 'tabindex="0"' in html


def test_unverifiable_annotations_are_not_silently_rewritten():
    assert statement_html({'text':'unverified <value>','numeric_claims':[]},{}) == 'unverified &lt;value&gt;'
