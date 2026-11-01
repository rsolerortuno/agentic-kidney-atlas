import hashlib
import pytest
from ccrcc.literature_review import validate


def test_reviewed_abstract_text_and_quote_are_both_checked():
    text='An assay reported a specific value in a defined context.'
    entry={'pmid':1,'abstract_sha256':hashlib.sha256(text.encode()).hexdigest(),
           'source_quotes':['a specific value']}
    validate(entry,text)
    with pytest.raises(ValueError,match='Abstract changed'):
        validate(entry,text.replace('specific','different'))
    with pytest.raises(ValueError,match='quotation'):
        validate({**entry,'source_quotes':['invented result']},text)
