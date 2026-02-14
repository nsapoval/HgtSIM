from hgtsim.utils import translate, CODON_TABLE


def test_translate_known_codons():
    assert translate('ATG') == 'M'
    assert translate('TAA') == '*'
    assert translate('TAG') == '*'
    assert translate('TGA') == '*'
    assert translate('TTT') == 'F'
    assert translate('GGG') == 'G'


def test_translate_full_sequence():
    # ATG=M, TTT=F, GGC=G, TAA=*
    result = translate('ATGTTTGGCTAA')
    assert result == 'MFG*'


def test_translate_ignores_trailing_partial_codon():
    result = translate('ATGTT')
    assert result == 'M'


def test_translate_empty():
    assert translate('') == ''


def test_translate_preserves_case():
    # Our translate uppercases internally
    assert translate('atg') == 'M'
    assert translate('Atg') == 'M'


def test_codon_table_complete():
    assert len(CODON_TABLE) == 64
    # All 4-letter combinations present
    bases = 'ACGT'
    for a in bases:
        for b in bases:
            for c in bases:
                assert a + b + c in CODON_TABLE
