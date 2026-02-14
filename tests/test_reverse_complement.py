from hgtsim.utils import reverse_complement


def test_basic():
    assert reverse_complement('ATGC') == 'GCAT'


def test_palindrome():
    assert reverse_complement('AATT') == 'AATT'


def test_single_base():
    assert reverse_complement('A') == 'T'
    assert reverse_complement('C') == 'G'


def test_lowercase():
    assert reverse_complement('atgc') == 'gcat'


def test_empty():
    assert reverse_complement('') == ''


def test_longer_sequence():
    seq = 'ATGCATGCATGC'
    rc = reverse_complement(seq)
    # Reverse complement of reverse complement should be original
    assert reverse_complement(rc) == seq
