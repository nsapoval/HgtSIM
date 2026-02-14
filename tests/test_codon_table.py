import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from hgtsim.utils import FORWARD_TABLE, START_CODONS, SYNONYMOUS_CODONS


def test_forward_table_size():
    assert len(FORWARD_TABLE) == 61


def test_forward_table_no_stops():
    for codon, aa in FORWARD_TABLE.items():
        assert aa != '*'


def test_start_codons():
    assert 'ATG' in START_CODONS
    assert 'TTG' in START_CODONS
    assert 'CTG' in START_CODONS


def test_synonymous_codons_cover_all_amino_acids():
    all_aas = set(FORWARD_TABLE.values())
    assert set(SYNONYMOUS_CODONS.keys()) == all_aas


def test_synonymous_codons_lists_correct():
    # All codons in a synonymous group should code for the same amino acid
    for aa, codons in SYNONYMOUS_CODONS.items():
        for codon in codons:
            assert FORWARD_TABLE[codon] == aa


def test_get_synonymous_codons():
    # Import the function from bin/HgtSIM
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'bin'))

    # Inline test: TGG has no synonymous codons (only codon for W)
    from hgtsim.utils import SYNONYMOUS_CODONS
    w_codons = SYNONYMOUS_CODONS['W']
    assert w_codons == ['TGG']

    # ATG is the only Met codon
    m_codons = SYNONYMOUS_CODONS['M']
    assert m_codons == ['ATG']

    # Leucine has 6 codons
    l_codons = SYNONYMOUS_CODONS['L']
    assert len(l_codons) == 6
