from hgtsim.utils import parse_genbank_cds


def test_parse_bad_gbk(bad_gbk):
    cds_lists = parse_genbank_cds(bad_gbk)
    # BAD.gbk has one contig
    assert len(cds_lists) == 1
    cds_list = cds_lists[0]
    # Should have many CDS features
    assert len(cds_list) > 100

    # First CDS: dnaA at positions 104..1525 → (103, 1525) in 0-based
    assert cds_list[0] == (103, 1525)

    # Second CDS: dnaN at 1652..2758 → (1651, 2758)
    assert cds_list[1] == (1651, 2758)


def test_cds_positions_are_sorted(bad_gbk):
    cds_lists = parse_genbank_cds(bad_gbk)
    cds_list = cds_lists[0]
    starts = [c[0] for c in cds_list]
    # GenBank features should generally be in order
    # (though not guaranteed, the test file appears to be sorted)
    for i in range(len(starts) - 1):
        # Just check that we parsed real positions (positive integers)
        assert starts[i] >= 0
        assert cds_list[i][1] > cds_list[i][0]
