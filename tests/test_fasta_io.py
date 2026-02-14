import os
import tempfile
from hgtsim.utils import parse_fasta, read_fasta, write_fasta


def test_parse_fasta_multi_record(gene_transfers_fasta):
    records = list(parse_fasta(gene_transfers_fasta))
    assert len(records) == 100
    assert records[0].id == 'AAM_00175'
    assert records[0].description == 'Electron transfer flavoprotein subunit alpha'
    assert records[0].seq.startswith('ATGGCAGAG')


def test_parse_fasta_sequence_no_whitespace(gene_transfers_fasta):
    for rec in parse_fasta(gene_transfers_fasta):
        assert '\n' not in rec.seq
        assert ' ' not in rec.seq


def test_read_fasta_single(recipient_genomes_dir):
    fasta_path = os.path.join(recipient_genomes_dir, 'BAD.fna')
    rec = read_fasta(fasta_path)
    assert rec is not None
    assert rec.id == 'BAD'
    assert len(rec.seq) > 0


def test_write_fasta_roundtrip(tmp_path):
    out_file = str(tmp_path / 'test.fasta')
    with open(out_file, 'w') as fh:
        write_fasta(fh, 'seq1', 'a description', 'ATGCATGCATGC' * 10)
        write_fasta(fh, 'seq2', '', 'GGGG')

    records = list(parse_fasta(out_file))
    assert len(records) == 2
    assert records[0].id == 'seq1'
    assert records[0].description == 'a description'
    assert records[0].seq == 'ATGCATGCATGC' * 10
    assert records[1].id == 'seq2'
    assert records[1].description == ''
    assert records[1].seq == 'GGGG'


def test_write_fasta_wrapping(tmp_path):
    out_file = str(tmp_path / 'wrap.fasta')
    seq = 'A' * 150
    with open(out_file, 'w') as fh:
        write_fasta(fh, 'test', '', seq, wrap=60)

    with open(out_file) as fh:
        lines = fh.readlines()
    # header + 3 sequence lines (60+60+30)
    assert lines[0].startswith('>')
    assert len(lines[1].rstrip()) == 60
    assert len(lines[2].rstrip()) == 60
    assert len(lines[3].rstrip()) == 30


def test_write_fasta_no_wrap(tmp_path):
    out_file = str(tmp_path / 'nowrap.fasta')
    seq = 'A' * 200
    with open(out_file, 'w') as fh:
        write_fasta(fh, 'test', '', seq, wrap=0)

    with open(out_file) as fh:
        lines = fh.readlines()
    assert len(lines) == 2
    assert lines[1].rstrip() == seq
