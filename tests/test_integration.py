import os
import shutil
import subprocess
import pytest


def blast_available():
    try:
        result = subprocess.run(
            ['blastn', '-version'],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=5
        )
        return result.returncode == 0
    except (FileNotFoundError, subprocess.TimeoutExpired, OSError):
        return False


@pytest.mark.skipif(not blast_available(), reason='BLAST+ not available')
def test_hgtsim_end_to_end(testdata_dir, tmp_path):
    """Run HgtSIM with test data at mutation level 10 and verify outputs."""
    bin_dir = os.path.join(os.path.dirname(__file__), '..', 'bin')
    hgtsim_script = os.path.join(bin_dir, 'HgtSIM')

    gene_transfers = os.path.join(testdata_dir, 'sequences_of_gene_transfers.fasta')
    distribution = os.path.join(testdata_dir, 'distribution_of_transfers.txt')
    genomes = os.path.join(testdata_dir, 'selected_10_Betaproteobacteria')

    # Ensure hgtsim package is importable by the subprocess
    project_root = os.path.join(os.path.dirname(__file__), '..')
    env = os.environ.copy()
    env['PYTHONPATH'] = os.path.abspath(project_root) + os.pathsep + env.get('PYTHONPATH', '')

    # Run from tmp_path so output goes there
    result = subprocess.run(
        [
            'python3', hgtsim_script,
            '-t', gene_transfers,
            '-d', distribution,
            '-f', genomes,
            '-r', '1-0-1-1',
            '-x', 'fna',
            '-i', '10',
            '-p', 'HgtSIM',
            '-quiet',
        ],
        cwd=str(tmp_path),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        env=env,
        timeout=300,
    )

    assert result.returncode == 0, (
        'HgtSIM failed:\nstdout: %s\nstderr: %s' % (result.stdout, result.stderr)
    )

    output_dir = tmp_path / 'HgtSIM_outputs_10_1-0-1-1'
    assert output_dir.is_dir()

    # Check key output files exist
    assert (output_dir / 'input_sequence_mutant_nc.fasta').is_file()
    assert (output_dir / 'input_sequence_mutant_aa.fasta').is_file()
    assert (output_dir / 'input_sequence_aa.fasta').is_file()
    assert (output_dir / 'Step_1_mutation_report.txt').is_file()
    assert (output_dir / 'Step_2_insertion_report.txt').is_file()

    genomes_dir = output_dir / 'Genomes_with_transfers'
    assert genomes_dir.is_dir()

    # Should have 10 output genome files
    genome_files = list(genomes_dir.glob('*.fna'))
    assert len(genome_files) == 10

    # Verify mutant FASTA has 100 sequences (same as input)
    from hgtsim.utils import parse_fasta
    mutant_records = list(parse_fasta(str(output_dir / 'input_sequence_mutant_nc.fasta')))
    assert len(mutant_records) == 100

    aa_records = list(parse_fasta(str(output_dir / 'input_sequence_mutant_aa.fasta')))
    assert len(aa_records) == 100
