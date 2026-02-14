import os
import pytest

TESTDATA_DIR = os.path.join(
    os.path.dirname(__file__), '..',
    'Test_datasets_Assessment_scripts_and_Commands', 'HgtSIM_working_directory'
)


@pytest.fixture
def testdata_dir():
    return TESTDATA_DIR


@pytest.fixture
def gene_transfers_fasta(testdata_dir):
    return os.path.join(testdata_dir, 'sequences_of_gene_transfers.fasta')


@pytest.fixture
def distribution_file(testdata_dir):
    return os.path.join(testdata_dir, 'distribution_of_transfers.txt')


@pytest.fixture
def recipient_genomes_dir(testdata_dir):
    return os.path.join(testdata_dir, 'selected_10_Betaproteobacteria')


@pytest.fixture
def gbk_dir(testdata_dir):
    return os.path.join(testdata_dir, 'selected_10_Betaproteobacteria_gbk')


@pytest.fixture
def bad_gbk(gbk_dir):
    return os.path.join(gbk_dir, 'BAD.gbk')


@pytest.fixture
def reference_output_dir(testdata_dir):
    return os.path.join(testdata_dir, 'outputs_10_1-0-1-1')
