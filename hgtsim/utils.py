"""Lightweight replacements for Biopython functionality used by HgtSIM."""

from collections import namedtuple


FastaRecord = namedtuple('FastaRecord', ['id', 'description', 'seq'])


def parse_fasta(filepath):
    """Yield FastaRecord(id, description, seq) from a multi-FASTA file."""
    with open(filepath) as fh:
        header = None
        chunks = []
        for line in fh:
            line = line.rstrip('\n').rstrip('\r')
            if line.startswith('>'):
                if header is not None:
                    seq = ''.join(chunks)
                    parts = header.split(None, 1)
                    sid = parts[0]
                    desc = parts[1] if len(parts) > 1 else ''
                    yield FastaRecord(sid, desc, seq)
                header = line[1:]
                chunks = []
            else:
                chunks.append(line)
        if header is not None:
            seq = ''.join(chunks)
            parts = header.split(None, 1)
            sid = parts[0]
            desc = parts[1] if len(parts) > 1 else ''
            yield FastaRecord(sid, desc, seq)


def read_fasta(filepath):
    """Read a single-record FASTA file, return a FastaRecord."""
    for record in parse_fasta(filepath):
        return record


def write_fasta(handle, seq_id, description, sequence, wrap=60):
    """Write a single FASTA record to an open file handle."""
    if description:
        handle.write('>%s %s\n' % (seq_id, description))
    else:
        handle.write('>%s\n' % seq_id)
    if wrap and wrap > 0:
        for i in range(0, len(sequence), wrap):
            handle.write(sequence[i:i + wrap] + '\n')
    else:
        handle.write(sequence + '\n')


# Standard genetic code (NCBI translation table 1)
CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

# Non-stop codons (61 sense codons) mapping codon -> amino acid
FORWARD_TABLE = {k: v for k, v in CODON_TABLE.items() if v != '*'}

START_CODONS = ['TTG', 'CTG', 'ATG']

# Synonymous codons grouped by amino acid
SYNONYMOUS_CODONS = {}
for _codon, _aa in FORWARD_TABLE.items():
    SYNONYMOUS_CODONS.setdefault(_aa, []).append(_codon)


def translate(dna_sequence):
    """Translate a DNA sequence to protein using the standard genetic code."""
    protein = []
    for i in range(0, len(dna_sequence) - 2, 3):
        codon = dna_sequence[i:i + 3].upper()
        aa = CODON_TABLE.get(codon, 'X')
        protein.append(aa)
    return ''.join(protein)


_COMPLEMENT = str.maketrans('ACGTacgt', 'TGCAtgca')


def reverse_complement(dna_sequence):
    """Return the reverse complement of a DNA sequence."""
    return dna_sequence.translate(_COMPLEMENT)[::-1]


def parse_genbank_cds(filepath):
    """Parse a GenBank file and yield a list of (start_0based, end) CDS tuples per record.

    Each call to next() returns the CDS list for one LOCUS/record in the file.
    """
    import re
    cds_lists = []
    current_cds = []
    in_record = False

    with open(filepath) as fh:
        in_features = False
        for line in fh:
            if line.startswith('LOCUS'):
                if in_record:
                    cds_lists.append(current_cds)
                    current_cds = []
                in_record = True
                in_features = False
            elif line.startswith('FEATURES'):
                in_features = True
            elif line.startswith('ORIGIN') or line.startswith('CONTIG'):
                in_features = False
            elif in_features and line.startswith('     CDS'):
                # Extract location from the CDS feature line
                location_str = line.strip().split(None, 1)[1].strip()
                # Handle complement, join, etc. — extract the outer range
                # Remove complement() wrapper
                inner = location_str
                if inner.startswith('complement('):
                    inner = inner[len('complement('):-1]
                if inner.startswith('join('):
                    inner = inner[len('join('):-1]
                # Find all numeric ranges
                ranges = re.findall(r'(\d+)\.\.(\d+)', inner)
                if ranges:
                    # Use the overall span (min start, max end)
                    start = min(int(r[0]) for r in ranges) - 1  # convert to 0-based
                    end = max(int(r[1]) for r in ranges)
                    current_cds.append((start, end))
            elif line.startswith('//'):
                if in_record:
                    cds_lists.append(current_cds)
                    current_cds = []
                    in_record = False

    if in_record:
        cds_lists.append(current_cds)

    return cds_lists
