# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

HgtSIM (Horizontal Gene Transfer Simulator) simulates horizontal gene transfer (HGT) in microbial communities. It mutates input gene sequences at configurable levels and inserts them into recipient genomes.

Publication: Song W, Steensen K, Thomas T. (2017) PeerJ 5:e4015 https://doi.org/10.7717/peerj.4015

## Installation & Dependencies

```bash
pip3 install HgtSIM
# or from source:
pip3 install .
```

**Python dependencies:** None (FASTA/GenBank I/O and codon tables are in `hgtsim/utils.py`).

**External tools:** BLAST+ (blastn, blastp, makeblastdb) must be on PATH or specified via `-blastn`/`-blastp` flags.

## Running

```bash
# Fixed mutation level (e.g. 10%)
HgtSIM -t genes.fasta -d distribution.txt -f input_genomes/ -r 1-0-1-1 -x fna -i 10

# Mixed/random mutation level (e.g. 5-25%)
HgtSIM -t genes.fasta -d distribution.txt -f input_genomes/ -r 1-0-1-1 -x fna -mixed 5-25

# With flanking sequences and keep_cds mode
HgtSIM -t genes.fasta -d distribution.txt -f input_genomes/ -i 10 -r 1-0-1-1 -x fna \
  -lf lf.fasta -rf rf.fasta -keep_cds -a annotation_gbk_folder/ -l 10
```

Test datasets are in `Test_datasets_Assessment_scripts_and_Commands/HgtSIM_working_directory/`.

## Architecture

The main tool is a Python 3 script **`bin/HgtSIM`** (~640 lines), with shared utilities in **`hgtsim/utils.py`** (FASTA/GenBank I/O, translation, codon tables).

**Two-step workflow:**
1. **Mutation simulation** — Reads input gene FASTA, applies codon-level mutations based on the specified level and ratio (`-r same_sense-non_same_sense-2bp-3bp`), then validates with BLAST (nucleotide and protein identity).
2. **Gene insertion** — Reads a distribution file mapping genes to recipient genomes, randomly inserts mutated genes into genome contigs (optionally with flanking sequences, optionally restricted to intergenic regions via `-keep_cds`).

**Key functions in `bin/HgtSIM`:**
- `get_codon_differences()` / `split_sequence()` / `get_synonymous_codons()` — codon mutation logic
- `get_mutant_codon_number()` — distributes mutations across the 4 mutation types based on ratio
- `get_random_insertion()` — core logic for inserting transfers into recipient genome sequences
- `get_flanking_seqs_dict()` — parses dynamic flanking sequence FASTA files

**Assessment scripts** in `Test_datasets_Assessment_scripts_and_Commands/Assessment_scripts/` validate recovery of inserted transfers from assembled reads (used for the paper's benchmarking sections).

## Testing

```bash
pytest tests/
```

Unit tests cover FASTA I/O, translation, codon tables, GenBank parsing, and reverse complement. Integration test runs a full HgtSIM simulation (requires BLAST+ on PATH).

## Important Notes

- Uses `os.system()` for BLAST calls — commands are constructed from user-provided paths.
- Non-deterministic by design (random mutations and insertion positions, no seed option).
