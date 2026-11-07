# BioFileCheck

**A local DNA FASTA/FASTQ validation toolkit for reproducible bioinformatics workflows.**

Status: **v0.1.0 working prototype**. This repository was created with AI coding assistance; it is not an established research package or a clinically validated tool.

## Why this project?
Malformed sequence files can interrupt analysis or lead to confusing downstream errors. BioFileCheck provides explicit structural checks and machine-readable results before a file enters an analysis pipeline. It separates file validity from biological interpretation.

## Implemented
- Multiline DNA FASTA and strict four-line FASTQ; plain text and gzip.
- IUPAC DNA alphabet checks, empty sequences, duplicate identifiers.
- FASTQ structure, sequence/quality length, optional repeated-header consistency, and printable quality-character checks.
- Sequence counts, lengths, ambiguous-base counts, canonical-base GC percentage.
- JSON reports, stable exit codes, synthetic fixtures, regression tests, and GitHub Actions configuration.
- Python standard library only at runtime; files stay on the local machine.

## Quick start
Requires Python 3.10 or newer. From this repository:

```bash
python -m venv .venv
# Linux/macOS:
source .venv/bin/activate
# Windows PowerShell instead: .venv\Scripts\Activate.ps1
python -m pip install -e .
biofilecheck examples/valid.fasta
biofilecheck examples/invalid.fastq --output report.json
python -m unittest discover -s tests -v
```

Without installation, Linux/macOS users can run:
```bash
PYTHONPATH=src python -m biofilecheck.cli examples/valid.fasta
```

Exit codes: `0` valid; `1` validation findings; `2` input/output or configuration error. Use `--format fasta` or `--format fastq` for unrecognized file extensions.

## Scientific conventions and limits
FASTA identifiers use the first whitespace-delimited token. Lowercase bases are accepted. Allowed DNA symbols are ACGTRYSWKMBDHVN; U and alignment gaps are rejected. Blank FASTA lines are ignored, but embedded sequence whitespace is flagged. GC percentage uses only A/C/G/T as its denominator; all-ambiguous input returns JSON null.

FASTQ support is deliberately limited to four-line records. Multiline FASTQ, protein FASTA, legacy quality encodings, paired-end consistency, VCF, and clinical interpretation are not supported. Printable ASCII quality checks do not establish the true quality encoding or read quality. Reports retain findings and identifier sets in memory, so very large or highly malformed inputs can consume substantial memory. A `.gz` suffix selects gzip decoding.

See [example report](docs/example-report.json), [development roadmap](docs/ROADMAP.md), and [application notes](docs/APPLICATION.md).

## Data and provenance
All bundled sequences are synthetic demonstration data, not patient or research-subject data. No legacy company source code or BioMystery evaluation assets were imported. Do not commit credentials, private biological datasets or third-party proprietary code.

## License
No open-source license is granted in this initial version. Rights are reserved pending the repository owner's licensing decision. Public visibility alone does not grant a reuse license.
