"""Streaming checks for DNA FASTA and strict four-line FASTQ files."""
import gzip
from pathlib import Path

DNA = frozenset("ACGTRYSWKMBDHVN")

def validate(path, file_format=None):
    """Return a JSON-compatible report; coordinates are one-based lines.

    GC percentage is calculated over canonical A/C/G/T bases only. Ambiguous
    IUPAC DNA symbols are allowed and counted separately. No inference about
    biological validity or sequencing quality is made.
    """
    path = Path(path)
    suffixes = path.name.lower().removesuffix(".gz")
    if file_format is None:
        ext = Path(suffixes).suffix
        file_format = {".fa":"fasta", ".fasta":"fasta", ".fna":"fasta",
                       ".fq":"fastq", ".fastq":"fastq"}.get(ext)
    if file_format not in ("fasta", "fastq"):
        raise ValueError("Specify --format fasta or fastq for an unrecognized extension")
    issues, seen = [], set()
    lengths = []
    counts = dict.fromkeys("ACGT", 0)
    ambiguous = 0
    def issue(line, code, message):
        issues.append({"line": line, "code": code, "message": message})
    def record(header, sequence, line):
        nonlocal ambiguous
        identifier = header.split()[0] if header.split() else ""
        if not identifier:
            issue(line, "EMPTY_ID", "Record identifier is empty")
        elif identifier in seen:
            issue(line, "DUPLICATE_ID", "Record identifier is repeated: " + identifier)
        seen.add(identifier)
        sequence = sequence.upper()
        lengths.append(len(sequence))
        if not sequence:
            issue(line, "EMPTY_SEQUENCE", "Record has no sequence")
        invalid = sorted(set(sequence) - DNA)
        if invalid:
            issue(line, "INVALID_BASE", "Non-IUPAC DNA symbols: " + repr("".join(invalid)))
        for base in sequence:
            if base in counts:
                counts[base] += 1
            elif base in DNA:
                ambiguous += 1
    opener = gzip.open if path.suffix.lower() == ".gz" else open
    with opener(path, "rt", encoding="utf-8-sig") as handle:
        if file_format == "fasta":
            header, chunks, start = None, [], 0
            for number, raw in enumerate(handle, 1):
                line = raw.rstrip("\r\n")
                if line.startswith(">"):
                    if header is not None:
                        record(header, "".join(chunks), start)
                    header, chunks, start = line[1:], [], number
                elif not line:
                    continue
                elif header is None:
                    issue(number, "MISSING_HEADER", "Sequence appears before a FASTA header")
                else:
                    chunks.append(line)
            if header is not None:
                record(header, "".join(chunks), start)
        else:
            number = 0
            while True:
                raw = handle.readline()
                if not raw:
                    break
                start = number + 1
                block = [raw] + [handle.readline() for _ in range(3)]
                number += 4
                if any(part == "" for part in block):
                    issue(start, "TRUNCATED_RECORD", "Expected four lines per FASTQ record")
                    break
                header, sequence, plus, quality = [s.rstrip("\r\n") for s in block]
                if not header.startswith("@"):
                    issue(start, "INVALID_HEADER", "FASTQ header must start with @")
                if not plus.startswith("+"):
                    issue(start+2, "INVALID_SEPARATOR", "FASTQ separator must start with +")
                elif plus[1:] and plus[1:] != header[1:]:
                    issue(start+2, "HEADER_MISMATCH", "Repeated separator header differs from record header")
                if len(sequence) != len(quality):
                    issue(start+3, "QUALITY_LENGTH", "Sequence and quality lengths differ")
                if any(not 33 <= ord(c) <= 126 for c in quality):
                    issue(start+3, "QUALITY_CHARACTER", "Quality must use printable Phred+33-range characters")
                record(header[1:] if header.startswith("@") else header, sequence, start)
    if not lengths:
        issue(1, "NO_RECORDS", "No complete sequence records found")
    canonical = sum(counts.values())
    return {"schema_version":"1.0", "format":file_format,
            "valid":not issues, "record_count":len(lengths),
            "total_bases":sum(lengths), "min_length":min(lengths, default=0),
            "max_length":max(lengths, default=0), "canonical_bases":canonical,
            "ambiguous_bases":ambiguous,
            "gc_percent_canonical":round(100*(counts["G"]+counts["C"])/canonical, 4) if canonical else None,
            "issues":issues}
