"""Command-line interface and stable exit codes."""
import argparse
import json
import sys
from pathlib import Path
from .core import validate

def main():
    parser = argparse.ArgumentParser(description="Validate local DNA FASTA or four-line FASTQ files")
    parser.add_argument("input", type=Path)
    parser.add_argument("--format", choices=("fasta", "fastq"))
    parser.add_argument("--output", type=Path, help="Write JSON report instead of stdout")
    args = parser.parse_args()
    if args.output and args.output.resolve() == args.input.resolve():
        parser.error("Output must differ from input")
    try:
        report = validate(args.input, args.format)
        content = json.dumps(report, indent=2) + "\n"
        if args.output:
            args.output.write_text(content, encoding="utf-8")
        else:
            print(content, end="")
    except (OSError, ValueError, EOFError) as exc:
        print("biofilecheck: " + str(exc), file=sys.stderr)
        return 2
    return 0 if report["valid"] else 1

if __name__ == "__main__":
    raise SystemExit(main())
