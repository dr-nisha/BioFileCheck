import gzip
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from biofilecheck.core import validate

class ValidationTests(unittest.TestCase):
    def check_text(self, text, fmt='fasta', zipped=False):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / ('input.' + fmt + ('.gz' if zipped else ''))
            if zipped:
                with gzip.open(path, 'wt') as handle:
                    handle.write(text)
            else:
                path.write_text(text)
            return validate(path)

    def codes(self, report):
        return {i['code'] for i in report['issues']}

    def test_multiline_gc_and_ambiguity(self):
        r = self.check_text('>one\nacg\ntNN\n>two\nGC\n')
        self.assertTrue(r['valid'])
        self.assertEqual((r['record_count'], r['total_bases']), (2, 8))
        self.assertEqual(r['ambiguous_bases'], 2)
        self.assertAlmostEqual(r['gc_percent_canonical'], 66.6667)

    def test_all_ambiguous_gc_is_null(self):
        self.assertIsNone(self.check_text('>one\nNNRY\n')['gc_percent_canonical'])

    def test_duplicate_first_token_empty_and_invalid(self):
        r = self.check_text('>one desc\nACX\n>one other\n')
        self.assertEqual(self.codes(r), {'DUPLICATE_ID', 'EMPTY_SEQUENCE', 'INVALID_BASE'})

    def test_missing_header_and_empty_id(self):
        self.assertEqual(self.codes(self.check_text('ACGT\n>\nACGT\n')), {'MISSING_HEADER', 'EMPTY_ID'})

    def test_empty_input(self):
        self.assertIn('NO_RECORDS', self.codes(self.check_text('')))

    def test_gzip_crlf(self):
        self.assertTrue(self.check_text('>one\r\nACGT\r\n', zipped=True)['valid'])

    def test_fastq_quality_can_start_with_at(self):
        r = self.check_text('@one\nACGT\n+one\n@@@@\n', 'fastq')
        self.assertTrue(r['valid'])

    def test_fastq_length_separator_and_header(self):
        r = self.check_text('one\nACGT\nx\nIII\n', 'fastq')
        self.assertEqual(self.codes(r), {'INVALID_HEADER','INVALID_SEPARATOR','QUALITY_LENGTH'})

    def test_fastq_repeated_header(self):
        self.assertIn('HEADER_MISMATCH', self.codes(self.check_text('@a\nAC\n+b\nII\n', 'fastq')))

    def test_fastq_truncation(self):
        self.assertIn('TRUNCATED_RECORD', self.codes(self.check_text('@a\nAC\n+\n', 'fastq')))

    def test_quality_character(self):
        self.assertIn('QUALITY_CHARACTER', self.codes(self.check_text('@a\nAC\n+\nI \n', 'fastq')))

    def test_unknown_extension(self):
        with self.assertRaises(ValueError):
            validate('file.xyz')

    def test_cli_exit_codes_and_json(self):
        root = Path(__file__).resolve().parents[1]
        env = dict(os.environ, PYTHONPATH=str(root/'src'))
        for name, expected in [('valid.fasta',0), ('invalid.fastq',1), ('missing.fasta',2)]:
            result = subprocess.run([sys.executable,'-m','biofilecheck.cli',str(root/'examples'/name)],env=env,capture_output=True,text=True)
            self.assertEqual(result.returncode, expected, result.stderr)
            if expected != 2:
                self.assertEqual(json.loads(result.stdout)['valid'], expected == 0)

    def test_cli_protects_input(self):
        root = Path(__file__).resolve().parents[1]
        env = dict(os.environ, PYTHONPATH=str(root/'src'))
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)/'test.fa'
            p.write_text('>one\nAC\n')
            result = subprocess.run([sys.executable,'-m','biofilecheck.cli',str(p),'--output',str(p)],env=env,capture_output=True)
            self.assertEqual(result.returncode,2)
            self.assertEqual(p.read_text(),'>one\nAC\n')

if __name__ == '__main__':
    unittest.main()
