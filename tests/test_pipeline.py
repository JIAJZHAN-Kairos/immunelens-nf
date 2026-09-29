import csv
import gzip
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


def load_module(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / 'bin' / (name + '.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


extract = load_module('extract_coverage')
gather = load_module('summarise_cohort')


class GenomeTests(unittest.TestCase):
    def test_hg38_with_chr_prefix(self):
        self.assertEqual(extract.detect_genome('@SQ\tSN:chr7\tLN:159345973\n@SQ\tSN:chr14\tLN:107043718'), ('hg38', 'chr'))

    def test_hg19_without_chr_prefix(self):
        self.assertEqual(extract.detect_genome('@SQ\tSN:7\tLN:159138663\n@SQ\tSN:14\tLN:107349540'), ('hg19', ''))

    def test_mixed_build_is_rejected(self):
        with self.assertRaises(RuntimeError):
            extract.detect_genome('@SQ\tSN:7\tLN:159345973\n@SQ\tSN:14\tLN:107349540')

    def test_missing_contig_is_rejected(self):
        with self.assertRaises(RuntimeError):
            extract.detect_genome('@SQ\tSN:7\tLN:159345973')


@unittest.skipUnless(shutil.which('samtools'), 'samtools required for indexed coverage test')
class CoverageTests(unittest.TestCase):
    def test_base_and_mapping_quality_match_upstream(self):
        with tempfile.TemporaryDirectory() as directory:
            work = Path(directory)
            sam = work / 'input.sam'
            sam.write_text('@HD\tVN:1.6\tSO:coordinate\n'
                           '@SQ\tSN:chr7\tLN:159345973\n@SQ\tSN:chr14\tLN:107043718\n'
                           'good\t0\tchr7\t101\t60\t20M\t*\t0\t0\t' + 'A' * 20 + '\t' + 'I' * 20 + '\n'
                           'low_mapq\t0\tchr7\t101\t19\t20M\t*\t0\t0\t' + 'A' * 20 + '\t' + 'I' * 20 + '\n'
                           'low_baseq\t0\tchr7\t101\t60\t20M\t*\t0\t0\t' + 'A' * 20 + '\t' + '!' * 20 + '\n')
            bam = work / 'input.bam'
            subprocess.run(['samtools', 'view', '-b', '-o', str(bam), str(sam)], check=True)
            # Test the alternate file.bai index name, not just file.bam.bai.
            subprocess.run(['samtools', 'index', str(bam), str(work / 'input.bai')], check=True)
            regions = work / 'regions.tsv'
            regions.write_text('genome\tlocus\tchromosome\tstart\tend\n'
                               'hg38\tTCRA\t14\t101\t140\n'
                               'hg38\tTCRB\t7\t101\t140\n'
                               'hg38\tTCRG\t7\t101\t140\n'
                               'hg38\tIGH\t14\t101\t140\n')
            subprocess.run(['python3', str(ROOT / 'bin/extract_coverage.py'), '--sample', 'S',
                            '--bam', str(bam), '--regions', str(regions)], cwd=work, check=True)
            output = work / 'S.coverage'
            with gzip.open(output / 'TCRB.txt.gz', 'rt') as handle:
                rows = [line.strip().split('\t') for line in handle]
            self.assertEqual(len(rows), 20)
            self.assertEqual({row[2] for row in rows}, {'1'})
            self.assertEqual([int(row[1]) for row in rows], list(range(101, 121)))
            metadata = json.loads((output / 'sample.json').read_text())
            self.assertEqual(metadata['coverage_positions']['TCRA'], 0)
            self.assertEqual(metadata['genome'], 'hg38')


class CohortTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.work = Path(self.temp.name)
        self.old_cwd = Path.cwd()
        os.chdir(self.work)
        self.manifest = self.work / 'manifest.csv'
        self.manifest.write_text('sample,bam\nB,s3://bucket/B.bam\nA,s3://bucket/A.bam\n')
        self.results = self.work / 'collected'
        for sample in ('A', 'B'):
            directory = self.results / (sample + '.immunelens')
            directory.mkdir(parents=True)
            with (directory / 'estimates.tsv').open('w', newline='') as handle:
                fields = ['sample', 'locus', 'genome', 'status', 'cell_fraction', 'correction', 'message']
                writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t')
                writer.writeheader()
                for locus in ('TCRA', 'TCRB', 'TCRG', 'IGH'):
                    writer.writerow(dict(zip(fields, [sample, locus, 'hg38', 'no_coverage', 'NA', 'unadjusted', 'No coverage'])))

    def tearDown(self):
        os.chdir(self.old_cwd)
        self.temp.cleanup()

    def test_all_na_estimates_are_retained_in_manifest_order(self):
        gather.summarise(self.manifest, self.results)
        with open('cohort_estimates.tsv') as handle:
            rows = list(csv.DictReader(handle, delimiter='\t'))
        self.assertEqual(len(rows), 8)
        self.assertEqual([row['sample'] for row in rows], ['B'] * 4 + ['A'] * 4)
        self.assertEqual(Path('samplesheet.csv').read_bytes(), self.manifest.read_bytes())
        self.assertEqual(json.loads(Path('cohort_qc.json').read_text())['status_counts'], {'no_coverage': 8})

    def test_missing_sample_cannot_pass_gather(self):
        shutil.rmtree(self.results / 'A.immunelens')
        with self.assertRaisesRegex(ValueError, 'Incomplete cohort'):
            gather.summarise(self.manifest, self.results)

    def test_duplicate_locus_cannot_pass_gather(self):
        path = self.results / 'A.immunelens' / 'estimates.tsv'
        with path.open('a') as handle:
            handle.write(path.read_text().splitlines()[1] + '\n')
        with self.assertRaisesRegex(ValueError, 'Duplicate estimate'):
            gather.summarise(self.manifest, self.results)


if __name__ == '__main__':
    unittest.main()
