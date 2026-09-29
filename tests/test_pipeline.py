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
prepare = load_module('prepare_manifest')


class ManifestTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.work = Path(self.temp.name)
        self.manifest = self.work / 'manifest.csv'
        self.manifest.write_text('sample,bam\nB,s3://bucket/B.bam\nA,s3://bucket/A.bam\n')
        self.output = self.work / 'immunelens.csv'
        for sample in ('A', 'B'):
            (self.work / (sample + '.purple.purity.tsv')).write_text('purity\tstatus\n0.5\tNORMAL\n')
            self.write_cn(sample)

    def tearDown(self):
        self.temp.cleanup()

    def write_cn(self, sample, heterogeneous=False):
        text = 'chromosome\tstart\tend\tgene\tisCanonical\tminCopyNumber\tmaxCopyNumber\n'
        for symbol, chromosome, start, end in prepare.PROXIES.values():
            high = 5 if heterogeneous and symbol == 'OR10G3' else 3
            text += f'chr{chromosome}\t{start}\t{end}\t{symbol}\ttrue\t3\t{high}\n'
        (self.work / (sample + '.purple.cnv.gene.tsv')).write_text(text)

    def test_proxy_cn_and_purity_preserve_source_and_sample_order(self):
        original = self.manifest.read_bytes()
        qc = prepare.prepare(self.manifest, self.work, self.work, self.output)
        with self.output.open() as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual([r['sample'] for r in rows], ['B', 'A'])
        self.assertEqual(self.manifest.read_bytes(), original)
        self.assertEqual(qc['samples_with_all_four_cn'], 2)
        self.assertTrue(all(r['purity'] == '0.5' and r['IGH_cn'] == '3' for r in rows))

    def test_heterogeneous_proxy_is_na_without_dropping_the_sample(self):
        self.write_cn('B', heterogeneous=True)
        qc = prepare.prepare(self.manifest, self.work, self.work, self.output)
        with self.output.open() as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[0]['TCRA_cn'], 'NA')
        self.assertEqual(rows[0]['TCRA_cn_status'], 'proxy_gene_heterogeneous')
        self.assertEqual(rows[0]['TCRB_cn'], '3')
        self.assertEqual(qc['samples_with_all_four_cn'], 1)

    def test_no_tumor_cannot_be_silently_reintroduced(self):
        (self.work / 'B.purple.purity.tsv').write_text('purity\tstatus\n0.08\tNO_TUMOR\n')
        with self.assertRaisesRegex(ValueError, 'NO_TUMOR'):
            prepare.prepare(self.manifest, self.work, self.work, self.output)
        self.assertFalse(self.output.exists())


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

    def test_mixed_corrections_keep_raw_and_adjusted_columns(self):
        path = self.results / 'A.immunelens' / 'estimates.tsv'
        with path.open() as handle:
            rows = list(csv.DictReader(handle, delimiter='\t'))
        for row in rows:
            row.update(raw_cell_fraction='NA', adjusted_cell_fraction='NA', purity='0.5',
                       local_cn='NA', cn_status='proxy_gene_heterogeneous', high_cell_fraction_flag='NA')
        rows[0].update(cell_fraction='0.25', raw_cell_fraction='0.2', adjusted_cell_fraction='0.25',
                       correction='purity_local_cn', local_cn='3', cn_status='proxy_gene_uniform', status='ok')
        with path.open('w') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t')
            writer.writeheader()
            writer.writerows(rows)
        gather.summarise(self.manifest, self.results)
        with open('cohort_estimates_wide.tsv') as handle:
            wide = list(csv.DictReader(handle, delimiter='\t'))
        self.assertEqual(wide[1]['correction'], 'mixed')
        self.assertEqual(wide[1]['TCRA_raw_fraction'], '0.2')
        self.assertEqual(wide[1]['TCRA_adjusted_fraction'], '0.25')
        self.assertEqual(wide[1]['TCRB_adjusted_fraction'], 'NA')
        self.assertEqual(json.loads(Path('cohort_qc.json').read_text())['correction_counts'],
                         {'unadjusted': 7, 'purity_local_cn': 1})


if __name__ == '__main__':
    unittest.main()
