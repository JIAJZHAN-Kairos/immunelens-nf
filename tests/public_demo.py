"""Run all four loci using the upstream publication's public coverage examples."""
import concurrent.futures
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import urllib.request
import csv
import math

root = Path(__file__).resolve().parents[1]
destination = Path(sys.argv[1]).resolve()
destination.mkdir(parents=True, exist_ok=True)
coverage = destination / 'DEMO.coverage'
coverage.mkdir(exist_ok=True)
files = {
    'TCRA': ('test_cov_TCRA.txt', '58141c038a2fd174dce839c6151c6956'),
    'TCRB': ('test_cov_TCRB.txt', '7dd431fb44a55aa0412a0a24d4f8ec69'),
    'TCRG': ('test_cov_TCRG.txt', 'bfd0cca7cb67273099ab1f8c065c5bce'),
    'IGH': ('tumour_test_IGH.txt', 'a96a8c92dfffef9bbaf33b5e60600d0e'),
    'NORMAL': ('germline_test_IGH.txt', 'f4e4349d191b88b6fb617448b421c278'),
}


def fetch(item):
    locus, (filename, checksum) = item
    data = urllib.request.urlopen('https://zenodo.org/records/11094087/files/' + filename + '?download=1', timeout=120).read()
    if hashlib.md5(data).hexdigest() != checksum:
        raise ValueError(filename + ': published checksum mismatch')
    with gzip.open(coverage / (locus + '.txt.gz'), 'wb') as handle:
        handle.write(data)
    return locus, data.count(b'\n')


with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
    counts = dict(pool.map(fetch, files.items()))
normal = destination / 'DEMO.normal.coverage'
normal.mkdir(exist_ok=True)
(coverage / 'NORMAL.txt.gz').replace(normal / 'IGH.txt.gz')
(normal / 'sample.json').write_text(json.dumps({'genome': 'hg38', 'provided': True,
    'coverage_positions': {'IGH': counts.pop('NORMAL')}}))
(coverage / 'sample.json').write_text(json.dumps({'sample': 'DEMO', 'genome': 'hg38',
    'coverage_positions': counts, 'correction': 'unadjusted', 'source': 'Zenodo 11094087 upstream tutorial'}))
corrections = destination / 'corrections.json'
corrections.write_text(json.dumps({'genome': 'hg38', 'purity': 0.9,
    'TCRA_cn': 3, 'TCRB_cn': 3, 'TCRG_cn': None, 'IGH_cn': 3,
    'TCRG_cn_status': 'proxy_gene_heterogeneous'}))
subprocess.run([str(root / 'bin/run_immunelens.R'), 'DEMO', str(coverage), str(corrections), str(normal)],
               cwd=destination, check=True)
with (destination / 'DEMO.immunelens/estimates.tsv').open() as handle:
    rows = list(csv.DictReader(handle, delimiter='\t'))
assert len(rows) == 4
assert all(row['status'] in ('ok', 'high_cell_fraction') for row in rows), rows
assert all(0 <= float(row['cell_fraction']) <= 1 for row in rows), rows
for row in rows:
    raw = float(row['raw_cell_fraction'])
    if row['locus'] == 'TCRG':
        assert row['correction'] == 'unadjusted'
        assert row['adjusted_cell_fraction'] == 'NA'
        assert float(row['cell_fraction']) == raw
        assert row['cn_status'] == 'proxy_gene_heterogeneous'
    else:
        assert row['correction'] == 'purity_local_cn'
        assert math.isclose(float(row['adjusted_cell_fraction']), raw * 1.45, abs_tol=1e-10)
        assert row['cell_fraction'] == row['adjusted_cell_fraction']
        assert (row['high_cell_fraction_flag'] == 'TRUE') == (raw * 1.45 > 0.1)
igh = next(row for row in rows if row['locus'] == 'IGH')
assert igh['igh_correction'] in ('germline_only', 'germline_somatic'), igh
germline = float(igh['igh_germline_fraction'])
if igh['igh_correction'] == 'germline_somatic':
    assert igh['igh_somatic_qc'] == 'TRUE'
    assert float(igh['igh_combined_fraction']) < germline
    assert float(igh['raw_cell_fraction']) == float(igh['igh_combined_fraction'])
else:
    assert float(igh['raw_cell_fraction']) == germline
assert not math.isclose(float(igh['igh_uncorrected_fraction']), float(igh['raw_cell_fraction']), abs_tol=1e-6)
assert (destination / 'DEMO.immunelens/IGH.germline_regions.tsv').is_file()
assert (destination / 'DEMO.immunelens/IGH.somatic_regions.tsv').is_file()
subprocess.run(['Rscript', str(root / 'tests/igh_correction_qc.R'), str(root / 'bin/correct_igh.R'), str(destination)],
               cwd=destination, check=True)
# A supplied mismatched genome must fail before fitting any models.
bad_normal = destination / 'bad.normal.coverage'
bad_normal.mkdir(exist_ok=True)
(bad_normal / 'sample.json').write_text(json.dumps({'provided': True, 'genome': 'hg19'}))
bad = subprocess.run([str(root / 'bin/run_immunelens.R'), 'BAD', str(coverage), str(corrections), str(bad_normal)],
                     cwd=destination, capture_output=True, text=True)
assert bad.returncode != 0 and 'Matched normal and tumour BAM genomes differ' in bad.stderr, bad.stderr
print(json.dumps(rows, indent=2))
