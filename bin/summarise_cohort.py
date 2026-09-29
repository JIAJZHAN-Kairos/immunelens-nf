#!/usr/bin/env python3
"""Gather all manifest samples and all four loci without dropping failed estimates."""
import collections
import csv
import json
from pathlib import Path
import shutil
import sys


def summarise(manifest, results):
    with open(manifest, newline='') as handle:
        samples = [row['sample'].strip() for row in csv.DictReader(handle)]
    loci = ('TCRA', 'TCRB', 'TCRG', 'IGH')
    gathered = {}
    for file in Path(results).glob('*.immunelens/estimates.tsv'):
        with file.open(newline='') as handle:
            for row in csv.DictReader(handle, delimiter='\t'):
                key = (row['sample'], row['locus'])
                if key in gathered:
                    raise ValueError(f'Duplicate estimate: {key}')
                gathered[key] = row
    expected = {(sample, locus) for sample in samples for locus in loci}
    if set(gathered) != expected:
        raise ValueError(f'Incomplete cohort: missing={sorted(expected - set(gathered))[:10]}, '
                         f'unexpected={sorted(set(gathered) - expected)[:10]}')
    rows = [gathered[(sample, locus)] for sample in samples for locus in loci]
    with open('cohort_estimates.tsv', 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    with open('cohort_estimates_wide.tsv', 'w', newline='') as handle:
        fields = ['sample', 'genome', 'correction']
        for locus in loci:
            fields += [locus + '_fraction', locus + '_status']
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        for sample in samples:
            row = {'sample': sample, 'genome': gathered[(sample, 'TCRA')]['genome'], 'correction': 'unadjusted'}
            for locus in loci:
                result = gathered[(sample, locus)]
                row[locus + '_fraction'] = result['cell_fraction']
                row[locus + '_status'] = result['status']
            writer.writerow(row)
    qc = {'manifest_samples': len(samples), 'sample_locus_rows': len(rows),
          'expected_sample_locus_rows': len(samples) * len(loci),
          'status_counts': dict(collections.Counter(row['status'] for row in rows)),
          'genome_counts': dict(collections.Counter(gathered[(sample, 'TCRA')]['genome'] for sample in samples)),
          'correction': 'unadjusted', 'stub_run': any(row['status'] == 'stub' for row in rows),
          'note': 'Fractions are DNA-based estimates without purity/local-CN or IGH haplotype/CNA correction.'}
    Path('cohort_qc.json').write_text(json.dumps(qc, indent=2) + '\n')
    shutil.copyfile(manifest, 'samplesheet.csv')
    print(json.dumps(qc, indent=2))


if __name__ == '__main__':
    summarise(sys.argv[1], sys.argv[2])
