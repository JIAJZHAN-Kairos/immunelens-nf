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
    for row in rows:
        for field, value in {'raw_cell_fraction': row['cell_fraction'],
                             'adjusted_cell_fraction': 'NA', 'purity': 'NA', 'local_cn': 'NA',
                             'cn_status': 'not_supplied', 'high_cell_fraction_flag': 'NA',
                             'igh_correction': 'not_provided' if row['locus'] == 'IGH' else 'not_applicable',
                             'igh_uncorrected_fraction': 'NA', 'igh_germline_fraction': 'NA',
                             'igh_combined_fraction': 'NA', 'igh_somatic_qc': 'NA',
                             'igh_somatic_selected': 'NA'}.items():
            row.setdefault(field, value)
    with open('cohort_estimates.tsv', 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    with open('cohort_estimates_wide.tsv', 'w', newline='') as handle:
        fields = ['sample', 'genome', 'purity', 'correction']
        igh_fields = ['igh_correction', 'igh_uncorrected_fraction', 'igh_germline_fraction',
                      'igh_combined_fraction', 'igh_somatic_qc', 'igh_somatic_selected']
        fields += igh_fields
        for locus in loci:
            fields += [locus + suffix for suffix in ('_fraction', '_raw_fraction', '_adjusted_fraction',
                       '_status', '_correction', '_local_cn', '_cn_status', '_high_cell_fraction_flag')]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        for sample in samples:
            corrections = {gathered[(sample, locus)]['correction'] for locus in loci}
            row = {'sample': sample, 'genome': gathered[(sample, 'TCRA')]['genome'],
                   'purity': gathered[(sample, 'TCRA')]['purity'],
                   'correction': next(iter(corrections)) if len(corrections) == 1 else 'mixed'}
            row.update({field: gathered[(sample, 'IGH')][field] for field in igh_fields})
            for locus in loci:
                result = gathered[(sample, locus)]
                row[locus + '_fraction'] = result['cell_fraction']
                row[locus + '_raw_fraction'] = result['raw_cell_fraction']
                row[locus + '_adjusted_fraction'] = result['adjusted_cell_fraction']
                row[locus + '_status'] = result['status']
                row[locus + '_correction'] = result['correction']
                row[locus + '_local_cn'] = result['local_cn']
                row[locus + '_cn_status'] = result['cn_status']
                row[locus + '_high_cell_fraction_flag'] = result['high_cell_fraction_flag']
            writer.writerow(row)
    correction_counts = dict(collections.Counter(row['correction'] for row in rows))
    qc = {'manifest_samples': len(samples), 'sample_locus_rows': len(rows),
          'expected_sample_locus_rows': len(samples) * len(loci),
          'status_counts': dict(collections.Counter(row['status'] for row in rows)),
          'genome_counts': dict(collections.Counter(gathered[(sample, 'TCRA')]['genome'] for sample in samples)),
          'correction': next(iter(correction_counts)) if len(correction_counts) == 1 else 'mixed',
          'correction_counts': correction_counts,
          'high_cell_fraction_flags': sum(row['high_cell_fraction_flag'] == 'TRUE' for row in rows),
          'stub_run': any(row['status'] == 'stub' for row in rows),
          'igh_correction_counts': dict(collections.Counter(gathered[(sample, 'IGH')]['igh_correction'] for sample in samples)),
          'note': 'cell_fraction uses purity/local-CN adjustment where supplied. IGH raw_cell_fraction is after the selected locus CNV correction and before purity adjustment; uncorrected and candidate IGH fractions are separate. Inspect IGH correction state and QC.'}
    Path('cohort_qc.json').write_text(json.dumps(qc, indent=2) + '\n')
    shutil.copyfile(manifest, 'samplesheet.csv')
    print(json.dumps(qc, indent=2))


if __name__ == '__main__':
    summarise(sys.argv[1], sys.argv[2])
