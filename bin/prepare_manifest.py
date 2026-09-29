#!/usr/bin/env python3
"""Add PURPLE purity and the paper's nearby-gene CN proxies to a BAM manifest."""
import argparse
import collections
import csv
import hashlib
import json
import math
from pathlib import Path


PROXIES = {
    'TCRA': ('OR10G3', '14', 21568520, 21580076),
    'TCRB': ('PRSS58', '7', 142252143, 142258058),
    'TCRG': ('STARD3NL', '7', 38178245, 38230670),
    'IGH': ('TMEM121', '14', 105526583, 105530202),
}


def prepare(manifest, purity_dir, gene_cn_dir, output):
    manifest, output = Path(manifest), Path(output)
    if manifest.resolve() == output.resolve():
        raise ValueError('Write the ImmuneLENS manifest to a separate file.')
    with manifest.open(newline='') as handle:
        samples = list(csv.DictReader(handle))
    if not samples or len({r['sample'] for r in samples}) != len(samples):
        raise ValueError('The source manifest must contain unique samples.')
    if len({r['bam'] for r in samples}) != len(samples):
        raise ValueError('The source manifest must contain unique BAMs.')
    prepared, audit = [], []
    for sample in samples:
        name = sample['sample']
        purity_path = Path(purity_dir) / (name + '.purple.purity.tsv')
        cn_path = Path(gene_cn_dir) / (name + '.purple.cnv.gene.tsv')
        with purity_path.open() as handle:
            purity_rows = list(csv.DictReader(handle, delimiter='\t'))
        if len(purity_rows) != 1:
            raise ValueError(name + ': expected one PURPLE purity row.')
        purity = purity_rows[0]
        if purity['status'] == 'NO_TUMOR':
            raise ValueError(name + ': NO_TUMOR samples must be excluded first.')
        value = float(purity['purity'])
        if not math.isfinite(value) or not 0 <= value <= 1:
            raise ValueError(name + ': invalid purity.')
        cn_bytes = cn_path.read_bytes()
        gene_rows = collections.defaultdict(list)
        lines = cn_bytes.decode().splitlines()
        fields = lines[0].split('\t')
        for line in lines[1:]:
            if not any('\t' + definition[0] + '\t' in line for definition in PROXIES.values()):
                continue
            gene = dict(zip(fields, line.split('\t')))
            if gene['isCanonical'] == 'true':
                gene_rows[gene['gene']].append(gene)
        row = {'sample': name, 'bam': sample['bam'], 'genome': 'hg38',
               'purity': purity['purity'], 'purple_status': purity['status']}
        for locus, (symbol, chromosome, start, end) in PROXIES.items():
            hits = gene_rows[symbol]
            if len(hits) != 1:
                raise ValueError(f'{name}: expected one canonical {symbol} CN record.')
            gene = hits[0]
            if (gene['chromosome'].removeprefix('chr'), int(gene['start']), int(gene['end'])) != (chromosome, start, end):
                raise ValueError(f'{name}: {symbol} coordinates do not match the hg38 annotation.')
            low, high = float(gene['minCopyNumber']), float(gene['maxCopyNumber'])
            if not all(math.isfinite(v) and v >= 0 for v in (low, high)) or low > high:
                status = 'invalid_proxy_cn'
            elif low != high:
                status = 'proxy_gene_heterogeneous'
            else:
                status = 'proxy_gene_uniform'
            cn = gene['minCopyNumber'] if status == 'proxy_gene_uniform' else 'NA'
            row[locus + '_cn'] = cn
            row[locus + '_cn_status'] = status
            audit.append({'sample': name, 'locus': locus, 'proxy_gene': symbol,
                          'genome': 'hg38', 'chromosome': chromosome, 'start': start, 'end': end,
                          'purity': purity['purity'], 'purple_status': purity['status'],
                          'local_cn': cn, 'cn_status': status,
                          'minCopyNumber': gene['minCopyNumber'], 'maxCopyNumber': gene['maxCopyNumber'],
                          'purity_file': str(purity_path.resolve()),
                          'purity_sha256': hashlib.sha256(purity_path.read_bytes()).hexdigest(),
                          'gene_cn_file': str(cn_path.resolve()),
                          'gene_cn_sha256': hashlib.sha256(cn_bytes).hexdigest()})
        prepared.append(row)
    fields = ['sample', 'bam', 'genome', 'purity', 'purple_status']
    fields += [locus + '_cn' for locus in PROXIES]
    fields += [locus + '_cn_status' for locus in PROXIES]
    for path, rows, columns, delimiter in [(output, prepared, fields, ','),
            (Path(str(output) + '.audit.tsv'), audit, list(audit[0]), '\t')]:
        with path.open('w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=columns, delimiter=delimiter)
            writer.writeheader()
            writer.writerows(rows)
    summary = {'manifest_samples': len(prepared), 'sample_locus_cn_records': len(audit),
               'samples_with_all_four_cn': sum(all(r[locus + '_cn'] != 'NA' for locus in PROXIES) for r in prepared),
               'cn_status_counts': dict(collections.Counter(r['cn_status'] for r in audit)),
               'source_manifest': str(manifest.resolve()),
               'source_manifest_sha256': hashlib.sha256(manifest.read_bytes()).hexdigest(),
               'manifest_sha256': hashlib.sha256(output.read_bytes()).hexdigest(),
               'method': 'ImmuneLENS paper nearby-gene proxies; use PURPLE canonical gene CN only when minCopyNumber equals maxCopyNumber.',
               'paper': 'https://doi.org/10.1038/s41588-025-02086-5',
               'note': 'NA CN skips adjustment for that locus, retaining its raw estimate. IGH matched-normal correction is not supplied.'}
    Path(str(output) + '.qc.json').write_text(json.dumps(summary, indent=2) + '\n')
    return summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for argument in ('input', 'purity-dir', 'gene-cn-dir', 'output'):
        parser.add_argument('--' + argument, required=True)
    args = parser.parse_args()
    print(json.dumps(prepare(args.input, args.purity_dir, args.gene_cn_dir, args.output), indent=2))
