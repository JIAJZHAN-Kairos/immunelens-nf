#!/usr/bin/env python3
"""Read only ImmuneLENS loci from an indexed local/S3 BAM."""
import argparse
import csv
import gzip
import json
import os
from pathlib import Path
import subprocess
import sys
from urllib.parse import urlparse


def run(command, **kwargs):
    result = subprocess.run(command, stderr=subprocess.PIPE, **kwargs)
    if result.returncode:
        # samtools can print signed URLs; keep credentials out of task logs.
        message = result.stderr.decode(errors='replace')
        if 'https://' in message:
            message = 'Remote BAM read failed; check S3 permissions, index and network access.'
        raise RuntimeError(f'{command[0]} exited {result.returncode}: {message.strip()}')
    return result


def bam_source(bam):
    if bam.startswith('s3://'):
        parsed = urlparse(bam)
        key = parsed.path.lstrip('/')
        index = None
        for candidate in (key + '.bai', key[:-4] + '.bai', key + '.csi', key[:-4] + '.csi'):
            check = subprocess.run(['aws', 's3api', 'head-object', '--bucket', parsed.netloc,
                                    '--key', candidate], capture_output=True)
            if check.returncode == 0:
                index = f's3://{parsed.netloc}/{candidate}'
                break
        if index is None:
            raise RuntimeError(f'No accessible BAI/CSI index for {bam}; an index is required.')
        local_index = Path('input' + Path(index).suffix)
        run(['aws', 's3', 'cp', index, str(local_index), '--only-show-errors'], stdout=subprocess.PIPE)
        url = run(['aws', 's3', 'presign', bam, '--expires-in', '43200'], stdout=subprocess.PIPE)
        return url.stdout.decode().strip(), local_index, index
    path = Path(bam).resolve()
    for index in (Path(str(path) + '.bai'), path.with_suffix('.bai'),
                  Path(str(path) + '.csi'), path.with_suffix('.csi')):
        if index.is_file():
            return str(path), index, str(index)
    raise RuntimeError(f'No BAI/CSI index for {path}; an index is required.')


def detect_genome(header):
    contigs = {}
    for line in header.splitlines():
        if line.startswith('@SQ\t'):
            fields = dict(field.split(':', 1) for field in line.split('\t')[1:])
            contigs[fields['SN']] = int(fields['LN'])
    signatures = {'hg38': {'7': 159345973, '14': 107043718},
                  'hg19': {'7': 159138663, '14': 107349540}}
    for genome, lengths in signatures.items():
        for prefix in ('chr', ''):
            if all(contigs.get(prefix + chromosome) == length for chromosome, length in lengths.items()):
                return genome, prefix
    raise RuntimeError('BAM chr7/chr14 lengths do not match hg19 or hg38; refusing to guess the genome.')


def main():
    os.environ.setdefault('AWS_PAGER', '')
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sample', required=True)
    parser.add_argument('--bam', required=True)
    parser.add_argument('--original-bam', help='Original manifest path retained in provenance after local staging.')
    parser.add_argument('--regions', required=True)
    args = parser.parse_args()
    source, index, index_uri = bam_source(args.bam)
    header = run(['samtools', 'view', '-H', source], stdout=subprocess.PIPE).stdout.decode()
    genome, prefix = detect_genome(header)
    output = Path(args.sample + '.coverage')
    output.mkdir()
    (output / 'bam_header.sam').write_text(header)
    counts = {}
    with open(args.regions) as handle:
        regions = [row for row in csv.DictReader(handle, delimiter='\t') if row['genome'] == genome]
    expected = {'TCRA', 'TCRB', 'TCRG', 'IGH'} if genome == 'hg38' else {'TCRA', 'TCRB', 'TCRG'}
    if {row['locus'] for row in regions} != expected:
        raise RuntimeError('The ImmuneLENS region table is incomplete.')
    for row in regions:
        region = f"{prefix}{row['chromosome']}:{row['start']}-{row['end']}"
        command = ['samtools', 'depth', '-q', '20', '-Q', '20', '-r', region,
                   '-X', source, str(index)]
        raw = output / (row['locus'] + '.txt')
        with raw.open('wb') as handle:
            run(command, stdout=handle)
        count = 0
        with raw.open('rb') as handle, gzip.open(str(raw) + '.gz', 'wb') as compressed:
            for line in handle:
                compressed.write(line)
                count += 1
        raw.unlink()
        counts[row['locus']] = count
        print(f'{args.sample}: {genome} {row["locus"]}, {count} covered positions', flush=True)
    (output / 'sample.json').write_text(json.dumps({
        'sample': args.sample, 'bam': args.original_bam or args.bam, 'bam_index': index_uri,
        'genome': genome, 'coverage_positions': counts,
        'base_quality': 20, 'mapping_quality': 20,
        'zero_coverage_positions': 'omitted, matching upstream getCovFromBam_WGS',
        'correction': 'unadjusted',
    }, indent=2) + '\n')


if __name__ == '__main__':
    try:
        main()
    except (RuntimeError, OSError, ValueError) as error:
        print(f'ERROR: {error}', file=sys.stderr)
        sys.exit(1)
