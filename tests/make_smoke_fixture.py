"""Make two tiny, indexed hg38 BAMs for a real end-to-end QC smoke run."""
import csv
import shutil
from pathlib import Path
import subprocess
import sys

regions, destination = sys.argv[1:]
output = Path(destination).resolve()
output.mkdir(parents=True, exist_ok=True)
header = '@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr7\tLN:159345973\n@SQ\tSN:chr14\tLN:107043718\n'
with open(regions) as handle:
    positions = sorted({(int(row['chromosome']), int(row['start']))
                        for row in csv.DictReader(handle, delimiter='\t') if row['genome'] == 'hg38'})
for sample in ('EMPTY', 'LOWCOV'):
    sam = output / (sample + '.sam')
    text = header
    if sample == 'LOWCOV':
        for i, (chromosome, position) in enumerate(positions):
            text += f'read{i}\t0\tchr{chromosome}\t{position}\t60\t100M\t*\t0\t0\t' + 'A' * 100 + '\t' + 'I' * 100 + '\n'
    sam.write_text(text)
    bam = output / (sample + '.bam')
    subprocess.run(['samtools', 'view', '-b', '-o', str(bam), str(sam)], check=True)
    subprocess.run(['samtools', 'index', str(bam)], check=True)
(output / 'normal').mkdir(exist_ok=True)
normal = output / 'normal/LOWCOV.bam'
# Same basename as tumour, to check staging isolation; empty normal is a QC case.
shutil.copyfile(output / 'EMPTY.bam', normal)
shutil.copyfile(output / 'EMPTY.bam.bai', str(normal) + '.bai')
(output / 'samplesheet.csv').write_text('sample,bam,normal_bam,genome,purity,TCRA_cn,TCRB_cn,TCRG_cn,IGH_cn\n' +
    ''.join(f'{s},{output / (s + ".bam")},{normal if s == "LOWCOV" else ""},hg38,0.5,3,3,{"NA" if s == "LOWCOV" else 3},3\n'
            for s in ('EMPTY', 'LOWCOV')))
