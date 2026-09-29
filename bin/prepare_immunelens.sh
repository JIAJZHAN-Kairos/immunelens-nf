#!/usr/bin/env bash
set -euo pipefail
export AWS_PAGER=''

revision=df34702b46c9a8e0e201469e1e80e65c9996716d
curl --fail --location --retry 3 "https://codeload.github.com/McGranahanLab/ImmuneLENS/tar.gz/${revision}" -o ImmuneLENS.tar.gz
curl --fail --location --retry 3 https://cran.r-project.org/src/contrib/Archive/restriktor/restriktor_0.6-50.tar.gz -o restriktor.tar.gz || \
    curl --fail --location --retry 3 https://cran.r-project.org/src/contrib/restriktor_0.6-50.tar.gz -o restriktor.tar.gz
python3 - <<'PY'
import hashlib
checksums = {
    'ImmuneLENS.tar.gz': 'daea052cbdd1a03e190c353a7fab332f136a2bc5a8f4a64516efd4ef7eecfe9f',
    'restriktor.tar.gz': 'd0d0f056e1db2d767065582733b296ea61db668cb3da88cc273ddbbb1343bdfb',
}
for filename, expected in checksums.items():
    with open(filename, 'rb') as handle:
        observed = hashlib.file_digest(handle, 'sha256').hexdigest()
    if observed != expected:
        raise SystemExit(f'{filename}: SHA-256 mismatch; refusing installation')
PY
mkdir immunelens-library
export R_LIBS_USER="$PWD/immunelens-library"
R CMD INSTALL --library="$R_LIBS_USER" restriktor.tar.gz
R CMD INSTALL --library="$R_LIBS_USER" ImmuneLENS.tar.gz
Rscript - <<'RS'
library(ImmuneLENS)
segments <- get('vdj_seg_list', envir = asNamespace('ImmuneLENS'))
chromosomes <- c(TCRA = '14', TCRB = '7', TCRG = '7', IGH = '14')
rows <- list()
for (genome in c('hg19', 'hg38')) {
    for (locus in names(chromosomes)) {
        if (genome == 'hg19' && locus == 'IGH') next
        definition <- segments[[paste0(locus, '_', genome)]]
        region <- definition[definition$segName == 'all', ]
        stopifnot(nrow(region) == 1L)
        rows[[length(rows) + 1L]] <- data.frame(genome, locus,
            chromosome = chromosomes[[locus]], start = region$start, end = region$end)
    }
}
write.table(do.call(rbind, rows), 'regions.tsv', sep = '\t', quote = FALSE, row.names = FALSE)
writeLines(c('ImmuneLENS revision: df34702b46c9a8e0e201469e1e80e65c9996716d',
    paste('ImmuneLENS version:', packageVersion('ImmuneLENS')),
    paste('restriktor version:', packageVersion('restriktor')),
    capture.output(sessionInfo())), 'software_versions.txt')
RS
samtools --version >> software_versions.txt
aws --version >> software_versions.txt 2>&1
tar -czf immunelens-library.tar.gz immunelens-library
