#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(ImmuneLENS))
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
sample <- args[[1]]
coverage_dir <- args[[2]]
metadata <- jsonlite::fromJSON(file.path(coverage_dir, 'sample.json'))
output_dir <- paste0(sample, '.immunelens')
dir.create(output_dir)
summaries <- list()
estimates <- list()

for (locus in c('TCRA', 'TCRB', 'TCRG', 'IGH')) {
    status <- 'ok'
    note <- ''
    fit <- NULL
    fraction <- NA_real_
    if (metadata$genome == 'hg19' && locus == 'IGH') {
        status <- 'unsupported_genome'
        note <- 'Upstream IGH analysis supports hg38 only.'
    } else if (metadata$coverage_positions[[locus]] == 0L) {
        status <- 'no_coverage'
        note <- 'No positions meet base-quality 20 and mapping-quality 20.'
    } else {
        coverage <- suppressMessages(loadCov(file.path(coverage_dir, paste0(locus, '.txt.gz'))))
        fit <- tryCatch(
            runImmuneLENS(coverage, vdj.gene = locus, hg19_or_38 = metadata$genome,
                GC_correct = TRUE, removed_flag = TRUE, sample_name = sample),
            error = function(error) {
                if (grepl('All positions have been removed due to low coverage', conditionMessage(error), fixed = TRUE)) {
                    status <<- 'insufficient_coverage'
                    note <<- conditionMessage(error)
                    return(NULL)
                }
                stop(paste(sample, locus, conditionMessage(error)), call. = FALSE)
            })
        if (status == 'ok') {
            column <- paste0(locus, if (locus == 'IGH') '.bcell.fraction' else '.tcell.fraction')
            if (!is.list(fit) || length(fit) != 3L || !is.data.frame(fit[[1]]) ||
                nrow(fit[[1]]) != 1L || !column %in% names(fit[[1]]) ||
                !is.finite(fit[[1]][[column]][[1]])) {
                status <- 'no_estimate'
                note <- 'The upstream model did not return a finite cell-fraction estimate.'
            } else {
                fraction <- fit[[1]][[column]][[1]]
                if (fraction < 0 || fraction > 1) {
                    status <- 'out_of_range'
                    note <- 'Raw upstream fraction is outside [0,1]; inspect local copy number and model fit.'
                }
                fit[[1]]$locus <- locus
                fit[[1]]$genome <- metadata$genome
                fit[[1]]$correction <- 'unadjusted'
                summaries[[locus]] <- fit[[1]]
                saveRDS(fit, file.path(output_dir, paste0(locus, '.rds')))
                write.table(fit[[2]], file.path(output_dir, paste0(locus, '.segments.tsv')),
                    sep = '\t', quote = FALSE, row.names = FALSE, na = 'NA')
                write.table(fit[[3]], file.path(output_dir, paste0(locus, '.model.tsv')),
                    sep = '\t', quote = FALSE, row.names = FALSE, na = 'NA')
            }
        }
    }
    if (locus == 'IGH') {
        note <- paste(note, 'Tumour-only IGH: no matched-normal haplotype or somatic-CNA correction applied.')
    }
    estimates[[locus]] <- data.frame(sample, locus, genome = metadata$genome, status,
        cell_fraction = fraction, correction = 'unadjusted', message = trimws(note))
    message(sample, ' ', locus, ': ', status)
}
write.table(do.call(rbind, estimates), file.path(output_dir, 'estimates.tsv'),
    sep = '\t', quote = FALSE, row.names = FALSE, na = 'NA')
if (length(summaries)) {
    write.table(dplyr::bind_rows(summaries), file.path(output_dir, 'upstream_summaries.tsv'),
        sep = '\t', quote = FALSE, row.names = FALSE, na = 'NA')
}
file.copy(file.path(coverage_dir, 'sample.json'), file.path(output_dir, 'sample.json'))
writeLines(capture.output(sessionInfo()), file.path(output_dir, 'sessionInfo.txt'))
