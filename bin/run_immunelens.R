#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(ImmuneLENS))
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) %in% c(2L, 3L, 4L))
script <- sub('^--file=', '', grep('^--file=', commandArgs(), value = TRUE)[[1]])
source(file.path(dirname(script), 'correct_igh.R'))
sample <- args[[1]]
coverage_dir <- args[[2]]
metadata <- jsonlite::fromJSON(file.path(coverage_dir, 'sample.json'))
corrections <- if (length(args) >= 3L) jsonlite::fromJSON(args[[3]]) else list()
normal_dir <- if (length(args) == 4L) args[[4]] else NULL
normal_metadata <- if (!is.null(normal_dir)) jsonlite::fromJSON(file.path(normal_dir, 'sample.json')) else NULL
if (isTRUE(normal_metadata$provided) && normal_metadata$genome != metadata$genome) {
    stop('Matched normal and tumour BAM genomes differ.', call. = FALSE)
}
if (length(corrections) && corrections$genome != metadata$genome) {
    stop('Manifest CN genome differs from the BAM genome.', call. = FALSE)
}
output_dir <- paste0(sample, '.immunelens')
dir.create(output_dir)
summaries <- list()
estimates <- list()

for (locus in c('TCRA', 'TCRB', 'TCRG', 'IGH')) {
    purity <- if (is.null(corrections$purity)) NA_real_ else corrections$purity
    local_cn <- corrections[[paste0(locus, '_cn')]]
    if (is.null(local_cn)) local_cn <- NA_real_
    cn_status <- corrections[[paste0(locus, '_cn_status')]]
    if (is.null(cn_status)) cn_status <- if (is.finite(local_cn)) 'user_supplied' else 'not_supplied'
    can_adjust <- is.finite(purity) && is.finite(local_cn)
    correction <- if (can_adjust) 'purity_local_cn' else 'unadjusted'
    status <- 'ok'
    note <- ''
    fit <- NULL
    fraction <- raw_fraction <- adjusted_fraction <- NA_real_
    high_fraction <- NA
    igh_correction <- if (locus != 'IGH') 'not_applicable' else
        if (isTRUE(normal_metadata$provided)) 'not_run' else 'not_provided'
    igh_uncorrected_fraction <- igh_germline_fraction <- igh_combined_fraction <- NA_real_
    igh_somatic_qc <- igh_somatic_selected <- NA
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
                if (locus == 'IGH') {
                    corrected <- correct_igh(coverage, normal_dir, normal_metadata, fit, sample, output_dir)
                    fit <- corrected$fit
                    audit <- corrected$audit
                    igh_correction <- audit$state
                    igh_uncorrected_fraction <- audit$uncorrected_fraction
                    igh_germline_fraction <- audit$germline_fraction
                    igh_combined_fraction <- audit$combined_fraction
                    igh_somatic_qc <- audit$somatic_qc
                    igh_somatic_selected <- audit$somatic_selected
                    note <- audit$note
                }
                raw_fraction <- fit[[1]][[column]][[1]]
                if (can_adjust) {
                    saveRDS(fit, file.path(output_dir, paste0(locus, '.raw.rds')))
                    fit <- lapply(fit, function(output) adjustImmuneLENS(output,
                        purity = purity, local.cn = local_cn, vdj.gene = locus))
                    adjusted_fraction <- fit[[1]][[paste0(column, '.adj')]][[1]]
                    high_fraction <- fit[[1]]$highTcellFlag[[1]]
                    if (!is.finite(adjusted_fraction)) stop('Non-finite adjusted fraction.', call. = FALSE)
                }
                fraction <- if (can_adjust) adjusted_fraction else raw_fraction
                if (fraction < 0 || fraction > 1) {
                    status <- 'out_of_range'
                    note <- paste(note, 'Estimated fraction is outside [0,1]; inspect local copy number and model fit.')
                } else if (isTRUE(high_fraction)) {
                    status <- 'high_cell_fraction'
                    note <- paste(note, 'Adjusted fraction exceeds the non-tumour fraction (1-purity).')
                }
                fit[[1]]$locus <- locus
                fit[[1]]$genome <- metadata$genome
                fit[[1]]$correction <- correction
                summaries[[locus]] <- fit[[1]]
                saveRDS(fit, file.path(output_dir, paste0(locus, '.rds')))
                write.table(fit[[2]], file.path(output_dir, paste0(locus, '.segments.tsv')),
                    sep = '\t', quote = FALSE, row.names = FALSE, na = 'NA')
                write.table(fit[[3]], file.path(output_dir, paste0(locus, '.model.tsv')),
                    sep = '\t', quote = FALSE, row.names = FALSE, na = 'NA')
            }
        }
    }
    if (length(corrections) && !can_adjust) {
        note <- paste(note, 'Purity/local-CN adjustment skipped:', cn_status)
    }
    if (locus == 'IGH' && igh_correction %in% c('not_run', 'not_provided')) {
        note <- paste(note, 'IGH matched-normal locus correction:', igh_correction)
    }
    estimates[[locus]] <- data.frame(sample, locus, genome = metadata$genome, status,
        cell_fraction = fraction, raw_cell_fraction = raw_fraction,
        adjusted_cell_fraction = adjusted_fraction, purity = purity, local_cn = local_cn,
        cn_status = cn_status, high_cell_fraction_flag = high_fraction,
        correction = correction, igh_correction = igh_correction,
        igh_uncorrected_fraction = igh_uncorrected_fraction,
        igh_germline_fraction = igh_germline_fraction, igh_combined_fraction = igh_combined_fraction,
        igh_somatic_qc = igh_somatic_qc, igh_somatic_selected = igh_somatic_selected,
        message = trimws(note))
    message(sample, ' ', locus, ': ', status)
}
write.table(do.call(rbind, estimates), file.path(output_dir, 'estimates.tsv'),
    sep = '\t', quote = FALSE, row.names = FALSE, na = 'NA')
if (length(summaries)) {
    write.table(dplyr::bind_rows(summaries), file.path(output_dir, 'upstream_summaries.tsv'),
        sep = '\t', quote = FALSE, row.names = FALSE, na = 'NA')
}
file.copy(file.path(coverage_dir, 'sample.json'), file.path(output_dir, 'sample.json'))
if (length(args) >= 3L) file.copy(args[[3]], file.path(output_dir, 'corrections.json'))
if (!is.null(normal_dir)) file.copy(file.path(normal_dir, 'sample.json'), file.path(output_dir, 'normal.json'))
writeLines(capture.output(sessionInfo()), file.path(output_dir, 'sessionInfo.txt'))
