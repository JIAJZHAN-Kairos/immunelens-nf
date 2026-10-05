# Coverage QC and matched-normal IGH correction using the pinned upstream implementation.
fit_with_coverage_qc <- function(coverage, locus, genome, sample) {
    tryCatch({
        fit <- withCallingHandlers(
            runImmuneLENS(coverage, vdj.gene = locus, hg19_or_38 = genome,
                GC_correct = TRUE, removed_flag = TRUE, sample_name = sample),
            warning = function(warning) {
                if (locus == 'IGH' && identical(conditionMessage(warning),
                    'Not enough bases with coverage in IGH class switch region')) {
                    # Upstream returns NULL constraints here, then lavaan fails parsing them.
                    stop(structure(list(message = conditionMessage(warning), call = NULL),
                        class = c('igh_coverage_error', 'error', 'condition')))
                }
            })
        list(fit = fit, status = 'ok', note = '')
    }, igh_coverage_error = function(error) {
        list(fit = NULL, status = 'insufficient_class_switch_coverage', note = conditionMessage(error))
    }, error = function(error) {
        if (grepl('All positions have been removed due to low coverage',
            conditionMessage(error), fixed = TRUE)) {
            return(list(fit = NULL, status = 'insufficient_coverage', note = conditionMessage(error)))
        }
        stop(error)
    })
}

igh_fraction <- function(fit) {
    if (!is.list(fit) || length(fit) != 3L || !is.data.frame(fit[[1]]) ||
        nrow(fit[[1]]) != 1L || !'IGH.bcell.fraction' %in% names(fit[[1]])) return(NA_real_)
    fit[[1]]$IGH.bcell.fraction[[1]]
}

fit_igh <- function(coverage, sample) {
    result <- fit_with_coverage_qc(coverage, 'IGH', 'hg38', sample)
    if (result$status != 'ok') message(sample, ' corrected IGH: ', result$status)
    result$fit
}

select_igh_correction <- function(germline_fraction, combined_fraction, somatic_qc) {
    # Supplementary Methods: accept combined correction only when it reduces B fraction.
    isTRUE(somatic_qc) && is.finite(combined_fraction) && combined_fraction < germline_fraction
}

correct_igh <- function(coverage, normal_dir, normal_metadata, raw_fit, sample, output_dir) {
    audit <- list(state = 'not_provided', uncorrected_fraction = igh_fraction(raw_fit),
        germline_fraction = NA_real_, combined_fraction = NA_real_,
        somatic_qc = FALSE, somatic_selected = FALSE,
        normal_median_covered_depth = NA_real_, tumor_median_covered_depth = median(coverage$reads),
        normal_low_bcell_assumed = TRUE,
        note = 'Matched normal not provided; IGH locus CNV correction was skipped.')
    result <- function() {
        jsonlite::write_json(audit, file.path(output_dir, 'IGH.correction.json'),
            pretty = TRUE, auto_unbox = TRUE, na = 'null')
        list(fit = raw_fit, audit = audit)
    }
    if (is.null(normal_metadata) || !isTRUE(normal_metadata$provided)) return(result())
    audit$state <- 'normal_no_coverage'
    audit$note <- 'Matched-normal IGH has no qualifying coverage; locus CNV correction was skipped.'
    if (normal_metadata$coverage_positions$IGH == 0L) return(result())
    normal <- suppressMessages(loadCov(file.path(normal_dir, 'IGH.txt.gz')))
    audit$normal_median_covered_depth <- median(normal$reads)
    if (!is.finite(audit$normal_median_covered_depth) || audit$normal_median_covered_depth <= 10) {
        audit$state <- 'normal_insufficient_depth'
        audit$note <- 'Normal median covered-position depth is <=10; IGH CNV correction was skipped.'
        return(result())
    }
    germline <- IGH_haplotype_norm_fun(normal)
    write.table(germline[[2]], file.path(output_dir, 'IGH.germline_regions.tsv'),
        sep = '\t', quote = FALSE, row.names = FALSE, na = 'NA')
    # Tumour normalisation uses observed coverage ratios, rather than diploid scaling.
    paired <- IGH_haplotype_norm_fun_tumour(normal, coverage, germline[[2]])
    audit$somatic_qc <- isTRUE(paired$somatic_correction_applied)
    write.table(paired$somatic_regions_df, file.path(output_dir, 'IGH.somatic_regions.tsv'),
        sep = '\t', quote = FALSE, row.names = FALSE, na = 'NA')
    germline_fit <- fit_igh(paired$cov.df_update2a, sample)
    audit$germline_fraction <- igh_fraction(germline_fit)
    if (!is.finite(audit$germline_fraction)) {
        audit$state <- 'germline_no_estimate'
        audit$note <- 'Germline-corrected model has no finite estimate; uncorrected IGH retained.'
        return(result())
    }
    saveRDS(raw_fit, file.path(output_dir, 'IGH.uncorrected.rds'))
    saveRDS(germline_fit, file.path(output_dir, 'IGH.germline.rds'))
    raw_fit <- germline_fit
    audit$state <- 'germline_only'
    audit$note <- 'Germline correction selected. Normal caller assumes low B-cell content (<10%).'
    if (audit$somatic_qc && audit$tumor_median_covered_depth > 20) {
        combined_fit <- fit_igh(paired$cov.df_update2b, sample)
        audit$combined_fraction <- igh_fraction(combined_fit)
        if (is.finite(audit$combined_fraction)) {
            saveRDS(combined_fit, file.path(output_dir, 'IGH.germline_somatic.rds'))
        }
        audit$somatic_selected <- select_igh_correction(audit$germline_fraction,
            audit$combined_fraction, audit$somatic_qc)
        if (audit$somatic_selected) {
            raw_fit <- combined_fit
            audit$state <- 'germline_somatic'
            audit$note <- 'Germline+somatic correction passed QC and reduced the B fraction. Normal caller assumes low B-cell content (<10%).'
        } else {
            audit$note <- paste(audit$note, 'Combined correction did not yield a finite lower B fraction.')
        }
    } else {
        audit$note <- paste(audit$note, 'Somatic correction skipped: CNA QC failed or tumour median covered-position depth <=20.')
    }
    result()
}
