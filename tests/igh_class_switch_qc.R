suppressPackageStartupMessages(library(ImmuneLENS))
args <- commandArgs(trailingOnly = TRUE)
source(args[[1]])
coverage <- loadCov(file.path(args[[2]], 'IGH.txt.gz'))
segments <- get('vdj_seg_list', asNamespace('ImmuneLENS'))$IGH_hg38
ighm <- segments[segments$segName == 'IGHM', ]
stopifnot(nrow(ighm) == 1L)
gap <- coverage[!(coverage$pos >= ighm$start & coverage$pos <= ighm$end), ]
stopifnot(nrow(gap) > 0L, nrow(gap) < nrow(coverage))
# Reproduce the real pinned upstream failure before checking its QC conversion.
warning_seen <- FALSE
error <- tryCatch(withCallingHandlers(
    runImmuneLENS(gap, vdj.gene = 'IGH', hg19_or_38 = 'hg38',
        GC_correct = TRUE, removed_flag = TRUE, sample_name = 'CS_GAP'),
    warning = function(warning) {
        if (conditionMessage(warning) == 'Not enough bases with coverage in IGH class switch region')
            warning_seen <<- TRUE
        invokeRestart('muffleWarning')
    }), error = identity)
stopifnot(warning_seen, inherits(error, 'error'),
    grepl('Model syntax is empty', conditionMessage(error), fixed = TRUE))
qc <- fit_with_coverage_qc(gap, 'IGH', 'hg38', 'CS_GAP')
stopifnot(is.null(qc$fit), qc$status == 'insufficient_class_switch_coverage',
    is.null(fit_igh(gap, 'CS_GAP')))
connection <- gzfile(file.path(args[[3]], 'IGH.txt.gz'), 'wt')
write.table(data.frame(chromosome = 'chr14', gap), connection,
    sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
close(connection)
