# Exercise fallback with the real public uncorrected tumour fit.
suppressPackageStartupMessages(library(ImmuneLENS))
args <- commandArgs(trailingOnly = TRUE)
source(args[[1]])
demo <- args[[2]]
raw <- readRDS(file.path(demo, 'DEMO.immunelens/IGH.uncorrected.rds'))
coverage <- loadCov(file.path(demo, 'DEMO.coverage/IGH.txt.gz'))
out <- file.path(demo, 'qc-checks')
dir.create(out)
missing <- correct_igh(coverage, NULL, NULL, raw, 'DEMO', out)
stopifnot(missing$audit$state == 'not_provided', identical(missing$fit, raw))
empty <- correct_igh(coverage, NULL, list(provided = TRUE,
    coverage_positions = list(IGH = 0L)), raw, 'DEMO', out)
stopifnot(empty$audit$state == 'normal_no_coverage', identical(empty$fit, raw))
low_dir <- file.path(out, 'low.normal.coverage')
dir.create(low_dir)
handle <- gzfile(file.path(low_dir, 'IGH.txt.gz'), 'wt')
write.table(data.frame(chr = 'chr14', pos = 105526000 + 1:100, reads = 1),
    handle, sep = '\t', quote = FALSE,
    row.names = FALSE, col.names = FALSE)
close(handle)
low <- correct_igh(coverage, low_dir, list(provided = TRUE,
    coverage_positions = list(IGH = 100L)), raw, 'DEMO', out)
stopifnot(low$audit$state == 'normal_insufficient_depth', identical(low$fit, raw))
cat('Missing, empty and low-depth normals retain the uncorrected model with explicit QC.\n')
