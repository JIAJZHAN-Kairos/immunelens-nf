source(commandArgs(trailingOnly = TRUE)[[1]])
class_switch_warning <- 'Not enough bases with coverage in IGH class switch region'
runImmuneLENS <- function(...) {
    warning(class_switch_warning)
    stop('Model syntax is empty.')
}
result <- fit_with_coverage_qc(NULL, 'IGH', 'hg38', 'TEST')
stopifnot(is.null(result$fit), result$status == 'insufficient_class_switch_coverage',
    result$note == class_switch_warning, is.null(fit_igh(NULL, 'TEST')))
# The same parser error without this specific IGH warning must remain fatal.
runImmuneLENS <- function(...) stop('Model syntax is empty.')
error <- tryCatch(fit_with_coverage_qc(NULL, 'IGH', 'hg38', 'TEST'), error = identity)
stopifnot(inherits(error, 'error'), conditionMessage(error) == 'Model syntax is empty.')
# Do not reinterpret unrelated warnings or failures from TCR models.
runImmuneLENS <- function(...) {
    warning(class_switch_warning)
    stop('Model syntax is empty.')
}
error <- suppressWarnings(tryCatch(fit_with_coverage_qc(NULL, 'TCRA', 'hg38', 'TEST'), error = identity))
stopifnot(inherits(error, 'error'))
runImmuneLENS <- function(...) {
    warning('Unrelated model warning')
    stop('Unexpected model failure')
}
error <- suppressWarnings(tryCatch(fit_with_coverage_qc(NULL, 'IGH', 'hg38', 'TEST'), error = identity))
stopifnot(inherits(error, 'error'), conditionMessage(error) == 'Unexpected model failure')
runImmuneLENS <- function(...) stop('All positions have been removed due to low coverage')
stopifnot(fit_with_coverage_qc(NULL, 'TCRB', 'hg38', 'TEST')$status == 'insufficient_coverage')
