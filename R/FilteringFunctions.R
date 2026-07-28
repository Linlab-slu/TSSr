################################################################################################
.isNormalized <- function(object) {
    sample_labels <- object@sampleLabelsMerged
    library_sizes <- object@librarySizes
    if (!is.null(names(library_sizes)) &&
        all(sample_labels %in% names(library_sizes))) {
        library_sizes <- library_sizes[sample_labels]
    } else {
        library_sizes <- library_sizes[seq_along(sample_labels)]
    }
    column_sums <- vapply(
        sample_labels,
        function(sample_label) {
            sum(object@TSSprocessedMatrix[[sample_label]], na.rm = TRUE)
        },
        numeric(1)
    )

    distance_to_raw <- abs(column_sums - library_sizes)
    distance_to_tpm <- abs(column_sums - 1e6)
    distance_to_tpm <= distance_to_raw
}

################################################################################################
.filterWithPoisson <- function(data, coverageDepth, genomeSize, pVal) {
    # calculate lambda value (average)
    lambda <- coverageDepth / (genomeSize * 2)
    # get cutoff value
    cutoff <- qpois(pVal, lambda, lower.tail = FALSE, log.p = FALSE)
    # fiter tss table
    data[data < cutoff, ] <- 0
    return(data)
}
