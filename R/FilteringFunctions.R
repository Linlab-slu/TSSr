################################################################################################
.normalizationStatus <- function(object) {
    status <- tryCatch(
        object@normalizationStatus,
        error = function(e) character()
    )
    if (length(status) == 1 && !is.na(status)) {
        return(status)
    }

    sample_labels <- object@sampleLabelsMerged
    if (length(sample_labels) == 0 ||
        nrow(object@TSSprocessedMatrix) == 0) {
        return(NA_character_)
    }

    processed <- as.data.frame(object@TSSprocessedMatrix)
    values <- unlist(processed[sample_labels], use.names = FALSE)
    values <- values[is.finite(values)]
    if (length(values) == 0) {
        return(NA_character_)
    }

    tolerance <- .Machine$double.eps^0.5
    if (all(abs(values - round(values)) < tolerance)) "raw" else "normalized"
}

.isNormalized <- function(object) {
    identical(.normalizationStatus(object), "normalized")
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
