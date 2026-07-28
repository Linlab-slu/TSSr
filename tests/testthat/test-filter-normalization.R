make_filter_test_object <- function(counts) {
    sample_labels <- names(counts)
    raw <- data.table::data.table(
        chr = rep("chrI", length(counts[[1]])),
        pos = seq_along(counts[[1]]),
        strand = rep("+", length(counts[[1]]))
    )
    for (sample_label in sample_labels) {
        raw[[sample_label]] <- counts[[sample_label]]
    }

    object <- TSSr(
        genomeName = "BSgenome.Scerevisiae.UCSC.sacCer3",
        sampleLabels = sample_labels,
        sampleLabelsMerged = sample_labels,
        mergeIndex = seq_along(sample_labels)
    )
    object@TSSrawMatrix <- raw
    object@TSSprocessedMatrix <- data.table::copy(raw)
    object@librarySizes <- vapply(counts, sum, numeric(1))
    object
}

test_that("TPM filtering accepts normalized libraries below one million", {
    object <- make_filter_test_object(list(small = c(1, 2, 3)))
    object <- normalizeTSS(object)

    expect_no_error(result <- filterTSS(object, method = "TPM", tpmLow = 0.1))
    expect_s4_class(result, "TSSr")
})

test_that("TPM filtering accepts normalized libraries above one million", {
    object <- make_filter_test_object(
        list(large = c(1, 1, 1999998))
    )
    object <- normalizeTSS(object)

    expect_no_error(result <- filterTSS(object, method = "TPM", tpmLow = 0.1))
    expect_s4_class(result, "TSSr")
})

test_that("poisson filtering rejects normalized libraries below one million", {
    object <- make_filter_test_object(list(small = c(1, 2, 3)))
    object <- normalizeTSS(object)

    expect_error(
        filterTSS(object, method = "poisson"),
        "Raw count data required for poisson method.",
        fixed = TRUE
    )
})

test_that("TPM filtering rejects mixed normalized and raw samples", {
    counts <- list(normalized = c(1, 2, 3), raw = c(1, 1, 1999998))
    object <- make_filter_test_object(counts)
    object <- normalizeTSS(object)
    object@TSSprocessedMatrix[["raw"]] <- counts$raw

    expect_error(
        filterTSS(object, method = "TPM", tpmLow = 0.1),
        "Data must be normalized.",
        fixed = TRUE
    )
})

test_that("poisson filtering rejects mixed normalized and raw samples", {
    counts <- list(normalized = c(1, 2, 3), raw = c(1, 1, 1999998))
    object <- make_filter_test_object(counts)
    object <- normalizeTSS(object)
    object@TSSprocessedMatrix[["raw"]] <- counts$raw

    expect_error(
        filterTSS(object, method = "poisson"),
        "Raw count data required for poisson method.",
        fixed = TRUE
    )
})

test_that("poisson filtering accepts raw count data", {
    skip_if_not_installed("BSgenome.Scerevisiae.UCSC.sacCer3")
    object <- make_filter_test_object(list(raw = c(1, 2, 3)))

    expect_no_error(
        result <- filterTSS(object, method = "poisson", normalization = FALSE)
    )
    expect_s4_class(result, "TSSr")
})

test_that("normalization detection aligns named library sizes", {
    counts <- list(small = c(1, 2, 3), large = c(1, 1, 1999998))
    object <- make_filter_test_object(counts)
    object <- normalizeTSS(object)
    object@TSSprocessedMatrix[["small"]] <- counts$small
    object@librarySizes <- object@librarySizes[c("large", "small")]

    expect_error(
        filterTSS(object, method = "TPM", tpmLow = 0.1),
        "Data must be normalized.",
        fixed = TRUE
    )
})
