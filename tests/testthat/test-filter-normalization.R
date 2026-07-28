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
    skip_if_not_installed("BSgenome.Scerevisiae.UCSC.sacCer3")
    object <- make_filter_test_object(list(small = c(1, 2, 3)))
    normalizeTSS(object)

    expect_no_error(filterTSS(object, method = "TPM", tpmLow = 0.1))
})

test_that("TPM filtering accepts normalized libraries above one million", {
    skip_if_not_installed("BSgenome.Scerevisiae.UCSC.sacCer3")
    object <- make_filter_test_object(
        list(large = c(1, 1, 1999998))
    )
    normalizeTSS(object)

    expect_no_error(filterTSS(object, method = "TPM", tpmLow = 0.1))
})

test_that("poisson filtering rejects normalized libraries below one million", {
    skip_if_not_installed("BSgenome.Scerevisiae.UCSC.sacCer3")
    object <- make_filter_test_object(list(small = c(1, 2, 3)))
    normalizeTSS(object)

    expect_error(
        filterTSS(object, method = "poisson"),
        "Raw count data required for poisson method.",
        fixed = TRUE
    )
})

test_that("normalization detection checks every merged sample", {
    skip_if_not_installed("BSgenome.Scerevisiae.UCSC.sacCer3")
    object <- make_filter_test_object(list(
        small = c(1, 2, 3),
        large = c(1, 1, 1999998)
    ))
    normalizeTSS(object)

    expect_no_error(filterTSS(object, method = "TPM", tpmLow = 0.1))
})
