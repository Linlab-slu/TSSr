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
    object@normalizationStatus <- "raw"
    object
}

test_that("TPM filtering accepts normalized libraries below one million", {
    object <- make_filter_test_object(list(small = c(1, 2, 3)))
    object <- normalizeTSS(object)

    expect_identical(object@normalizationStatus, "normalized")
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

test_that("TPM filtering preserves normalized state after rows are removed", {
    counts <- c(rep(1, 9000), 973703)
    object <- make_filter_test_object(list(sample = counts))
    object <- normalizeTSS(object)
    object <- filterTSS(object, method = "TPM", tpmLow = 2)

    expect_lt(
        sum(object@TSSprocessedMatrix$sample),
        1e6
    )
    expect_no_error(
        result <- filterTSS(object, method = "TPM", tpmLow = 2)
    )
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

test_that("TPM filtering rejects objects recorded as raw", {
    object <- make_filter_test_object(list(raw = c(1, 2, 3)))

    expect_error(
        filterTSS(object, method = "TPM", tpmLow = 0.1),
        "Data must be normalized.",
        fixed = TRUE
    )
})

test_that("normalizeTSS rejects repeated normalization", {
    object <- make_filter_test_object(list(sample = c(1, 2, 3)))
    object <- normalizeTSS(object)

    expect_error(
        normalizeTSS(object),
        "data is already normalized",
        fixed = TRUE
    )
})

test_that("poisson filtering accepts raw count data", {
    object <- make_filter_test_object(list(raw = c(1, 2, 3)))

    expect_no_error(
        result <- filterTSS(object, method = "poisson", normalization = FALSE)
    )
    expect_s4_class(result, "TSSr")
    expect_identical(result@normalizationStatus, "raw")
})

test_that("poisson filtering records normalization of its result", {
    object <- make_filter_test_object(list(raw = c(10, 20, 30)))

    expect_no_error(
        result <- filterTSS(object, method = "poisson", normalization = TRUE)
    )
    expect_identical(result@normalizationStatus, "normalized")
})

test_that("unknown state falls back to integer-valued legacy data", {
    object <- make_filter_test_object(list(sample = c(1, 2, 3)))
    object@normalizationStatus <- NA_character_

    expect_no_error(result <- normalizeTSS(object))
    expect_identical(result@normalizationStatus, "normalized")
})

test_that("normalization status validity accepts only known scalar states", {
    object <- TSSr()

    expect_true(validObject(object))

    object@normalizationStatus <- character()
    expect_true(validObject(object))

    object@normalizationStatus <- c("raw", "normalized")
    expect_error(validObject(object), "normalizationStatus")

    object@normalizationStatus <- "other"
    expect_error(validObject(object), "normalizationStatus")
})
