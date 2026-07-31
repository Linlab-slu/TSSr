test_that("mergeSamples returns changes without modifying caller symbols", {
    data(exampleTSSr)
    object <- exampleTSSr
    before <- tssr_content(object)

    result <- mergeSamples(object)

    expect_tssr_content_equal(object, before)
    expect_s4_class(result, "TSSr")
    expect_equal(
        result@TSSprocessedMatrix$control,
        rowSums(result@TSSrawMatrix[, c("SL01", "SL02")])
    )
    expect_equal(
        result@TSSprocessedMatrix$treat,
        rowSums(result@TSSrawMatrix[, c("SL03", "SL04")])
    )
    expect_equal(
        unname(result@librarySizes),
        unname(colSums(result@TSSprocessedMatrix[, c("control", "treat")]))
    )
    expect_identical(result@normalizationStatus, "raw")

    toto <- function(x) mergeSamples(x)
    wrapped_input <- exampleTSSr
    wrapped_before <- tssr_content(wrapped_input)
    wrapped_result <- toto(wrapped_input)

    expect_tssr_content_equal(wrapped_input, wrapped_before)
    expect_s4_class(wrapped_result, "TSSr")
    wrapped_input <- toto(wrapped_input)
    expect_equal(
        wrapped_input@TSSprocessedMatrix,
        wrapped_result@TSSprocessedMatrix
    )
})

test_that("normalizeTSS returns normalized data without modifying its input", {
    data(exampleTSSr)
    object <- mergeSamples(exampleTSSr)
    before <- tssr_content(object)

    result <- normalizeTSS(object)

    expect_tssr_content_equal(object, before)
    expect_s4_class(result, "TSSr")
    normalized_sums <- colSums(
        result@TSSprocessedMatrix[
            , .SD, .SDcols = result@sampleLabelsMerged
        ]
    )
    expect_equal(unname(normalized_sums), rep(1e6, length(normalized_sums)))
    expect_false(isTRUE(all.equal(
        result@TSSprocessedMatrix,
        object@TSSprocessedMatrix
    )))
})

test_that("filterTSS returns filtered data without modifying its input", {
    data(exampleTSSr)
    object <- normalizeTSS(mergeSamples(exampleTSSr))
    before <- tssr_content(object)

    result <- filterTSS(object, method = "TPM", tpmLow = 10)

    expect_tssr_content_equal(object, before)
    expect_s4_class(result, "TSSr")
    expect_gt(nrow(result@TSSprocessedMatrix), 0)
    expect_lt(
        nrow(result@TSSprocessedMatrix),
        nrow(object@TSSprocessedMatrix)
    )
})

test_that("filterTSS returns its unchanged input for an unknown method", {
    data(exampleTSSr)
    object <- normalizeTSS(mergeSamples(exampleTSSr))
    before <- tssr_content(object)

    expect_message(
        result <- filterTSS(object, method = "unknown"),
        "No filtering method is defined",
        fixed = TRUE
    )

    expect_tssr_content_equal(object, before)
    expect_s4_class(result, "TSSr")
    expect_tssr_content_equal(result, before)
})

test_that("getTSS imports data without modifying its input", {
    fixture <- system.file(
        "extdata", "example-tss-table.tsv",
        package = "TSSr",
        mustWork = TRUE
    )
    sample_labels <- paste0("SL0", seq_len(4L))
    object <- TSSr(
        genomeName = "BSgenome.Scerevisiae.UCSC.sacCer3",
        inputFiles = fixture,
        inputFilesType = "TSStable",
        sampleLabels = sample_labels,
        sampleLabelsMerged = c("control", "treat"),
        mergeIndex = c(1, 1, 2, 2)
    )
    before <- tssr_content(object)

    result <- getTSS(object)

    expect_tssr_content_equal(object, before)
    expect_s4_class(result, "TSSr")
    imported <- TSSmatrix(result, data = "raw")
    expected <- data.table::fread(fixture)
    expected[, pos := as.integer(pos)]
    data.table::setorder(expected, strand, chr, pos)

    expect_equal(nrow(imported), 100L)
    expect_identical(
        names(imported),
        c("chr", "pos", "strand", sample_labels)
    )
    expect_setequal(unique(imported$strand), c("+", "-"))
    expect_equal(imported, as.data.frame(expected))
    expect_equal(
        unname(librarySizes(result)),
        unname(colSums(expected[, ..sample_labels]))
    )
    expect_identical(result@normalizationStatus, "raw")
})
