make_unmerged_import <- function() {
    input <- system.file(
        "extdata", "example-tss-table.tsv",
        package = "TSSr", mustWork = TRUE
    )
    TSSr(
        genomeName = "BSgenome.Scerevisiae.UCSC.sacCer3",
        inputFiles = input,
        inputFilesType = "TSStable",
        sampleLabels = c("SL01", "SL02", "SL03", "SL04")
    )
}

test_that("mergeSamples is a one-to-one operation when grouping is omitted", {
    object <- getTSS(make_unmerged_import())
    before <- tssr_content(object)

    result <- mergeSamples(object)

    expect_tssr_content_equal(object, before)
    expect_identical(
        TSSmatrix(result, data = "processed"),
        TSSmatrix(object, data = "raw")
    )
    expect_identical(librarySizes(result), librarySizes(object))
})

test_that("workflow runs with separate samples when mergeSamples is skipped", {
    object <- getTSS(make_unmerged_import())
    object <- normalizeTSS(object)
    result <- clusterTSS(
        object,
        method = "peakclu",
        peakDistance = 100,
        clusterThreshold = 1,
        useMultiCore = FALSE
    )

    expect_identical(names(tagClusters(result)), object@sampleLabels)
    expect_length(tagClusters(result), length(object@sampleLabels))
    expect_true(all(vapply(tagClusters(result), nrow, integer(1)) > 0))
})
