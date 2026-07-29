# Test the core TSSr workflow: mergeSamples -> normalizeTSS -> filterTSS ->
# clusterTSS -> consensusCluster
# Run workflow once, test each step's output

data(exampleTSSr)

test_that("mergeSamples merges raw TSS data correctly", {
    result <- mergeSamples(exampleTSSr)

    processed <- result@TSSprocessedMatrix
    expect_true(nrow(processed) > 0)
    merged_labels <- result@sampleLabelsMerged
    expect_true(all(merged_labels %in% names(processed)))
})

test_that("normalizeTSS normalizes to TPM", {
    object <- mergeSamples(exampleTSSr)
    result <- normalizeTSS(object)

    processed <- result@TSSprocessedMatrix
    merged_labels <- result@sampleLabelsMerged
    first_col <- processed[[merged_labels[1]]]
    nonzero <- first_col[first_col > 0]
    ## TPM values should have decimals
    expect_true(any(nonzero != round(nonzero)))
})

test_that("filterTSS with TPM method reduces rows", {
    object <- normalizeTSS(mergeSamples(exampleTSSr))
    rows_before <- nrow(object@TSSprocessedMatrix)
    result <- filterTSS(object, method = "TPM", tpmLow = 2)
    rows_after <- nrow(result@TSSprocessedMatrix)

    expect_true(rows_after > 0)
    expect_lt(rows_after, rows_before)
})

test_that("clusterTSS produces tagClusters", {
    object <- exampleTSSr
    object@tagClusters <- list()
    before <- tssr_content(object)
    result <- clusterTSS(object,
        method = "peakclu", clusterThreshold = 1,
        useMultiCore = FALSE
    )

    expect_tssr_content_equal(object, before)
    expect_s4_class(result, "TSSr")
    tc <- result@tagClusters
    expect_type(tc, "list")
    expect_true(length(tc) > 0)

    first_tc <- tc[[1]]
    expect_true(is.data.frame(first_tc) || inherits(first_tc, "data.table"))
    expect_true(nrow(first_tc) > 0)
    expect_true(all(c("chr", "start", "end", "strand") %in% names(first_tc)))
})

test_that("consensusCluster produces consensus clusters", {
    object <- exampleTSSr
    object@consensusClusters <- list()
    before <- tssr_content(object)
    result <- consensusCluster(object, useMultiCore = FALSE)

    expect_tssr_content_equal(object, before)
    expect_s4_class(result, "TSSr")
    cc <- result@consensusClusters
    expect_type(cc, "list")
    expect_true(length(cc) > 0)

    first_cc <- cc[[1]]
    expect_true(is.data.frame(first_cc) || inherits(first_cc, "data.table"))
    expect_true(all(c("chr", "start", "end", "strand") %in% names(first_cc)))
})
