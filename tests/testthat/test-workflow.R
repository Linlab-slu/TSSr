# Test the core TSSr workflow: mergeSamples -> normalizeTSS -> filterTSS ->
# clusterTSS -> consensusCluster
# Run workflow once, test each step's output

data(exampleTSSr)

test_that("mergeSamples merges raw TSS data correctly", {
    mergeSamples(exampleTSSr)

    processed <- exampleTSSr@TSSprocessedMatrix
    expect_true(nrow(processed) > 0)
    merged_labels <- exampleTSSr@sampleLabelsMerged
    expect_true(all(merged_labels %in% names(processed)))
})

test_that("normalizeTSS normalizes to TPM", {
    normalizeTSS(exampleTSSr)

    processed <- exampleTSSr@TSSprocessedMatrix
    merged_labels <- exampleTSSr@sampleLabelsMerged
    first_col <- processed[[merged_labels[1]]]
    nonzero <- first_col[first_col > 0]
    ## TPM values should have decimals
    expect_true(any(nonzero != round(nonzero)))
})

test_that("filterTSS with TPM method reduces rows", {
    rows_before <- nrow(exampleTSSr@TSSprocessedMatrix)
    filterTSS(exampleTSSr, method = "TPM", tpmLow = 0.1)
    rows_after <- nrow(exampleTSSr@TSSprocessedMatrix)

    expect_true(rows_after > 0)
    expect_true(rows_after <= rows_before)
})

test_that("clusterTSS produces tagClusters", {
    clusterTSS(exampleTSSr,
        method = "peakclu", clusterThreshold = 1,
        useMultiCore = FALSE
    )

    tc <- exampleTSSr@tagClusters
    expect_type(tc, "list")
    expect_true(length(tc) > 0)

    first_tc <- tc[[1]]
    expect_true(is.data.frame(first_tc) || inherits(first_tc, "data.table"))
    expect_true(all(c("chr", "start", "end", "strand") %in% names(first_tc)))
})

test_that("consensusCluster produces consensus clusters", {
    consensusCluster(exampleTSSr, useMultiCore = FALSE)

    cc <- exampleTSSr@consensusClusters
    expect_type(cc, "list")
    expect_true(length(cc) > 0)

    first_cc <- cc[[1]]
    expect_true(is.data.frame(first_cc) || inherits(first_cc, "data.table"))
    expect_true(all(c("chr", "start", "end", "strand") %in% names(first_cc)))
})
