# Test using exampleTSSr data

test_that("exampleTSSr data loads correctly", {
    data(exampleTSSr)
    expect_s4_class(exampleTSSr, "TSSr")
})

test_that("exampleTSSr has expected structure", {
    data(exampleTSSr)

    # Essential slots must be populated
    expect_true(length(exampleTSSr@sampleLabels) > 0)
    expect_true(nrow(exampleTSSr@TSSrawMatrix) > 0)
    expect_true(nrow(exampleTSSr@TSSprocessedMatrix) > 0)
})

test_that("exampleTSSr sampleLabels are valid", {
    data(exampleTSSr)

    expect_type(exampleTSSr@sampleLabels, "character")
    expect_equal(length(exampleTSSr@sampleLabels), 4)
    expect_equal(exampleTSSr@sampleLabels, c("SL01", "SL02", "SL03", "SL04"))
})

test_that("exampleTSSr TSSrawMatrix has correct columns", {
    data(exampleTSSr)

    raw_matrix <- exampleTSSr@TSSrawMatrix
    expect_true(nrow(raw_matrix) > 100)
    expect_true(ncol(raw_matrix) >= 3)
})

test_that("exampleTSSr has all expected slots populated", {
    data(exampleTSSr)

    expect_equal(exampleTSSr@inputFilesType, "bam")
    expect_true(length(exampleTSSr@tagClusters) > 0)
    expect_true(length(exampleTSSr@consensusClusters) > 0)
    expect_true(length(exampleTSSr@clusterShape) > 0)
    expect_true(length(exampleTSSr@assignedClusters) > 0)
    expect_true(length(exampleTSSr@unassignedClusters) > 0)
    expect_true(length(exampleTSSr@DEtables) > 0)
    expect_true(length(exampleTSSr@PromoterShift) > 0)
})
