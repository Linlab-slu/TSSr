# Test clustering functionality

test_that("clusterTSS produces tagClusters", {
    data(exampleTSSr)

    expect_type(exampleTSSr@tagClusters, "list")
    expect_true(length(exampleTSSr@tagClusters) > 0)

    first_cluster <- exampleTSSr@tagClusters[[1]]
    expect_true(is.data.frame(first_cluster) ||
                    inherits(first_cluster, "data.table"))
    expect_true(nrow(first_cluster) > 0)
})

test_that("consensusClusters has valid structure", {
    data(exampleTSSr)

    expect_type(exampleTSSr@consensusClusters, "list")
    expect_true(length(exampleTSSr@consensusClusters) > 0)

    cc <- exampleTSSr@consensusClusters[[1]]
    expect_true(is.data.frame(cc))
    expect_true(nrow(cc) > 0)
})

test_that("tagClusters contain required columns", {
    data(exampleTSSr)

    tc <- exampleTSSr@tagClusters[[1]]
    col_names <- names(tc)

    has_position_info <- any(grepl("chr|start|end|strand|pos", col_names, ignore.case = TRUE))
    expect_true(has_position_info)
    expect_true(ncol(tc) >= 3)
})
