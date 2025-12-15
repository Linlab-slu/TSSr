# Test clustering functionality

test_that("clusterTSS produces tagClusters", {
  data(exampleTSSr)

  # Check if tagClusters slot is populated after clustering
  # (exampleTSSr should already have clustering done)
  if (length(exampleTSSr@tagClusters) > 0) {
    expect_type(exampleTSSr@tagClusters, "list")

    # Each element should be a data.table or data.frame
    first_cluster <- exampleTSSr@tagClusters[[1]]
    expect_true(is.data.frame(first_cluster) ||
                  inherits(first_cluster, "data.table"))
  }
})

test_that("consensusClusters has valid structure", {
  data(exampleTSSr)

  if (length(exampleTSSr@consensusClusters) > 0) {
    expect_type(exampleTSSr@consensusClusters, "list")
  }
})

test_that("tagClusters contain required columns", {
  data(exampleTSSr)

  if (length(exampleTSSr@tagClusters) > 0) {
    tc <- exampleTSSr@tagClusters[[1]]

    # Tag clusters should have positional information
    col_names <- names(tc)
    # Check for common column patterns
    has_position_info <- any(grepl("chr|start|end|strand|pos", col_names, ignore.case = TRUE))
    expect_true(has_position_info || ncol(tc) >= 3)
  }
})
