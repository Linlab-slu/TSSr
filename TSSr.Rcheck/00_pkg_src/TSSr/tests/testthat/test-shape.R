# Test promoter shape functionality

test_that("clusterShape slot has valid structure", {
  data(exampleTSSr)

  if (length(exampleTSSr@clusterShape) > 0) {
    expect_type(exampleTSSr@clusterShape, "list")

    # Check first shape result
    first_shape <- exampleTSSr@clusterShape[[1]]

    if (is.data.frame(first_shape) && nrow(first_shape) > 0) {
      # Shape data should have shape-related columns
      col_names <- names(first_shape)
      # Common shape metrics: PSS, SI, IQR, etc.
      has_shape_metrics <- any(grepl("shape|PSS|SI|IQR|width", col_names, ignore.case = TRUE)) ||
        ncol(first_shape) >= 4
      expect_true(has_shape_metrics)
    }
  }
})

test_that("shape values are numeric", {
  data(exampleTSSr)

  if (length(exampleTSSr@clusterShape) > 0) {
    first_shape <- exampleTSSr@clusterShape[[1]]

    if (is.data.frame(first_shape) && nrow(first_shape) > 0) {
      # At least some columns should be numeric (shape scores)
      numeric_cols <- sapply(first_shape, is.numeric)
      expect_true(any(numeric_cols))
    }
  }
})
