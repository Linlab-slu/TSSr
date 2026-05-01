# Test promoter shape functionality

test_that("clusterShape slot has valid structure", {
    data(exampleTSSr)

    expect_type(exampleTSSr@clusterShape, "list")
    expect_true(length(exampleTSSr@clusterShape) > 0)

    first_shape <- exampleTSSr@clusterShape[[1]]
    expect_true(is.data.frame(first_shape))
    expect_true(nrow(first_shape) > 0)

    col_names <- names(first_shape)
    has_shape_metrics <- any(grepl("shape|PSS|SI|IQR|width", col_names, ignore.case = TRUE))
    expect_true(has_shape_metrics)
})

test_that("shape values are numeric", {
    data(exampleTSSr)

    first_shape <- exampleTSSr@clusterShape[[1]]
    numeric_cols <- vapply(first_shape, is.numeric, logical(1))
    expect_true(any(numeric_cols))
})
