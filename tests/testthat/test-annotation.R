# Test annotation functionality

test_that("assignedClusters have gene annotations", {
    data(exampleTSSr)

    expect_type(exampleTSSr@assignedClusters, "list")
    expect_true(length(exampleTSSr@assignedClusters) > 0)

    first_assigned <- exampleTSSr@assignedClusters[[1]]
    expect_true(is.data.frame(first_assigned))
    expect_true(nrow(first_assigned) > 0)

    col_names <- names(first_assigned)
    has_gene_col <- any(grepl("gene", col_names, ignore.case = TRUE))
    expect_true(has_gene_col)
})

test_that("unassignedClusters slot is populated", {
    data(exampleTSSr)

    expect_type(exampleTSSr@unassignedClusters, "list")
    expect_true(length(exampleTSSr@unassignedClusters) > 0)
})
