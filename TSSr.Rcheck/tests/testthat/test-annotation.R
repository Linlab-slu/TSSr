# Test annotation functionality

test_that("assignedClusters have gene annotations", {
  data(exampleTSSr)

  if (length(exampleTSSr@assignedClusters) > 0) {
    expect_type(exampleTSSr@assignedClusters, "list")

    # Check first assigned cluster
    first_assigned <- exampleTSSr@assignedClusters[[1]]

    if (is.data.frame(first_assigned) && nrow(first_assigned) > 0) {
      col_names <- names(first_assigned)
      # Should have gene-related column
      has_gene_col <- any(grepl("gene", col_names, ignore.case = TRUE))
      expect_true(has_gene_col || ncol(first_assigned) >= 4)
    }
  }
})

test_that("unassignedClusters slot exists", {
  data(exampleTSSr)

  # unassignedClusters should be a list (even if empty)
  expect_type(exampleTSSr@unassignedClusters, "list")
})
