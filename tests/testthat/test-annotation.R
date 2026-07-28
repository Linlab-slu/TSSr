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

test_that("annotateCluster creates valid filtered cluster tables", {
    data(exampleTSSr)
    result <- annotateCluster(
        exampleTSSr,
        clusters = "consensusClusters",
        filterCluster = TRUE
    )

    for (sample_label in result@sampleLabelsMerged) {
        filtered <- result@filteredClusters[[sample_label]]
        expected_names <- names(
            result@assignedClusters[[sample_label]]
        )[seq_len(12)]

        expect_s3_class(filtered, "data.frame")
        expect_identical(names(filtered), expected_names)
        expect_false(any(grepl("^V[0-9]+$", names(filtered))))
        expect_gt(nrow(filtered), 2L)
    }
})
