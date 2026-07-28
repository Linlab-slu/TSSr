# Test annotation functionality

test_that("annotateCluster preserves reference and consensus cluster inputs", {
    data(exampleTSSr)
    object <- exampleTSSr

    object@refTable <- data.table::copy(object@refTable)
    if ("subset" %in% names(object@refTable)) {
        object@refTable[, subset := NULL]
    }
    data.table::setkey(object@refTable, NULL)

    object@consensusClusters <- lapply(
        object@consensusClusters,
        function(cluster_table) {
            cluster_table <- data.table::copy(cluster_table)
            if ("subset" %in% names(cluster_table)) {
                cluster_table[, subset := NULL]
            }
            data.table::setkey(cluster_table, NULL)
            cluster_table
        }
    )
    object@assignedClusters <- list()
    object@unassignedClusters <- list()
    object@filteredClusters <- list()
    object_before <- tssr_content(object)

    ref_before <- list(
        names = data.table::copy(names(object@refTable)),
        rows = nrow(object@refTable),
        gene_id = data.table::copy(object@refTable$gene_id),
        key = data.table::key(object@refTable)
    )
    consensus_before <- lapply(
        object@consensusClusters,
        function(cluster_table) {
            list(
                names = data.table::copy(names(cluster_table)),
                rows = nrow(cluster_table),
                cluster = data.table::copy(cluster_table$cluster),
                key = data.table::key(cluster_table)
            )
        }
    )

    result <- annotateCluster(
        object,
        clusters = "consensusClusters",
        filterCluster = TRUE
    )

    expect_tssr_content_equal(object, object_before)
    expect_true(length(result@assignedClusters) > 0)
    expect_true(length(result@unassignedClusters) > 0)
    expect_true(length(result@filteredClusters) > 0)
    expect_identical(names(result@refTable), ref_before$names)
    expect_identical(nrow(result@refTable), ref_before$rows)
    expect_identical(result@refTable$gene_id, ref_before$gene_id)
    expect_identical(data.table::key(result@refTable), ref_before$key)
    expect_false("subset" %in% names(result@refTable))

    for (sample_label in result@sampleLabelsMerged) {
        cluster_table <- result@consensusClusters[[sample_label]]
        snapshot <- consensus_before[[sample_label]]
        expect_identical(names(cluster_table), snapshot$names)
        expect_identical(nrow(cluster_table), snapshot$rows)
        expect_identical(cluster_table$cluster, snapshot$cluster)
        expect_identical(data.table::key(cluster_table), snapshot$key)
        expect_false("subset" %in% names(cluster_table))
    }
})

test_that("annotateCluster returns the same annotations when called twice", {
    data(exampleTSSr)

    first <- annotateCluster(
        exampleTSSr,
        clusters = "consensusClusters",
        filterCluster = TRUE
    )
    second <- annotateCluster(
        first,
        clusters = "consensusClusters",
        filterCluster = TRUE
    )

    expect_equal(second@assignedClusters, first@assignedClusters)
    expect_equal(second@unassignedClusters, first@unassignedClusters)
})

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
