test_that("callEnhancer preserves annotated cluster data across repeated calls", {
    data(exampleTSSr)
    object <- annotateCluster(
        exampleTSSr,
        clusters = "consensusClusters",
        filterCluster = TRUE
    )
    names_before <- lapply(
        object@unassignedClusters,
        function(cluster_table) data.table::copy(names(cluster_table))
    )
    expect_true(all(vapply(
        names_before,
        function(column_names) {
            all(c("gene", "inCoding") %in% column_names)
        },
        logical(1)
    )))

    first_warnings <- character()
    first <- withCallingHandlers(
        callEnhancer(object),
        warning = function(warning) {
            first_warnings <<- c(first_warnings, conditionMessage(warning))
            invokeRestart("muffleWarning")
        }
    )
    first_enhancers <- lapply(first@enhancers, data.table::copy)
    expect_length(first_warnings, 0L)
    expect_identical(lapply(object@unassignedClusters, names), names_before)
    expect_identical(lapply(first@unassignedClusters, names), names_before)

    second_warnings <- character()
    second <- withCallingHandlers(
        callEnhancer(first),
        warning = function(warning) {
            second_warnings <<- c(second_warnings, conditionMessage(warning))
            invokeRestart("muffleWarning")
        }
    )
    expect_length(second_warnings, 0L)
    expect_identical(lapply(first@unassignedClusters, names), names_before)
    expect_equal(second@enhancers, first_enhancers)
})
