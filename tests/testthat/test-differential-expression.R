make_three_group_de_object <- function() {
    data(exampleTSSr)
    object <- exampleTSSr
    object@TSSrawMatrix <- data.table::copy(exampleTSSr@TSSrawMatrix)
    object@TSSrawMatrix[["SL05"]] <- object@TSSrawMatrix[["SL03"]]
    object@TSSrawMatrix[["SL06"]] <- object@TSSrawMatrix[["SL04"]]
    object@sampleLabels <- c(object@sampleLabels, "SL05", "SL06")
    object@sampleLabelsMerged <- c("control", "treat", "third")
    object@mergeIndex <- c(1, 1, 2, 2, 3, 3)
    object@assignedClusters <- lapply(
        exampleTSSr@assignedClusters,
        data.table::copy
    )
    object@assignedClusters[["third"]] <- data.table::copy(
        object@assignedClusters[["treat"]]
    )
    object@TAGtables <- list()
    object@DEtables <- list()
    object
}

# Test differential expression functionality

test_that("DEtables slot has valid structure", {
    data(exampleTSSr)

    expect_type(exampleTSSr@DEtables, "list")
    expect_true(length(exampleTSSr@DEtables) > 0)

    first_de <- exampleTSSr@DEtables[[1]]
    expect_type(first_de, "list")
    expect_true(length(first_de) > 0)
    expect_true("DEtable" %in% names(first_de))
})

test_that("DEtable contains expected DESeq2 columns", {
    data(exampleTSSr)

    first_de <- exampleTSSr@DEtables[[1]]
    de_table <- first_de$DEtable

    expect_true(is.data.frame(de_table))
    expect_true(nrow(de_table) > 0)

    col_names <- names(de_table)
    deseq_cols <- c("baseMean", "log2FoldChange", "pvalue", "padj")
    has_deseq_cols <- all(deseq_cols %in% col_names)
    expect_true(has_deseq_cols)
})

test_that("deGene recomputes TAG tables after assigned clusters change", {
    data(exampleTSSr)
    object <- exampleTSSr
    object@assignedClusters <- lapply(
        object@assignedClusters,
        function(clusters) clusters[-1L]
    )
    object@DEtables <- list()

    result <- deGene(
        object,
        comparePairs = list(c("control", "treat")),
        pval = 0.01,
        useMultiCore = FALSE
    )

    expect_identical(names(result@TAGtables), c("control", "treat"))
    expect_equal(
        vapply(result@TAGtables, nrow, integer(1)),
        vapply(object@assignedClusters, nrow, integer(1))
    )
})

test_that("deGene retains TAG tables across multiple comparisons", {
    object <- make_three_group_de_object()

    result <- deGene(
        object,
        comparePairs = list(
            c("control", "treat"),
            c("control", "third")
        ),
        pval = 0.01,
        useMultiCore = FALSE
    )

    expect_identical(
        names(result@TAGtables),
        c("control", "treat", "third")
    )
    expect_equal(
        vapply(result@TAGtables, nrow, integer(1)),
        vapply(object@assignedClusters, nrow, integer(1))
    )
    expected_tag_summary <- paste(
        paste(
            names(result@TAGtables),
            vapply(result@TAGtables, nrow, integer(1))
        ),
        collapse = "; "
    )
    expect_true(any(grepl(
        paste0("TAG tables: ", expected_tag_summary),
        capture.output(show(result)),
        fixed = TRUE
    )))
})

test_that("combined and separate deGene comparisons return the same results", {
    object <- make_three_group_de_object()
    common_arguments <- list(pval = 0.01, useMultiCore = FALSE)

    combined <- do.call(deGene, c(list(
        object = object,
        comparePairs = list(
            c("control", "treat"),
            c("control", "third")
        )
    ), common_arguments))
    control_treat <- do.call(deGene, c(list(
        object = object,
        comparePairs = list(c("control", "treat"))
    ), common_arguments))
    control_third <- do.call(deGene, c(list(
        object = object,
        comparePairs = list(c("control", "third"))
    ), common_arguments))

    expect_equal(
        combined@DEtables[["control_VS_treat"]],
        control_treat@DEtables[["control_VS_treat"]]
    )
    expect_equal(
        combined@DEtables[["control_VS_third"]],
        control_third@DEtables[["control_VS_third"]]
    )
    expect_equal(
        combined@TAGtables[["control"]],
        control_treat@TAGtables[["control"]]
    )
    expect_equal(
        combined@TAGtables[["third"]],
        control_third@TAGtables[["third"]]
    )
})
