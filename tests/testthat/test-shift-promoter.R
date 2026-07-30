make_shift_test_object <- function() {
    object <- TSSr(
        sampleLabels = c("control", "treat"),
        sampleLabelsMerged = c("control", "treat"),
        mergeIndex = c(1, 2)
    )
    object@librarySizes <- c(control = 1e6, treat = 1e6)
    object@assignedClusters <- list(
        control = data.table::data.table(
            cluster = c("cluster_1", "cluster_2"),
            strand = c("+", "+"),
            dominant_tss = c(100, 200),
            tags = c(600, 400),
            gene = c("gene_1", "gene_1")
        ),
        treat = data.table::data.table(
            cluster = c("cluster_1", "cluster_2"),
            strand = c("+", "+"),
            dominant_tss = c(100, 200),
            tags = c(400, 600),
            gene = c("gene_1", "gene_1")
        )
    )
    object
}

test_that("shiftPromoter returns numeric statistics", {
    object <- make_shift_test_object()

    result <- shiftPromoter(
        object,
        comparePairs = list(c("control", "treat")),
        pval = 1
    )
    shifts <- result@PromoterShift[["control_VS_treat"]]

    expect_type(shifts$gene, "character")
    expect_type(shifts$Ds, "double")
    expect_type(shifts$pval, "double")
    expect_type(shifts$padj, "double")
})
