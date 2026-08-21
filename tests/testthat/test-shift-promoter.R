make_shift_test_object <- function() {
    object <- methods::new(
        "TSSr",
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

test_that("shiftPromoter tests reconstructed raw counts", {
    object <- make_shift_test_object()
    object@librarySizes <- c(control = 500000, treat = 800000)
    raw_counts <- matrix(
        c(200, 300, 480, 320),
        nrow = 2,
        dimnames = list(
            c("cluster_2", "cluster_1"),
            c("control", "treat")
        )
    )
    expected_pval <- chisq.test(raw_counts)$p.value

    result <- shiftPromoter(
        object,
        comparePairs = list(c("control", "treat")),
        pval = 1
    )
    shifts <- result@PromoterShift[["control_VS_treat"]]

    expect_equal(
        unname(shifts$pval),
        unname(expected_pval),
        tolerance = 1e-15
    )
})

test_that("shiftPromoter summarizes sparse chi-squared warnings once", {
    object <- make_shift_test_object()
    object@assignedClusters <- list(
        control = data.table::data.table(
            cluster = paste0("cluster_", 1:4),
            strand = rep("+", 4),
            dominant_tss = c(100, 200, 300, 400),
            tags = c(0, 1, 0, 2),
            gene = rep(c("gene_1", "gene_2"), each = 2)
        ),
        treat = data.table::data.table(
            cluster = paste0("cluster_", 1:4),
            strand = rep("+", 4),
            dominant_tss = c(100, 200, 300, 400),
            tags = c(1, 0, 2, 0),
            gene = rep(c("gene_1", "gene_2"), each = 2)
        )
    )
    warnings <- character()

    result <- withCallingHandlers(
        shiftPromoter(
            object,
            comparePairs = list(c("control", "treat")),
            pval = 1
        ),
        warning = function(warning_condition) {
            warnings <<- c(warnings, conditionMessage(warning_condition))
            invokeRestart("muffleWarning")
        }
    )

    expect_s4_class(result, "TSSr")
    expect_length(warnings, 1L)
    expect_match(warnings, "2 gene-level test.*p-values may be unreliable")
})
