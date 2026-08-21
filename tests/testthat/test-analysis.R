# Test downstream analysis: shapeCluster, shiftPromoter, callEnhancer
# Use the bundled upstream results and recompute each method under test.

data(exampleTSSr)

test_that("shapeCluster calculates shape scores with PSS method", {
    object <- exampleTSSr
    object@clusterShape <- list()
    before <- tssr_content(object)
    result <- shapeCluster(object, clusters = "consensusClusters", method = "PSS",
        useMultiCore = FALSE)

    expect_tssr_content_equal(object, before)
    expect_s4_class(result, "TSSr")
    cs <- result@clusterShape
    expect_type(cs, "list")
    expect_true(length(cs) > 0)

    first_shape <- cs[[1]]
    expect_true(is.data.frame(first_shape))
    expect_true("shape.score" %in% names(first_shape))
    expect_true(is.numeric(first_shape$shape.score))
})

test_that("deGene calculates differential expression tables", {
    object <- exampleTSSr
    object@DEtables <- list()
    object@TAGtables <- list()
    before <- tssr_content(object)
    result <- deGene(
        object,
        comparePairs = list(c("control", "treat")),
        pval = 0.01,
        useMultiCore = FALSE
    )

    expect_tssr_content_equal(object, before)
    expect_s4_class(result, "TSSr")
    expect_true(length(result@DEtables) > 0)
    expect_true(length(result@TAGtables) > 0)
})

test_that("deGene supports mean dispersion fitting for small data sets", {
    data(exampleTSSr)
    object <- exampleTSSr
    object@DEtables <- list()
    object@TAGtables <- list()

    expect_no_warning(
        result <- deGene(
            object,
            comparePairs = list(c("control", "treat")),
            pval = 0.01,
            useMultiCore = FALSE,
            fitType = "mean"
        )
    )

    expect_s4_class(result, "TSSr")
    expect_gt(length(result@DEtables), 0L)
})

test_that("shiftPromoter detects promoter shifts", {
    object <- exampleTSSr
    object@PromoterShift <- list()
    before <- tssr_content(object)
    result <- shiftPromoter(object,
        comparePairs = list(c("control", "treat")),
        pval = 0.01
    )

    expect_tssr_content_equal(object, before)
    expect_s4_class(result, "TSSr")
    ps <- result@PromoterShift
    expect_type(ps, "list")
    expect_true(length(ps) > 0)
    expect_true("control_VS_treat" %in% names(ps))
})

test_that("callEnhancer identifies enhancer candidates when data available", {
    expect_gt(length(exampleTSSr@assignedClusters), 0)
    expect_gt(length(exampleTSSr@unassignedClusters), 0)

    object <- exampleTSSr
    object@enhancers <- list()
    before <- tssr_content(object)
    result <- callEnhancer(object, flanking = 400, dis2gene = 2000)

    expect_tssr_content_equal(object, before)
    expect_s4_class(result, "TSSr")
    en <- result@enhancers
    expect_type(en, "list")
    expect_true(length(en) > 0)
})
