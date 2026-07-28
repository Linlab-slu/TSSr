# Test downstream analysis: shapeCluster, shiftPromoter, callEnhancer
# Prepare data once at the top to avoid repeating expensive workflow steps

# Top-level workflow setup mirrors .skip_on_bioc_spb(): on Bioconductor
# build machines we skip the heavy prep so the SPB 15-min budget holds;
# everywhere else (devtools::test, CI, user installs) it runs normally.
.on_bioc_spb <- identical(Sys.info()[["user"]], "biocbuild") ||
    grepl("bbs-", R.home(), fixed = TRUE) ||
    nzchar(Sys.getenv("BBS_HOME")) ||
    identical(Sys.getenv("IS_BIOC_BUILD_MACHINE"), "true")

data(exampleTSSr)
if (!.on_bioc_spb) {
    exampleTSSr <- mergeSamples(exampleTSSr)
    exampleTSSr <- normalizeTSS(exampleTSSr)
    exampleTSSr <- filterTSS(exampleTSSr, method = "TPM", tpmLow = 0.1)
    exampleTSSr <- clusterTSS(exampleTSSr,
        method = "peakclu", clusterThreshold = 1,
        useMultiCore = FALSE
    )
    exampleTSSr <- consensusCluster(exampleTSSr, useMultiCore = FALSE)
}

test_that("shapeCluster calculates shape scores with PSS method", {
    .skip_on_bioc_spb()
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
    .skip_on_bioc_spb()
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

test_that("shiftPromoter detects promoter shifts", {
    .skip_on_bioc_spb()
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
    .skip_on_bioc_spb()
    ## callEnhancer requires annotated clusters; skip if not available
    skip_if(length(exampleTSSr@assignedClusters) == 0,
        "No assignedClusters in exampleTSSr")
    skip_if(length(exampleTSSr@unassignedClusters) == 0,
        "No unassignedClusters in exampleTSSr")

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
