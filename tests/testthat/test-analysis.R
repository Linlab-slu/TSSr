# Test downstream analysis: shapeCluster, shiftPromoter, callEnhancer
# Prepare data once at the top to avoid repeating expensive workflow steps

data(exampleTSSr)
mergeSamples(exampleTSSr)
normalizeTSS(exampleTSSr)
filterTSS(exampleTSSr, method = "TPM", tpmLow = 0.1)
clusterTSS(exampleTSSr,
    method = "peakclu", clusterThreshold = 1,
    useMultiCore = FALSE
)
consensusCluster(exampleTSSr, useMultiCore = FALSE)

test_that("shapeCluster calculates shape scores with PSS method", {
    shapeCluster(exampleTSSr, clusters = "consensusClusters", method = "PSS",
        useMultiCore = FALSE)

    cs <- exampleTSSr@clusterShape
    expect_type(cs, "list")
    expect_true(length(cs) > 0)

    first_shape <- cs[[1]]
    expect_true("shape.score" %in% names(first_shape))
    expect_true(is.numeric(first_shape$shape.score))
})

test_that("shiftPromoter detects promoter shifts", {
    shiftPromoter(exampleTSSr,
        comparePairs = list(c("control", "treat")),
        pval = 0.01
    )

    ps <- exampleTSSr@PromoterShift
    expect_type(ps, "list")
    expect_true(length(ps) > 0)
    expect_true("control_VS_treat" %in% names(ps))
})

test_that("callEnhancer identifies enhancer candidates when data available", {
    ## callEnhancer requires annotated clusters; skip if not available
    skip_if(length(exampleTSSr@assignedClusters) == 0,
        "No assignedClusters in exampleTSSr")
    skip_if(length(exampleTSSr@unassignedClusters) == 0,
        "No unassignedClusters in exampleTSSr")

    callEnhancer(exampleTSSr, flanking = 400, dis2gene = 2000)

    en <- exampleTSSr@enhancers
    expect_type(en, "list")
    expect_true(length(en) > 0)
})
