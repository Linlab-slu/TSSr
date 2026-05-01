# Test that sapply->vapply/lapply replacements produce correct results
# These tests exercise the code paths where sapply was replaced.

test_that("clusterTSS works correctly after sapply->lapply in ClusteringFunctions", {
    data(exampleTSSr)
    # clusterTSS calls .clusterTSS which uses lapply (was sapply)
    result <- clusterTSS(exampleTSSr, method = "peakclu", peakDistance = 100,
                         extensionDistance = 30, localThreshold = 0.02,
                         clusterThreshold = 1, useMultiCore = FALSE,
                         numCores = 1)
    expect_s4_class(result, "TSSr")
    expect_true(length(result@tagClusters) > 0)
    tc <- result@tagClusters[[1]]
    expect_true(is.data.frame(tc))
    expect_true(nrow(tc) > 0)
})

test_that("consensusCluster works after sapply->lapply in ConsensusMethods", {
    data(exampleTSSr)
    # consensusCluster uses lapply (was sapply) for chr/strand conversion
    result <- consensusCluster(exampleTSSr, useMultiCore = FALSE)
    expect_s4_class(result, "TSSr")
    expect_true(length(result@consensusClusters) > 0)
    cc <- result@consensusClusters[[1]]
    expect_true(is.data.frame(cc))
})

test_that("shapeCluster works after sapply->vapply in ShapeMethods", {
    data(exampleTSSr)
    # shapeCluster uses vapply (was sapply) for entropy calculation
    result <- consensusCluster(exampleTSSr, useMultiCore = FALSE)
    result <- shapeCluster(result, clusters = "consensusClusters", method = "PSS",
                           useMultiCore = FALSE)
    expect_s4_class(result, "TSSr")
    expect_true(length(result@clusterShape) > 0)
    shape <- result@clusterShape[[1]]
    expect_true(is.data.frame(shape))
    # shape.score column should be numeric
    expect_true(is.numeric(shape$shape.score))
})

test_that("deGene works after sapply->vapply in ExpressionFunctions", {
    data(exampleTSSr)
    # deGene uses vapply (was sapply) for tag counting and D.names
    result <- deGene(exampleTSSr,
                     comparePairs = list(c("control", "treat")),
                     pval = 0.01)
    expect_s4_class(result, "TSSr")
    expect_true(length(result@DEtables) > 0)
})

test_that("shiftPromoter works after sapply->vapply in ShiftingMethods", {
    data(exampleTSSr)
    # shiftPromoter uses vapply (was sapply) for D.names
    # suppress chisq.test warnings about approximation
    result <- suppressWarnings(
        shiftPromoter(exampleTSSr,
                      comparePairs = list(c("control", "treat")),
                      pval = 0.01)
    )
    expect_s4_class(result, "TSSr")
    expect_true(length(result@PromoterShift) > 0)
})
