get_tssr_generic <- function(name) {
    methods::getGeneric(name, where = asNamespace("TSSr"))
}

tssr_generic_names <- c(
    "annotateCluster", "callEnhancer", "clusterTSS", "consensusCluster",
    "deGene", "exportClustersTable", "exportClustersToBed", "exportDETable",
    "exportEnhancerTable", "exportShapeTable", "exportShiftTable",
    "exportTSStable", "exportTSStoBedgraph", "filterTSS", "getTSS",
    "mergeSamples", "normalizeTSS", "plotCorrelation", "plotDE",
    "plotInterQuantile", "plotShape", "plotTSS", "plotTssPCA",
    "shapeCluster", "shiftPromoter"
)

test_that("all exported TSSr generics are standard S4 generics", {
    for (name in tssr_generic_names) {
        generic <- get_tssr_generic(name)
        expect_true(
            methods::is(generic, "standardGeneric"),
            info = name
        )
    }
})

test_that("all exported TSSr generics dispatch only on object", {
    for (name in tssr_generic_names) {
        generic <- get_tssr_generic(name)
        expect_identical(generic@signature, "object", info = name)
    }
})

test_that("exported TSSr generics retain compatible TSSr methods", {
    exports <- getNamespaceExports("TSSr")

    for (name in tssr_generic_names) {
        generic <- get_tssr_generic(name)
        method <- methods::selectMethod(name, "TSSr")

        expect_true(name %in% exports, info = name)
        expect_identical(
            names(formals(method)),
            names(formals(generic)),
            info = name
        )
        expect_identical(as.character(method@target), "TSSr", info = name)
        expect_identical(as.character(method@defined), "TSSr", info = name)
    }
})
