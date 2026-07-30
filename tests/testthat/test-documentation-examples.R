test_that("reviewer-requested examples call their documented functions", {
    documented_functions <- c(
        "callEnhancer", "clusterTSS", "exportClustersTable",
        "exportClustersToBed", "exportDETable", "exportEnhancerTable",
        "exportShapeTable", "exportShiftTable", "exportTSStable",
        "exportTSStoBedgraph", "getTSS", "plotCorrelation", "plotDE",
        "plotInterQuantile", "plotShape", "plotTssPCA", "plotTSS",
        "shapeCluster", "shiftPromoter"
    )
    documentation <- tools::Rd_db("TSSr")
    examples_tag <- paste0(intToUtf8(92), "examples")

    for (function_name in documented_functions) {
        rd <- documentation[[paste0(function_name, ".Rd")]]
        examples <- rd[vapply(
            rd,
            function(node) identical(attr(node, "Rd_tag"), examples_tag),
            logical(1)
        )]
        example_text <- paste(unlist(examples), collapse = "")

        expect_match(
            example_text,
            paste0("\\b", function_name, "\\s*\\("),
            info = paste(function_name, "must be called by its man example")
        )
    }
})
