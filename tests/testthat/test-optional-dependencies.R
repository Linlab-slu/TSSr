test_that("feature-specific heavy packages are not forced namespace imports", {
    optional_packages <- c("DESeq2", "Gviz", "ggfortify", "calibrate")
    namespace_imports <- names(getNamespaceImports("TSSr"))

    expect_length(intersect(optional_packages, namespace_imports), 0L)
})

test_that("feature-specific heavy packages are declared in Suggests", {
    optional_packages <- c("DESeq2", "Gviz", "ggfortify", "calibrate")
    suggests_field <- utils::packageDescription("TSSr")$Suggests
    suggested_packages <- trimws(sub(
        "[[:space:]]*\\(.*$", "", strsplit(suggests_field, ",")[[1]]
    ))

    expect_true(all(optional_packages %in% suggested_packages))
})

test_that("missing feature packages produce an actionable error", {
    expect_error(
        TSSr:::`.requireSuggestedPackage`(
            "TSSrPackageThatDoesNotExist", "exampleFeature()"
        ),
        "exampleFeature.*BiocManager::install"
    )
})
