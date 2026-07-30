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

test_that("the GFF taxonomy data dependency is declared in Suggests", {
    suggests_field <- utils::packageDescription("TSSr")$Suggests
    suggested_packages <- trimws(sub(
        "[[:space:]]*\\(.*$", "", strsplit(suggests_field, ",")[[1]]
    ))

    expect_true("GenomeInfoDbData" %in% suggested_packages)
})

test_that("missing feature packages produce an actionable error", {
    expect_error(
        TSSr:::`.requireSuggestedPackage`(
            "TSSrPackageThatDoesNotExist", "exampleFeature()"
        ),
        "exampleFeature.*BiocManager::install"
    )
})

test_that("GFF annotation checks its taxonomy data dependency", {
    data(exampleTSSr)
    object <- exampleTSSr
    object@refTable <- data.table::data.table()
    object@refSource <- "unused-test-annotation.gff"

    testthat::local_mocked_bindings(
        .getGenome = function(...) NULL,
        .requireSuggestedPackage = function(package, feature) {
            expect_identical(package, "GenomeInfoDbData")
            expect_identical(feature, "GFF-based annotation")
            stop("taxonomy dependency guard reached", call. = FALSE)
        },
        .package = "TSSr"
    )

    expect_error(
        annotateCluster(object, filterCluster = FALSE),
        "taxonomy dependency guard reached"
    )
})

test_that("prefilled reference tables do not require taxonomy data", {
    data(exampleTSSr)

    testthat::local_mocked_bindings(
        .getGenome = function(...) NULL,
        .requireSuggestedPackage = function(...) {
            stop("unexpected taxonomy dependency check", call. = FALSE)
        },
        .package = "TSSr"
    )

    expect_no_error(
        annotateCluster(exampleTSSr, filterCluster = FALSE)
    )
})
