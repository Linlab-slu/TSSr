# Test plot functions — all write PDF to tempdir

expect_valid_pdf <- function(path) {
    expect_true(file.exists(path))
    expect_gt(unname(file.info(path)$size), 1000)
    header <- readBin(path, what = "raw", n = 5L)
    expect_identical(rawToChar(header), "%PDF-")
}

test_that("plotCorrelation creates PDF", {
    data(exampleTSSr)
    tmpdir <- tempdir()
    withr::with_dir(tmpdir, {
        invalid_graphics_warning <- FALSE
        withCallingHandlers(
            plotCorrelation(exampleTSSr, samples = "all"),
            warning = function(warning) {
                if (grepl(
                    "argument 1 does not name a graphical parameter",
                    conditionMessage(warning),
                    fixed = TRUE
                )) {
                    invalid_graphics_warning <<- TRUE
                    invokeRestart("muffleWarning")
                }
            }
        )
        expect_false(invalid_graphics_warning)
        expect_valid_pdf("TSS_correlation_plot_of_all_samples.pdf")
        unlink("TSS_correlation_plot_of_all_samples.pdf")
    })
})

test_that("plotTssPCA creates PDF", {
    data(exampleTSSr)
    tmpdir <- tempdir()
    withr::with_dir(tmpdir, {
        plotTssPCA(exampleTSSr, TSS.threshold = 10)
        expect_valid_pdf("PCA_plot.pdf")
        unlink("PCA_plot.pdf")
    })
})

test_that("plotInterQuantile creates PDF", {
    data(exampleTSSr)
    tmpdir <- tempdir()
    expect_gt(length(exampleTSSr@clusterShape), 0L)
    withr::with_dir(tmpdir, {
        plotInterQuantile(exampleTSSr, samples = "all", tagsThreshold = 1)
        expect_valid_pdf("Interquantile_plot_of_ALL_samples.pdf")
        unlink("Interquantile_plot_of_ALL_samples.pdf")
    })
})

test_that("plotShape creates PDF", {
    data(exampleTSSr)
    tmpdir <- tempdir()
    expect_gt(length(exampleTSSr@clusterShape), 0L)
    withr::with_dir(tmpdir, {
        plotShape(exampleTSSr, samples = "all")
        expect_valid_pdf("Shape_plot_of_ALL_samples.pdf")
        unlink("Shape_plot_of_ALL_samples.pdf")
    })
})

test_that("plotDE creates PDF", {
    data(exampleTSSr)
    tmpdir <- tempdir()
    expect_gt(length(exampleTSSr@DEtables), 0L)
    withr::with_dir(tmpdir, {
        plotDE(exampleTSSr, withGeneName = TRUE)
        expect_valid_pdf("Volcano_plot.pdf")
        unlink("Volcano_plot.pdf")
    })
})
