# Test plot functions — all write PDF to tempdir

test_that("plotCorrelation creates PDF", {
    data(exampleTSSr)
    tmpdir <- tempdir()
    withr::with_dir(tmpdir, {
        plotCorrelation(exampleTSSr, samples = "all")
        expect_true(file.exists("TSS_correlation_plot_of_all_samples.pdf"))
        unlink("TSS_correlation_plot_of_all_samples.pdf")
    })
})

test_that("plotTssPCA creates PDF", {
    data(exampleTSSr)
    tmpdir <- tempdir()
    withr::with_dir(tmpdir, {
        plotTssPCA(exampleTSSr, TSS.threshold = 10)
        expect_true(file.exists("PCA_plot.pdf"))
        unlink("PCA_plot.pdf")
    })
})

test_that("plotInterQuantile creates PDF", {
    data(exampleTSSr)
    tmpdir <- tempdir()
    if (length(exampleTSSr@clusterShape) > 0) {
        withr::with_dir(tmpdir, {
            plotInterQuantile(exampleTSSr, samples = "all", tagsThreshold = 1)
            expect_true(file.exists("Interquantile_plot_of_ALL_samples.pdf"))
            unlink("Interquantile_plot_of_ALL_samples.pdf")
        })
    }
})

test_that("plotShape creates PDF", {
    data(exampleTSSr)
    tmpdir <- tempdir()
    if (length(exampleTSSr@clusterShape) > 0) {
        withr::with_dir(tmpdir, {
            plotShape(exampleTSSr, samples = "all")
            expect_true(file.exists("Shape_plot_of_ALL_samples.pdf"))
            unlink("Shape_plot_of_ALL_samples.pdf")
        })
    }
})

test_that("plotDE creates PDF", {
    data(exampleTSSr)
    tmpdir <- tempdir()
    if (length(exampleTSSr@DEtables) > 0) {
        withr::with_dir(tmpdir, {
            plotDE(exampleTSSr, withGeneName = "TRUE")
            expect_true(file.exists("Volcano_plot.pdf"))
            unlink("Volcano_plot.pdf")
        })
    }
})
