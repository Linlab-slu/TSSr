# Test print->plot and print->show changes in export methods

test_that("plotTssPCA creates PDF using plot() instead of print()", {
    data(exampleTSSr)
    # plotTssPCA was changed from print(autoplot(...)) to plot(autoplot(...))
    # It creates a PCA_plot.pdf in the working directory
    tmpdir <- tempdir()
    old_wd <- setwd(tmpdir)
    on.exit(setwd(old_wd))

    plotTssPCA(exampleTSSr)
    expect_true(file.exists("PCA_plot.pdf"))
    expect_gt(file.size("PCA_plot.pdf"), 0)
    unlink("PCA_plot.pdf")
})

test_that("exportClustersToBed runs without print() calls", {
    data(exampleTSSr)
    # ExportFunctions.R:149 changed print() to show() in a debug path
    # exportClustersToBed uses that code path
    tmpdir <- tempdir()
    old_wd <- setwd(tmpdir)
    on.exit(setwd(old_wd))

    result <- consensusCluster(exampleTSSr, useMultiCore = FALSE)
    # This exercises the BED export code including the show() path
    expect_no_error(
        exportClustersToBed(result, data = "consensusClusters")
    )
})
