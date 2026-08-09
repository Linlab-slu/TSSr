test_that("consensusCluster preserves the bundled consensus results", {
    data(exampleTSSr)
    expected <- lapply(exampleTSSr@consensusClusters, as.data.frame)
    object <- exampleTSSr
    object@consensusClusters <- list()
    original <- tssr_content(object)

    result <- consensusCluster(object, useMultiCore = FALSE)
    observed <- lapply(result@consensusClusters, as.data.frame)

    expect_equal(observed, expected, tolerance = 1e-12)
    expect_tssr_content_equal(object, original)
})

test_that("consensusCluster serial and multicore results agree", {
    skip_on_os("windows")
    data(exampleTSSr)
    object <- exampleTSSr
    object@consensusClusters <- list()

    serial <- consensusCluster(object, useMultiCore = FALSE)
    parallel <- consensusCluster(
        object,
        useMultiCore = TRUE,
        numCores = 2
    )

    expect_equal(
        lapply(parallel@consensusClusters, as.data.frame),
        lapply(serial@consensusClusters, as.data.frame),
        tolerance = 1e-12
    )
})
