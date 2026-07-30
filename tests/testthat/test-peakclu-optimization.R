test_that("clusterTSS rejects an empty peak window clearly", {
    data(exampleTSSr)

    expect_error(
        clusterTSS(
            exampleTSSr,
            method = "peakclu",
            peakDistance = 0,
            clusterThreshold = 1,
            useMultiCore = FALSE
        ),
        "peakDistance must be greater than zero"
    )
})
