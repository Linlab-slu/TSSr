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

test_that("optimized peakclu preserves the published tag clusters", {
    data(exampleTSSr)
    expected <- lapply(exampleTSSr@tagClusters, function(value) {
        as.data.frame(data.table::copy(value))
    })
    object <- exampleTSSr
    object@tagClusters <- list()
    before <- tssr_content(object)

    result <- clusterTSS(
        object,
        method = "peakclu",
        clusterThreshold = 1,
        useMultiCore = FALSE
    )

    expect_tssr_content_equal(object, before)
    observed <- lapply(result@tagClusters, function(value) {
        as.data.frame(data.table::copy(value))
    })
    expect_equal(observed, expected, tolerance = 0)
})

make_peakclu_test_object <- function(positions, tags, strand = "+") {
    object <- TSSr(
        sampleLabels = "sample",
        sampleLabelsMerged = "sample",
        mergeIndex = 1
    )
    object@TSSprocessedMatrix <- data.table::data.table(
        chr = "chrI",
        pos = as.integer(positions),
        strand = strand,
        sample = as.numeric(tags)
    )
    object
}

test_that("peakclu uses an open peak-distance window", {
    object <- make_peakclu_test_object(
        positions = c(100, 200, 201),
        tags = c(10, 20, 5)
    )

    result <- clusterTSS(
        object,
        peakDistance = 100,
        extensionDistance = 30,
        localThreshold = 0,
        clusterThreshold = 0,
        useMultiCore = FALSE
    )
    clusters <- as.data.frame(result@tagClusters$sample)

    expect_equal(clusters$start, c(100, 200))
    expect_equal(clusters$end, c(100, 200))
    expect_equal(clusters$dominant_tss, c(100, 200))
})

test_that("peakclu resolves tied window maxima to the lowest position", {
    object <- make_peakclu_test_object(
        positions = c(100, 150, 300),
        tags = c(10, 10, 5)
    )

    result <- clusterTSS(
        object,
        peakDistance = 100,
        extensionDistance = 30,
        localThreshold = 0,
        clusterThreshold = 0,
        useMultiCore = FALSE
    )
    clusters <- as.data.frame(result@tagClusters$sample)

    expect_equal(clusters$start, c(100, 300))
    expect_equal(clusters$dominant_tss, c(100, 300))
})

test_that("range maxima match brute-force first-maximum indices", {
    set.seed(20260730)

    for (iteration in seq_len(50L)) {
        size <- sample(1:500, 1L)
        values <- sample(c(1:20, rep(10, 10), rep(20, 5)),
            size,
            replace = TRUE
        )
        left <- sample(seq_len(size), 100L, replace = TRUE)
        right <- vapply(left, function(index) {
            candidates <- seq.int(index, size)
            candidates[[sample.int(length(candidates), 1L)]]
        }, integer(1))
        observed <- TSSr:::.rangeMaxFirstIndex(values)$query(left, right)
        expected <- vapply(seq_along(left), function(index) {
            window <- left[[index]]:right[[index]]
            window[[which.max(values[window])]]
        }, integer(1))

        expect_identical(observed, expected)
    }
})

test_that("peakclu multicore branch preserves single-core output", {
    object <- make_peakclu_test_object(
        positions = c(100, 150, 300, 340, 600),
        tags = c(10, 10, 5, 20, 8)
    )

    single <- clusterTSS(
        object,
        peakDistance = 100,
        extensionDistance = 30,
        localThreshold = 0.02,
        clusterThreshold = 0,
        useMultiCore = FALSE
    )
    multicore_branch <- clusterTSS(
        object,
        peakDistance = 100,
        extensionDistance = 30,
        localThreshold = 0.02,
        clusterThreshold = 0,
        useMultiCore = TRUE,
        numCores = 1
    )

    expect_equal(
        lapply(single@tagClusters, as.data.frame),
        lapply(multicore_branch@tagClusters, as.data.frame),
        tolerance = 0
    )
})
