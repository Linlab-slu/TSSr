make_peakclu_max_test_object <- function(positions, tags, strand = "+") {
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

test_that("peakcluMax keeps the next strongest unsuppressed peak", {
    object <- make_peakclu_max_test_object(
        positions = c(1, 91, 181),
        tags = c(10, 9, 8)
    )

    result <- clusterTSS(
        object,
        method = "peakcluMax",
        peakDistance = 100,
        extensionDistance = 30,
        localThreshold = 0,
        clusterThreshold = 0,
        useMultiCore = FALSE
    )
    clusters <- as.data.frame(result@tagClusters$sample)

    expect_identical(clusters$start, c(1L, 181L))
    expect_identical(clusters$end, c(1L, 181L))
    expect_identical(clusters$dominant_tss, c(1L, 181L))
})

test_that("peakcluMax suppresses positions exactly peakDistance away", {
    object <- make_peakclu_max_test_object(
        positions = c(100, 200, 301),
        tags = c(10, 9, 8)
    )

    result <- clusterTSS(
        object,
        method = "peakcluMax",
        peakDistance = 100,
        extensionDistance = 30,
        localThreshold = 0,
        clusterThreshold = 0,
        useMultiCore = FALSE
    )
    clusters <- as.data.frame(result@tagClusters$sample)

    expect_identical(clusters$dominant_tss, c(100L, 301L))
})

test_that("peakcluMax resolves equal signals to the lowest position", {
    object <- make_peakclu_max_test_object(
        positions = c(100, 150, 300),
        tags = c(10, 10, 5)
    )

    result <- clusterTSS(
        object,
        method = "peakcluMax",
        peakDistance = 100,
        extensionDistance = 30,
        localThreshold = 0,
        clusterThreshold = 0,
        useMultiCore = FALSE
    )
    clusters <- as.data.frame(result@tagClusters$sample)

    expect_identical(clusters$dominant_tss, c(100L, 300L))
})

test_that("peakcluMax selector matches a brute-force greedy reference", {
    brute_force <- function(positions, tags, peak_distance) {
        priority <- order(-tags, positions)
        suppressed <- rep(FALSE, length(positions))
        selected <- integer()
        for (index in priority) {
            if (suppressed[[index]]) next
            selected <- c(selected, index)
            suppressed[
                positions >= positions[[index]] - peak_distance &
                    positions <= positions[[index]] + peak_distance
            ] <- TRUE
        }
        sort(selected)
    }

    set.seed(20260806)
    for (iteration in seq_len(100L)) {
        size <- sample(1:500, 1L)
        positions <- sort(sample.int(size * 10L, size))
        tags <- sample(c(1:20, rep(10, 10), rep(20, 5)),
            size,
            replace = TRUE
        )
        peak_distance <- sample(1:200, 1L)

        expect_identical(
            TSSr:::.selectPeakMaxIndices(
                positions, tags, peak_distance
            ),
            brute_force(positions, tags, peak_distance),
            info = paste("iteration", iteration)
        )
    }
})

test_that("peakcluMax multicore branch preserves single-core output", {
    object <- make_peakclu_max_test_object(
        positions = c(1, 91, 181, 300, 340, 600),
        tags = c(10, 9, 8, 5, 20, 8)
    )

    single <- clusterTSS(
        object,
        method = "peakcluMax",
        peakDistance = 100,
        extensionDistance = 30,
        localThreshold = 0.02,
        clusterThreshold = 0,
        useMultiCore = FALSE
    )
    multicore_branch <- clusterTSS(
        object,
        method = "peakcluMax",
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

test_that("clusterTSS rejects unknown clustering methods clearly", {
    object <- make_peakclu_max_test_object(positions = 100, tags = 10)

    expect_error(
        clusterTSS(object, method = "unknown"),
        "must be one of.*peakclu.*peakcluMax"
    )
})
