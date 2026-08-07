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
    for (sample_label in names(observed)) {
        exact_columns <- setdiff(names(observed[[sample_label]]), "tags")
        expect_identical(
            observed[[sample_label]][, exact_columns],
            expected[[sample_label]][, exact_columns],
            info = sample_label
        )
        expect_equal(
            observed[[sample_label]]$tags,
            expected[[sample_label]]$tags,
            tolerance = 1e-12,
            info = sample_label
        )
    }
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

test_that("peakclu includes a TSS exactly on the quantile threshold", {
    object <- make_peakclu_test_object(
        positions = c(100, 103, 106),
        tags = c(81, 9, 1)
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

    expect_equal(clusters$q_0.9, 103)
    expect_equal(clusters$interquantile_width, 4)
})

test_that("peakclu uses the same quantiles after positive signal scaling", {
    scales <- c(1e-6, 1, 1e6)

    boundaries <- lapply(scales, function(scale) {
        object <- make_peakclu_test_object(
            positions = c(100, 103, 106),
            tags = c(9, 81, 1) * scale
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

        unname(c(clusters$q_0.1, clusters$q_0.9))
    })

    expect_identical(boundaries, rep(list(c(100L, 103L)), length(scales)))
})

test_that("consensus clusters include a TSS on the quantile threshold", {
    object <- make_peakclu_test_object(
        positions = c(100, 103, 106),
        tags = c(81, 9, 1)
    )
    object <- clusterTSS(
        object,
        peakDistance = 100,
        extensionDistance = 30,
        localThreshold = 0,
        clusterThreshold = 0,
        useMultiCore = FALSE
    )

    result <- consensusCluster(object, useMultiCore = FALSE)
    clusters <- as.data.frame(result@consensusClusters$sample)

    expect_equal(clusters$q_0.9, 103)
    expect_equal(clusters$interquantile_width, 4)
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

test_that("cluster assembly preserves directional local filtering", {
    expected <- list(
        `+` = data.frame(
            start = 100L, end = 100L, tags = 10,
            interquantile_width = 1
        ),
        `-` = data.frame(
            start = 100L, end = 120L, tags = 10.1,
            interquantile_width = 1
        )
    )

    for (strand in names(expected)) {
        object <- make_peakclu_test_object(
            positions = c(100, 120, 140),
            tags = c(10, 0.1, 1),
            strand = strand
        )
        result <- clusterTSS(
            object,
            method = "peakcluMax",
            peakDistance = 100,
            extensionDistance = 30,
            localThreshold = 0.02,
            clusterThreshold = 0,
            useMultiCore = FALSE
        )
        observed <- as.data.frame(result@tagClusters$sample)
        columns <- names(expected[[strand]])

        expect_identical(
            observed[, columns, drop = FALSE],
            expected[[strand]],
            info = strand
        )
    }
})
