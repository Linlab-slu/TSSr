make_shape_interval_object <- function() {
    object <- methods::new(
        "TSSr",
        sampleLabels = "S1",
        sampleLabelsMerged = "S1",
        mergeIndex = 1,
        normalizationStatus = "normalized"
    )
    object@TSSprocessedMatrix <- data.table::data.table(
        chr = c("chrI", "chrI", "chrI", "chrI", "chrII"),
        pos = c(10L, 11L, 12L, 11L, 11L),
        strand = c("+", "+", "+", "-", "+"),
        S1 = c(1, 3, 6, 50, 100)
    )
    object@consensusClusters <- list(
        S1 = data.table::data.table(
            cluster = 1:2,
            chr = c("chrI", "chrI"),
            start = c(10L, 11L),
            end = c(11L, 12L),
            strand = c("+", "+"),
            dominant_tss = c(11L, 12L),
            tags = c(4, 9),
            tags.dominant_tss = c(3, 6),
            q_0.1 = c(10L, 11L),
            q_0.9 = c(11L, 12L),
            interquantile_width = c(2L, 2L)
        )
    )
    object
}

shape_contribution <- function(tags) {
    proportions <- tags / sum(tags)
    sum(proportions * log(proportions, 2))
}

test_that("shapeCluster includes shared points in every overlapping interval", {
    object <- make_shape_interval_object()
    original <- tssr_content(object)

    pss <- shapeCluster(
        object,
        clusters = "consensusClusters",
        method = "PSS",
        useMultiCore = FALSE
    )
    si <- shapeCluster(
        object,
        clusters = "consensusClusters",
        method = "SI",
        useMultiCore = FALSE
    )

    expected.contributions <- c(
        shape_contribution(c(1, 3)),
        shape_contribution(c(3, 6))
    )
    expect_equal(
        pss@clusterShape$S1$shape.score,
        -expected.contributions * log(c(2, 2), 2),
        tolerance = 1e-12
    )
    expect_equal(
        si@clusterShape$S1$shape.score,
        2 + expected.contributions,
        tolerance = 1e-12
    )
    expect_identical(pss@clusterShape$S1$cluster, 1:2)
    expect_tssr_content_equal(object, original)
})

test_that("shapeCluster preserves the established empty-cluster result", {
    object <- make_shape_interval_object()
    object@consensusClusters$S1 <- object@consensusClusters$S1[0]

    result <- shapeCluster(
        object,
        clusters = "consensusClusters",
        method = "PSS",
        useMultiCore = FALSE
    )

    expect_s3_class(result@clusterShape$S1, "data.table")
    expect_identical(dim(result@clusterShape$S1), c(0L, 0L))
})
