make_example_bam_object <- function() {
    bam <- system.file(
        "extdata", "example-cigar.bam",
        package = "TSSr",
        mustWork = TRUE
    )
    TSSr(
        genomeName = "BSgenome.Scerevisiae.UCSC.sacCer3",
        inputFiles = bam,
        inputFilesType = "bam",
        sampleLabels = "sample",
        sampleLabelsMerged = "sample",
        mergeIndex = 1
    )
}

expected_example_bam_tss <- function(correct_uncoded_g) {
    plus_positions <- c(
        1000L, 1200L, 1400L, 1600L,
        if (correct_uncoded_g) 1801L else 1800L,
        1903L
    )
    minus_positions <- c(
        1109L, 1308L, 1510L, 1707L,
        if (correct_uncoded_g) 2008L else 2009L,
        2112L
    )
    data.frame(
        chr = rep("chrI", 12L),
        pos = c(plus_positions, minus_positions),
        strand = rep(c("+", "-"), each = 6L),
        sample = c(2L, rep(1L, 11L)),
        stringsAsFactors = FALSE
    )
}

test_that("getTSS imports CIGAR-aware BAM coordinates through its public API", {
    object <- make_example_bam_object()
    before <- tssr_content(object)

    aligned_boundary <- getTSS(
        object,
        sequencingQualityThreshold = 10,
        mappingQualityThreshold = 20,
        softclippingAllowed = TRUE
    )

    expect_tssr_content_equal(object, before)
    expect_equal(
        TSSmatrix(aligned_boundary, data = "raw"),
        expected_example_bam_tss(correct_uncoded_g = FALSE),
        tolerance = 0
    )
    expect_identical(librarySizes(aligned_boundary), c(sample = 13))
    expect_identical(
        typeof(TSSmatrix(aligned_boundary, data = "raw")$pos),
        "integer"
    )

    corrected <- getTSS(
        object,
        sequencingQualityThreshold = 10,
        mappingQualityThreshold = 20,
        softclippingAllowed = FALSE
    )

    expect_tssr_content_equal(object, before)
    expect_equal(
        TSSmatrix(corrected, data = "raw"),
        expected_example_bam_tss(correct_uncoded_g = TRUE),
        tolerance = 0
    )
    expect_identical(librarySizes(corrected), c(sample = 13))

    minus_positions <- TSSmatrix(corrected, data = "raw")$pos[
        TSSmatrix(corrected, data = "raw")$strand == "-"
    ]
    expect_true(all(c(1308L, 1510L) %in% minus_positions))
    expect_false(any(c(1309L, 1509L) %in% minus_positions))
})
