# Test TSSr() constructor function and validation

test_that("TSSr() constructor creates valid object", {
    obj <- TSSr(
        genomeName = "BSgenome.Scerevisiae.UCSC.sacCer3",
        inputFiles = c("s1.bam", "s2.bam"),
        inputFilesType = "bam",
        sampleLabels = c("S1", "S2"),
        sampleLabelsMerged = c("merged"),
        mergeIndex = c(1, 1)
    )
    expect_s4_class(obj, "TSSr")
    expect_equal(obj@genomeName, "BSgenome.Scerevisiae.UCSC.sacCer3")
    expect_equal(obj@inputFilesType, "bam")
    expect_equal(obj@sampleLabels, c("S1", "S2"))
})

test_that("TSSr() defaults to one group per sample when merging is omitted", {
    obj <- TSSr(
        genomeName = "BSgenome.Scerevisiae.UCSC.sacCer3",
        inputFiles = c("s1.bam", "s2.bam", "s3.bam"),
        inputFilesType = "bam",
        sampleLabels = c("S1", "S2", "S3")
    )

    expect_identical(obj@sampleLabelsMerged, c("S1", "S2", "S3"))
    expect_identical(obj@mergeIndex, 1:3)
    expect_true(validObject(obj))
})

test_that("TSSr() rejects partially specified sample grouping", {
    common_args <- list(
        inputFiles = c("s1.bam", "s2.bam"),
        inputFilesType = "bam",
        sampleLabels = c("S1", "S2")
    )

    condition <- expect_error(
        do.call(TSSr, c(common_args, list(sampleLabelsMerged = "merged"))),
        regexp = "sampleLabelsMerged.*mergeIndex.*both"
    )
    message <- conditionMessage(condition)
    expect_match(message, "sampleLabelsMerged.*sampleLabels")
    expect_match(message, "mergeIndex.*seq_along\\(sampleLabels\\)")
    expect_error(
        do.call(TSSr, c(common_args, list(mergeIndex = c(1, 1)))),
        regexp = "sampleLabelsMerged.*mergeIndex.*both"
    )
})

test_that("TSSr class validity rejects partially specified sample grouping", {
    condition <- expect_error(
        methods::new(
            "TSSr",
            sampleLabels = c("S1", "S2"),
            sampleLabelsMerged = "merged"
        ),
        regexp = "sampleLabelsMerged.*mergeIndex.*both"
    )
    expect_match(
        conditionMessage(condition),
        "sampleLabelsMerged.*sampleLabels.*mergeIndex.*seq_along"
    )
})

test_that("TSSr() constructor with empty arguments creates valid object", {
    obj <- TSSr()
    expect_s4_class(obj, "TSSr")
    expect_equal(length(obj@inputFiles), 0)
    expect_equal(length(obj@sampleLabels), 0)
    expect_identical(obj@normalizationStatus, NA_character_)
})

test_that("TSSr() rejects invalid inputFilesType", {
    expect_error(
        TSSr(inputFilesType = "invalid_type"),
        regexp = "inputFilesType.*must be one of"
    )
})

test_that("TSSr() rejects mismatched inputFiles and sampleLabels lengths", {
    expect_error(
        TSSr(
            inputFiles = c("a.bam", "b.bam", "c.bam"),
            inputFilesType = "bam",
            sampleLabels = c("S1", "S2")
        ),
        regexp = "inputFiles.*sampleLabels"
    )
})

test_that("TSSr() rejects mismatched mergeIndex and sampleLabelsMerged", {
    expect_error(
        TSSr(
            inputFiles = c("a.bam", "b.bam"),
            inputFilesType = "bam",
            sampleLabels = c("S1", "S2"),
            sampleLabelsMerged = c("G1", "G2", "G3"),
            mergeIndex = c(1, 1)
        ),
        regexp = "mergeIndex.*sampleLabelsMerged"
    )
})

test_that("TSSr() accepts all valid inputFilesType values", {
    valid_types <- c("bam", "bamPairedEnd", "bed", "tss", "TSStable", "BigWig")
    for (type in valid_types) {
        obj <- TSSr(inputFilesType = type)
        expect_s4_class(obj, "TSSr")
        expect_equal(obj@inputFilesType, type)
    }
})
