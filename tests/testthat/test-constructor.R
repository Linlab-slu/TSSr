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

test_that("TSSr() constructor with empty arguments creates valid object", {
    obj <- TSSr()
    expect_s4_class(obj, "TSSr")
    expect_equal(length(obj@inputFiles), 0)
    expect_equal(length(obj@sampleLabels), 0)
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
