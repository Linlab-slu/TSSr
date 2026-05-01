# Test getTSS error handling for missing input files

test_that("getTSS gives informative error when files don't exist", {
    obj <- TSSr(
        genomeName = "BSgenome.Scerevisiae.UCSC.sacCer3",
        inputFiles = c("nonexistent_file.bam"),
        inputFilesType = "bam",
        sampleLabels = c("S1"),
        sampleLabelsMerged = c("S1"),
        mergeIndex = c(1)
    )
    expect_error(
        getTSS(obj),
        regexp = "Input file.*not found"
    )
})

test_that("getTSS error message includes the missing file name", {
    obj <- TSSr(
        genomeName = "BSgenome.Scerevisiae.UCSC.sacCer3",
        inputFiles = c("missing_sample.bam"),
        inputFilesType = "bam",
        sampleLabels = c("S1"),
        sampleLabelsMerged = c("S1"),
        mergeIndex = c(1)
    )
    expect_error(
        getTSS(obj),
        regexp = "missing_sample\\.bam"
    )
})
