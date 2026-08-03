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

test_that("getTSS validates input files before loading the genome", {
    obj <- TSSr(
        genomeName = "BSgenome.ThisPackageDoesNotExist",
        inputFiles = "missing_before_genome.bam",
        inputFilesType = "bam",
        sampleLabels = "S1",
        sampleLabelsMerged = "S1",
        mergeIndex = 1
    )

    expect_error(
        getTSS(obj),
        regexp = "Input file.*not found.*missing_before_genome\\.bam"
    )
})

test_that("getTSS rejects an empty input file list", {
    obj <- TSSr(genomeName = "BSgenome.ThisPackageDoesNotExist")

    expect_error(
        getTSS(obj),
        regexp = "No input files were provided"
    )
})

test_that("getTSS reports every invalid input path before doing any work", {
    input_files <- c(
        "missing_sample_1.bam",
        "missing_sample_2.bam",
        "",
        NA_character_,
        tempdir()
    )
    obj <- TSSr(
        genomeName = "BSgenome.ThisPackageDoesNotExist",
        inputFiles = input_files,
        inputFilesType = "bam",
        sampleLabels = paste0("S", seq_along(input_files)),
        sampleLabelsMerged = "merged",
        mergeIndex = rep(1, length(input_files))
    )
    before <- tssr_content(obj)

    condition <- expect_error(getTSS(obj))
    message <- conditionMessage(condition)

    expect_match(message, "missing_sample_1\\.bam")
    expect_match(message, "missing_sample_2\\.bam")
    expect_match(message, "<empty>", fixed = TRUE)
    expect_match(message, "<NA>", fixed = TRUE)
    expect_match(message, tempdir(), fixed = TRUE)
    expect_tssr_content_equal(obj, before)
})

test_that("getTSS reports only invalid paths in a mixed file list", {
    existing_file <- withr::local_tempfile(fileext = ".bam")
    expect_true(file.create(existing_file))
    missing_file <- paste0(existing_file, "-missing")
    obj <- TSSr(
        genomeName = "BSgenome.ThisPackageDoesNotExist",
        inputFiles = c(existing_file, missing_file),
        inputFilesType = "bam",
        sampleLabels = c("S1", "S2"),
        sampleLabelsMerged = "merged",
        mergeIndex = c(1, 1)
    )

    condition <- expect_error(getTSS(obj))
    message <- conditionMessage(condition)

    expect_match(message, missing_file, fixed = TRUE)
    expect_false(grepl(sQuote(existing_file), message, fixed = TRUE))
})
