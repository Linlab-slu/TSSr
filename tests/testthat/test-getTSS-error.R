# Test getTSS error handling for input files removed after construction

test_that("getTSS revalidates input files before importing", {
    input_file <- withr::local_tempfile(
        pattern = "removed-input-",
        fileext = ".tsv",
        lines = ""
    )
    obj <- TSSr(
        genomeName = "BSgenome.DoesNotExist",
        inputFiles = input_file,
        inputFilesType = "TSStable",
        sampleLabels = "S1"
    )
    unlink(input_file)

    expect_error(
        getTSS(obj),
        regexp = "inputFiles.*existing regular file.*removed-input"
    )
})
