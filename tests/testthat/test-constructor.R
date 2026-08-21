# Test TSSr() constructor function and validation

local_input_files <- function(n, extension = ".bam",
                              .local_envir = parent.frame()) {
    input_dir <- withr::local_tempdir(.local_envir = .local_envir)
    paths <- file.path(
        input_dir,
        paste0("sample", seq_len(n), extension)
    )
    stopifnot(all(file.create(paths)))
    paths
}

test_that("TSSr() constructor creates valid object", {
    input_files <- local_input_files(2)
    obj <- TSSr(
        genomeName = "BSgenome.Scerevisiae.UCSC.sacCer3",
        inputFiles = input_files,
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

test_that("TSSr() defaults omitted grouping to one group per sample", {
    input_files <- local_input_files(2)
    obj <- TSSr(
        genomeName = "BSgenome.Scerevisiae.UCSC.sacCer3",
        inputFiles = input_files,
        inputFilesType = "bam",
        sampleLabels = c("S1", "S2")
    )

    expect_equal(obj@sampleLabelsMerged, c("S1", "S2"))
    expect_equal(obj@mergeIndex, c(1, 2))
})

test_that("TSSr() requires core input arguments", {
    expect_error(
        TSSr(),
        regexp = "is missing, with no default"
    )
})

test_that("TSSr() rejects empty required metadata", {
    input_file <- local_input_files(1)

    expect_error(
        TSSr(
            genomeName = character(),
            inputFiles = input_file,
            inputFilesType = "bam",
            sampleLabels = "S1"
        ),
        regexp = "genomeName.*non-empty character string"
    )
    expect_error(
        TSSr(
            genomeName = "BSgenome.Scerevisiae.UCSC.sacCer3",
            inputFiles = input_file,
            inputFilesType = character(),
            sampleLabels = "S1"
        ),
        regexp = "inputFilesType.*non-empty character string"
    )
    expect_error(
        TSSr(
            genomeName = "BSgenome.Scerevisiae.UCSC.sacCer3",
            inputFiles = input_file,
            inputFilesType = "bam",
            sampleLabels = character()
        ),
        regexp = "sampleLabels.*non-empty character"
    )
})

test_that("an empty TSSr object can be created for internal use", {
    obj <- methods::new("TSSr")
    expect_s4_class(obj, "TSSr")
    expect_equal(length(obj@inputFiles), 0)
    expect_equal(length(obj@sampleLabels), 0)
    expect_identical(obj@normalizationStatus, NA_character_)
})

test_that("TSSr() rejects invalid inputFilesType", {
    input_file <- local_input_files(1)
    expect_error(
        TSSr(
            genomeName = "BSgenome.Scerevisiae.UCSC.sacCer3",
            inputFiles = input_file,
            inputFilesType = "invalid_type",
            sampleLabels = "S1"
        ),
        regexp = "inputFilesType.*must be one of"
    )
})

test_that("TSSr() rejects invalid input file paths at construction", {
    missing_file <- file.path(tempdir(), "tssr-file-that-does-not-exist.bam")
    input_directory <- withr::local_tempdir()
    invalid_paths <- list(missing_file, "", NA_character_, input_directory)

    for (path in invalid_paths) {
        expect_error(
            TSSr(
                genomeName = "BSgenome.Scerevisiae.UCSC.sacCer3",
                inputFiles = path,
                inputFilesType = "bam",
                sampleLabels = "S1"
            ),
            regexp = "inputFiles.*existing regular file"
        )
    }
})

test_that("TSSr() validates refSource only when it is provided", {
    input_file <- local_input_files(1)
    annotation_file <- local_input_files(1, extension = ".gff3")
    common <- list(
        genomeName = "BSgenome.Scerevisiae.UCSC.sacCer3",
        inputFiles = input_file,
        inputFilesType = "bam",
        sampleLabels = "S1"
    )

    without_reference <- do.call(TSSr, common)
    expect_length(without_reference@refSource, 0L)

    with_reference <- do.call(
        TSSr,
        c(common, list(refSource = annotation_file))
    )
    expect_identical(with_reference@refSource, annotation_file)

    expect_error(
        do.call(
            TSSr,
            c(common, list(refSource = paste0(annotation_file, ".missing")))
        ),
        regexp = "refSource.*existing regular file"
    )
})

test_that("TSSr() requires grouping arguments as a pair", {
    input_file <- local_input_files(1)

    expect_error(
        TSSr(
            genomeName = "BSgenome.Scerevisiae.UCSC.sacCer3",
            inputFiles = input_file,
            inputFilesType = "bam",
            sampleLabels = "S1",
            sampleLabelsMerged = "group"
        ),
        regexp = "sampleLabelsMerged.*mergeIndex.*provided together"
    )
})

test_that("TSSr() rejects mismatched inputFiles and sampleLabels lengths", {
    input_files <- local_input_files(3)
    expect_error(
        TSSr(
            genomeName = "BSgenome.Scerevisiae.UCSC.sacCer3",
            inputFiles = input_files,
            inputFilesType = "bam",
            sampleLabels = c("S1", "S2")
        ),
        regexp = "inputFiles.*sampleLabels"
    )
})

test_that("TSSr() accepts exactly one TSStable input table", {
    input_files <- local_input_files(2, extension = ".tsv")

    expect_error(
        TSSr(
            genomeName = "BSgenome.Scerevisiae.UCSC.sacCer3",
            inputFiles = input_files,
            inputFilesType = "TSStable",
            sampleLabels = c("S1", "S2")
        ),
        regexp = "Exactly one.*TSStable"
    )
})

test_that("TSSr() rejects mismatched mergeIndex and sampleLabelsMerged", {
    input_files <- local_input_files(2)
    expect_error(
        TSSr(
            genomeName = "BSgenome.Scerevisiae.UCSC.sacCer3",
            inputFiles = input_files,
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
        input_file <- local_input_files(1)
        obj <- TSSr(
            genomeName = "BSgenome.Scerevisiae.UCSC.sacCer3",
            inputFiles = input_file,
            inputFilesType = type,
            sampleLabels = "S1"
        )
        expect_s4_class(obj, "TSSr")
        expect_equal(obj@inputFilesType, type)
    }
})
