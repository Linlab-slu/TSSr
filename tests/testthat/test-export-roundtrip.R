make_export_roundtrip_object <- function() {
    sample_labels <- paste0("SL0", 1:4)
    object <- methods::new(
        "TSSr",
        genomeName = "BSgenome.Scerevisiae.UCSC.sacCer3",
        inputFilesType = "TSStable",
        sampleLabels = sample_labels,
        sampleLabelsMerged = c("control", "treat"),
        mergeIndex = c(1, 1, 2, 2)
    )
    raw <- data.table::data.table(
        chr = rep("chrI", 3L),
        pos = as.numeric(c(100000L, 200000L, 1000000L)),
        strand = c("+", "-", "+"),
        SL01 = c(1, 0, 2),
        SL02 = c(0, 0.00025, 3),
        SL03 = c(4, 5, 6),
        SL04 = c(7, 8, 12.5)
    )
    object@TSSrawMatrix <- raw
    object@TSSprocessedMatrix <- data.table::copy(raw)
    object@librarySizes <- colSums(
        raw[, .SD, .SDcols = sample_labels]
    )
    object@normalizationStatus <- "raw"
    object
}

test_that("TSStable export and getTSS preserve coordinates and values", {
    object <- make_export_roundtrip_object()
    withr::with_tempdir({
        exportTSStable(object, data = "raw", merged = FALSE)
        exported <- file.path(getwd(), "ALL.samples.TSS.raw.txt")
        imported <- TSSr(
            genomeName = "BSgenome.Scerevisiae.UCSC.sacCer3",
            inputFiles = exported,
            inputFilesType = "TSStable",
            sampleLabels = object@sampleLabels,
            sampleLabelsMerged = object@sampleLabelsMerged,
            mergeIndex = object@mergeIndex
        )
        imported <- getTSS(imported)
        expected <- data.table::copy(object@TSSrawMatrix)
        expected[, pos := as.integer(pos)]
        data.table::setorder(expected, strand, chr, pos)

        expect_equal(imported@TSSrawMatrix, expected)
        expect_identical(typeof(imported@TSSrawMatrix$pos), "integer")
    })
})

test_that("TSStable text uses fixed notation for coordinates and values", {
    object <- make_export_roundtrip_object()
    withr::with_tempdir({
        scipen_before <- getOption("scipen")
        exportTSStable(object, data = "raw", merged = FALSE)
        text <- readLines("ALL.samples.TSS.raw.txt", warn = FALSE)

        expect_identical(getOption("scipen"), scipen_before)
        expect_false(any(grepl("[0-9][eE][+-][0-9]", text)))
        expect_true(any(grepl("100000", text, fixed = TRUE)))
        expect_true(any(grepl("0.00025", text, fixed = TRUE)))
    })
})

test_that("TSStable import accepts exponent notation and decimal values", {
    withr::with_tempdir({
        writeLines(
            c(
                "chr\tpos\tstrand\tsample",
                "chrI\t1e+05\t+\t2.5e-04"
            ),
            "external.tsv"
        )
        imported <- TSSr:::`.getTSS_from_TSStable`(
            "external.tsv", "sample"
        )

        expect_identical(imported$pos, 100000L)
        expect_identical(imported$sample, 0.00025)
    })
})

test_that("BAM-derived TSS coordinates are integers at production", {
    plus <- GenomicRanges::GRanges(
        "chrI", IRanges::IRanges(start = 100000L, width = 10L),
        strand = "+"
    )
    minus <- GenomicRanges::GRanges(
        "chrI", IRanges::IRanges(start = 199991L, width = 10L),
        strand = "-"
    )
    tss <- TSSr:::`.makeTSSFromGRanges`(plus, minus)

    expect_identical(tss$pos, c(100000L, 200000L))
    expect_identical(typeof(tss$pos), "integer")
})

test_that("all named export coordinates use integer storage", {
    table <- data.table::data.table(
        pos = 100000,
        start = 200000,
        end = 1000000,
        dominant_tss = 100000,
        dominant_tss.m = 200000,
        dominant_tss.p = 1000000,
        q_0.1 = 100000,
        q_0.9 = 200000
    )
    prepared <- TSSr:::`.prepareExportTable`(table)

    expect_true(all(vapply(prepared, is.integer, logical(1))))
    expect_true(all(vapply(table, is.double, logical(1))))
})
