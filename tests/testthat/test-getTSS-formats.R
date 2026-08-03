make_import_object <- function(files, type, labels) {
    TSSr(
        genomeName = "BSgenome.Scerevisiae.UCSC.sacCer3",
        inputFiles = files,
        inputFilesType = type,
        sampleLabels = labels,
        sampleLabelsMerged = labels,
        mergeIndex = seq_along(labels)
    )
}

test_that("getTSS supplies one-to-one grouping for legacy objects", {
    input <- system.file(
        "extdata", "example-tss-table.tsv",
        package = "TSSr", mustWork = TRUE
    )
    labels <- c("SL01", "SL02", "SL03", "SL04")
    legacy <- methods::new(
        "TSSr",
        genomeName = "BSgenome.Scerevisiae.UCSC.sacCer3",
        inputFiles = input,
        inputFilesType = "TSStable",
        sampleLabels = labels
    )

    result <- getTSS(legacy)

    expect_identical(result@sampleLabelsMerged, labels)
    expect_identical(result@mergeIndex, seq_along(labels))
    expect_true(validObject(result))
})

test_that("getTSS imports BED boundaries and aggregates duplicate reads", {
    bed <- withr::local_tempfile(fileext = ".bed")
    writeLines(
        c(
            "chrI\t999\t1009\tplus_1\t0\t+",
            "chrI\t999\t1009\tplus_2\t0\t+",
            "chrI\t1099\t1109\tminus_1\t0\t-",
            "chrUnknown\t1999\t2009\tunknown\t0\t+",
            "chrI\t230217\t230228\tout_of_bounds\t0\t-"
        ),
        bed
    )
    object <- make_import_object(bed, "bed", "sample")
    before <- tssr_content(object)

    result <- getTSS(object)

    expect_tssr_content_equal(object, before)
    expect_identical(
        TSSmatrix(result, data = "raw"),
        data.frame(
            chr = c("chrI", "chrI"),
            pos = c(1000L, 1109L),
            strand = c("+", "-"),
            sample = c(2L, 1L)
        )
    )
    expect_identical(librarySizes(result), c(sample = 3))
})

test_that("getTSS merges per-sample TSS files and fills absent sites", {
    first <- withr::local_tempfile(fileext = ".tss")
    second <- withr::local_tempfile(fileext = ".tss")
    writeLines(
        c(
            "chr\tpos\tstrand\ttags",
            "chrI\t100000\t+\t2",
            "chrI\t200000\t-\t3"
        ),
        first
    )
    writeLines(
        c(
            "chr\tpos\tstrand\ttags",
            "chrI\t100000\t+\t5",
            "chrI\t300000\t-\t7"
        ),
        second
    )
    object <- make_import_object(
        c(first, second),
        "tss",
        c("sample_1", "sample_2")
    )
    before <- tssr_content(object)

    result <- getTSS(object)
    observed <- TSSmatrix(result, data = "raw")

    expect_tssr_content_equal(object, before)
    expect_equal(
        observed,
        data.frame(
            chr = rep("chrI", 3L),
            pos = c(100000L, 200000L, 300000L),
            strand = c("+", "-", "-"),
            sample_1 = c(2, 3, 0),
            sample_2 = c(5, 0, 7)
        ),
        tolerance = 0
    )
    expect_identical(typeof(observed$pos), "integer")
    expect_identical(librarySizes(result), c(sample_1 = 5, sample_2 = 12))
})

test_that("getTSS decodes BigWig score signs as TSS strands", {
    bigwig <- withr::local_tempfile(fileext = ".bw")
    track <- GenomicRanges::GRanges(
        seqnames = c("chrI", "chrI"),
        ranges = IRanges::IRanges(
            start = c(1000L, 1100L),
            end = c(1009L, 1109L)
        ),
        strand = "*",
        score = c(2, -3),
        seqinfo = GenomeInfoDb::Seqinfo(
            seqnames = "chrI",
            seqlengths = 230218L
        )
    )
    rtracklayer::export(track, bigwig, format = "BigWig")
    object <- make_import_object(bigwig, "BigWig", "sample")
    before <- tssr_content(object)

    result <- getTSS(object)
    observed <- TSSmatrix(result, data = "raw")

    expect_tssr_content_equal(object, before)
    expect_equal(
        observed,
        data.frame(
            chr = c("chrI", "chrI"),
            pos = c(1000L, 1109L),
            strand = c("+", "-"),
            sample = c(2, 3)
        ),
        tolerance = 0
    )
    expect_identical(typeof(observed$pos), "integer")
    expect_identical(librarySizes(result), c(sample = 5))
})

test_that("getTSS counts only first mates from proper BAM pairs", {
    paired_bam <- system.file(
        "extdata", "example-paired.bam",
        package = "TSSr",
        mustWork = TRUE
    )
    object <- make_import_object(paired_bam, "bamPairedEnd", "sample")
    before <- tssr_content(object)

    result <- getTSS(
        object,
        sequencingQualityThreshold = 10,
        mappingQualityThreshold = 20,
        softclippingAllowed = TRUE
    )
    observed <- TSSmatrix(result, data = "raw")

    expect_tssr_content_equal(object, before)
    expect_identical(
        observed,
        data.frame(
            chr = c("chrI", "chrI"),
            pos = c(3000L, 3209L),
            strand = c("+", "-"),
            sample = c(1L, 1L)
        )
    )
    expect_identical(librarySizes(result), c(sample = 2))
    expect_false(any(observed$pos %in% c(3059L, 3150L, 3400L, 3600L)))
})
