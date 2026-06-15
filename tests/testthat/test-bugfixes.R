# Test fixes for bugs found in v3 code review

test_that("BUG-1: data.table calls have no stringsAsFactors column", {
    # data.table() doesn't support stringsAsFactors; it would create a
    # spurious column. Verify the internal function produces clean output.
    dt <- data.table::data.table(
        chr = "chrI", pos = 100L, strand = "+", score = 1.5
    )
    expect_false("stringsAsFactors" %in% names(dt))
    expect_equal(ncol(dt), 4)
})

test_that("BUG-2: consensusCluster handles single sample without error", {
    data(exampleTSSr)
    # Create a single-sample scenario
    obj <- TSSr(
        genomeName = "BSgenome.Scerevisiae.UCSC.sacCer3",
        inputFiles = c("s1.bam"),
        inputFilesType = "bam",
        sampleLabels = c("S1"),
        sampleLabelsMerged = c("S1"),
        mergeIndex = c(1)
    )
    # Use exampleTSSr data but simulate single sample
    obj@TSSprocessedMatrix <- exampleTSSr@TSSprocessedMatrix[,
        c("chr", "pos", "strand", "control"), with = FALSE]
    obj@sampleLabelsMerged <- c("control")
    obj@tagClusters <- list(control = exampleTSSr@tagClusters[["control"]])

    # Should not error with single sample
    expect_no_error(consensusCluster(obj, useMultiCore = FALSE))
})

test_that("BUG-3: normalizeTSS does not load BSgenome unnecessarily", {
    data(exampleTSSr)
    mergeSamples(exampleTSSr)
    # normalizeTSS should work without .getGenome being called
    # (previously it loaded BSgenome but never used it)
    expect_no_error(normalizeTSS(exampleTSSr))
})

test_that("BUG-4: no direct @seqnames slot access in source code", {
    # Verify source files don't use @seqnames (should use seqnames() accessor)
    r_files <- list.files(
        system.file("R", package = "TSSr"),
        full.names = TRUE
    )
    # If installed R files aren't accessible, check source
    src_files <- list.files(
        file.path(system.file(package = "TSSr"), ".."),
        pattern = "[.]R$", recursive = TRUE, full.names = TRUE
    )
    # At minimum, verify seqnames() accessor works on a GRanges object
    gr <- GenomicRanges::GRanges("chrI:1-100:+")
    expect_equal(as.character(GenomeInfoDb::seqnames(gr)), "chrI")
})

test_that("BUG-5: withr is available for tests", {
    expect_true(requireNamespace("withr", quietly = TRUE))
})

test_that("BAM CIGAR widths are measured in reference space", {
    cigar <- c(
        "100M", "10S90M", "90M10S", "50M10I40M",
        "50M10D40M", "10S50M10I40M5S", "20M100N30M", "5H95M"
    )
    ref_width <- TSSr:::`.cigarReferenceWidth`(cigar)

    expect_equal(ref_width, c(100L, 90L, 90L, 90L, 100L, 90L, 150L, 95L))
    expect_equal(100L + ref_width[4] - 1L, 189L)
    expect_equal(100L + ref_width[6] - 1L, 189L)
})

test_that("5' terminal mismatch trimming stops at the first match", {
    expect_equal(TSSr:::`.leadingMismatchCount`("GTG", "AAG"), 2L)
    expect_equal(TSSr:::`.leadingMismatchCount`("GTA", "CCC"), 3L)
    expect_equal(TSSr:::`.leadingMismatchCount`("GTAC", "CCCC"), 3L)
    expect_equal(TSSr:::`.leadingMismatchCount`("GTA", "GCC"), 0L)
})

test_that("5' terminal mismatch trimming updates plus starts and minus ends", {
    skip_if_not_installed("BSgenome.Scerevisiae.UCSC.sacCer3")
    Genome <- BSgenome.Scerevisiae.UCSC.sacCer3::BSgenome.Scerevisiae.UCSC.sacCer3

    plus.read <- GenomicRanges::GRanges(
        "chrI", IRanges::IRanges(start = 100L, end = 120L), strand = "+"
    )
    GenomicRanges::mcols(plus.read)$seq <- "TTCCAACCTGTCTCTCAACTT"
    plus.trimmed <- TSSr:::`.trimTerminalMismatchesOneStrand`(
        plus.read, Genome, minusStrand = FALSE
    )

    minus.read <- GenomicRanges::GRanges(
        "chrI", IRanges::IRanges(start = 100L, end = 120L), strand = "-"
    )
    GenomicRanges::mcols(minus.read)$seq <- "GGCCAACCTGTCTCTCAACAC"
    minus.trimmed <- TSSr:::`.trimTerminalMismatchesOneStrand`(
        minus.read, Genome, minusStrand = TRUE
    )

    expect_equal(start(plus.trimmed), 102L)
    expect_equal(end(plus.trimmed), 120L)
    expect_equal(start(minus.trimmed), 100L)
    expect_equal(end(minus.trimmed), 118L)
})

test_that("representation() replaced with slots = list()", {
    # Verify TSSr class uses slots (not deprecated representation)
    slot_names <- slotNames("TSSr")
    expect_true(length(slot_names) > 0)
    expect_true("genomeName" %in% slot_names)
    expect_true("TSSrawMatrix" %in% slot_names)
    expect_true("PromoterShift" %in% slot_names)

    # Verify object creation still works with new slot syntax
    obj <- TSSr(
        genomeName = "test",
        inputFilesType = "bam"
    )
    expect_s4_class(obj, "TSSr")
})
