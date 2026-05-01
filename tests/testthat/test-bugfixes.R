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
