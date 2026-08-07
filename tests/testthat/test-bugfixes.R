# Regression tests for package behavior fixed during earlier reviews.

test_that("consensusCluster returns clusters for a single sample", {
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

    before <- tssr_content(obj)
    result <- consensusCluster(obj, useMultiCore = FALSE)

    expect_tssr_content_equal(obj, before)
    expect_s4_class(result, "TSSr")
    expect_identical(names(result@consensusClusters), "control")
    expect_gt(nrow(result@consensusClusters$control), 0L)
    expect_true(all(
        c("cluster", "chr", "start", "end", "strand", "dominant_tss") %in%
            names(result@consensusClusters$control)
    ))
    expect_identical(
        unique(result@consensusClusters$control$strand),
        c("+", "-")
    )
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

test_that("uncoded G trimming only removes terminal mismatched Gs", {
    expect_equal(TSSr:::`.uncodedGTrimWidth`("GGG", "AAG"), 2L)
    expect_equal(TSSr:::`.uncodedGTrimWidth`("GGGG", "AAAA"), 4L)
    expect_equal(TSSr:::`.uncodedGTrimWidth`("GTA", "AAA"), 1L)
    expect_equal(TSSr:::`.uncodedGTrimWidth`("TGA", "AAA"), 0L)
    expect_equal(TSSr:::`.uncodedGTrimWidth`("GGA", "GAA"), 0L)
})

test_that("terminal G correction handles complete read batches", {
    read_seq <- c(
        "GGG", "GGGG", "GTA", "TGA", "GGA", "", "ggn", "G", "GG"
    )
    reference_seq <- c(
        "AAG", "AAAA", "AAA", "AAA", "GAA", "", "aan", "", "A"
    )

    expect_identical(
        TSSr:::`.uncodedGTrimWidths`(read_seq, reference_seq),
        c(2L, 4L, 1L, 0L, 0L, 0L, 2L, 0L, 1L)
    )
    expect_identical(
        TSSr:::`.reverseComplementText`(
            c(
                simple = "acgt", iupac = "ACGTRYSWKMBDHVN",
                sam_symbols = "=ac*", empty = ""
            )
        ),
        c("ACGT", "NBDHVKMWSRYACGT", "*GT=", "")
    )
})

test_that("uncoded G trimming updates plus starts and minus ends", {
    Genome <- BSgenome.Scerevisiae.UCSC.sacCer3::BSgenome.Scerevisiae.UCSC.sacCer3

    plus.read <- GenomicRanges::GRanges(
        "chrI", IRanges::IRanges(start = 118L, end = 140L), strand = "+"
    )
    GenomicRanges::mcols(plus.read)$seq <- "GGGACCCTCCATTACCCTGCCTC"
    plus.trimmed <- TSSr:::`.trimUncodedGOneStrand`(
        plus.read, Genome, minusStrand = FALSE
    )

    plus.non.g.read <- GenomicRanges::GRanges(
        "chrI", IRanges::IRanges(start = 118L, end = 140L), strand = "+"
    )
    GenomicRanges::mcols(plus.non.g.read)$seq <- "TTTACCCTCCATTACCCTGCCTC"
    plus.non.g.trimmed <- TSSr:::`.trimUncodedGOneStrand`(
        plus.non.g.read, Genome, minusStrand = FALSE
    )

    minus.read <- GenomicRanges::GRanges(
        "chrI", IRanges::IRanges(start = 100L, end = 120L), strand = "-"
    )
    GenomicRanges::mcols(minus.read)$seq <- "GGCCAACCTGTCTCTCAACCC"
    minus.trimmed <- TSSr:::`.trimUncodedGOneStrand`(
        minus.read, Genome, minusStrand = TRUE
    )

    expect_equal(start(plus.trimmed), 121L)
    expect_equal(end(plus.trimmed), 140L)
    expect_equal(start(plus.non.g.trimmed), 118L)
    expect_equal(end(plus.non.g.trimmed), 140L)
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
