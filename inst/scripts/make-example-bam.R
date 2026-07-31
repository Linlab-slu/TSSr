#!/usr/bin/env Rscript

# Rebuild the small BAM fixture used by the public getTSS() regression test.
# Every read is based on the real sacCer3 chromosome-I sequence. The fixture
# intentionally has no BAM index because getTSS() scans the complete file.

script_argument <- grep("^--file=", commandArgs(), value = TRUE)
if (length(script_argument) == 1L) {
    script_path <- sub("^--file=", "", script_argument)
    if (!grepl("^/", script_path)) {
        script_path <- file.path(Sys.getenv("PWD"), script_path)
    }
    script_path <- normalizePath(script_path)
    setwd(dirname(dirname(dirname(script_path))))
}

genome <- getExportedValue(
    "BSgenome.Scerevisiae.UCSC.sacCer3",
    "BSgenome.Scerevisiae.UCSC.sacCer3"
)
chromosome <- as.character(genome[["chrI"]])

reference_sequence <- function(start, end) {
    substr(chromosome, start, end)
}

insert_base <- function(start, base = "A") {
    paste0(
        reference_sequence(start, start + 4L),
        base,
        reference_sequence(start + 5L, start + 8L)
    )
}

delete_base <- function(start) {
    paste0(
        reference_sequence(start, start + 4L),
        reference_sequence(start + 6L, start + 10L)
    )
}

replace_last_base <- function(sequence, base) {
    paste0(substr(sequence, 1L, nchar(sequence) - 1L), base)
}

# Fixed reference checks make accidental genome/coordinate drift fail loudly.
stopifnot(
    reference_sequence(1800L, 1800L) == "T",
    reference_sequence(1903L, 1903L) == "G",
    reference_sequence(2009L, 2009L) == "T",
    reference_sequence(2112L, 2112L) == "C"
)

records <- data.frame(
    qname = c(
        "plus_plain_1", "plus_plain_2", "minus_plain",
        "plus_insertion", "minus_insertion",
        "plus_deletion", "minus_deletion",
        "plus_softclip", "minus_softclip",
        "plus_uncoded_g", "plus_matched_g",
        "minus_uncoded_g", "minus_matched_g",
        "low_mapq", "low_quality"
    ),
    flag = c(0L, 0L, 16L, 0L, 16L, 0L, 16L, 0L, 16L,
             0L, 0L, 16L, 16L, 0L, 0L),
    pos = c(1000L, 1000L, 1100L, 1200L, 1300L, 1400L, 1500L,
            1600L, 1700L, 1800L, 1903L, 2000L, 2103L, 2200L, 2300L),
    mapq = c(rep(60L, 13L), 5L, 60L),
    cigar = c(
        "10M", "10M", "10M", "5M1I4M", "5M1I4M",
        "5M1D5M", "5M1D5M", "2S8M", "8M2S",
        "10M", "10M", "10M", "10M", "10M", "10M"
    ),
    purpose = c(
        "plain plus and count aggregation",
        "plain plus and count aggregation",
        "plain minus",
        "plus insertion",
        "minus insertion distinguishes reference from query width",
        "plus deletion",
        "minus deletion distinguishes reference from query width",
        "plus soft clipping",
        "minus soft clipping",
        "plus terminal mismatched G",
        "plus terminal reference-matched G",
        "minus terminal mismatched transcript-sense G",
        "minus terminal reference-matched transcript-sense G",
        "MAPQ filtering",
        "sequencing-quality filtering"
    ),
    stringsAsFactors = FALSE
)

records$seq <- c(
    reference_sequence(1000L, 1009L),
    reference_sequence(1000L, 1009L),
    reference_sequence(1100L, 1109L),
    insert_base(1200L),
    insert_base(1300L),
    delete_base(1400L),
    delete_base(1500L),
    paste0("AA", reference_sequence(1600L, 1607L)),
    paste0(reference_sequence(1700L, 1707L), "AA"),
    paste0("G", reference_sequence(1801L, 1809L)),
    reference_sequence(1903L, 1912L),
    replace_last_base(reference_sequence(2000L, 2009L), "C"),
    reference_sequence(2103L, 2112L),
    reference_sequence(2200L, 2209L),
    reference_sequence(2300L, 2309L)
)
records$qual <- vapply(records$seq, function(sequence) {
    strrep("I", nchar(sequence))
}, character(1))
records$qual[records$qname == "low_quality"] <- strrep("!", 10L)

stopifnot(
    nrow(records) == 15L,
    all(nchar(records$seq) == 10L),
    all(nchar(records$qual) == nchar(records$seq))
)

sam_lines <- c(
    "@HD\tVN:1.6\tSO:unsorted",
    paste0(
        "@SQ\tSN:chrI\tLN:",
        GenomeInfoDb::seqlengths(genome)[["chrI"]]
    ),
    vapply(seq_len(nrow(records)), function(index) {
        paste(
            records$qname[index], records$flag[index], "chrI",
            records$pos[index], records$mapq[index], records$cigar[index],
            "*", 0L, 0L, records$seq[index], records$qual[index],
            sep = "\t"
        )
    }, character(1))
)

sam_path <- tempfile(fileext = ".sam")
on.exit(unlink(sam_path), add = TRUE)
writeLines(sam_lines, sam_path, useBytes = TRUE)

destination <- file.path("inst", "extdata", "example-cigar")
bam_path <- Rsamtools::asBam(
    sam_path,
    destination = destination,
    overwrite = TRUE,
    indexDestination = FALSE
)

stopifnot(
    identical(normalizePath(bam_path), normalizePath(paste0(destination, ".bam"))),
    !file.exists(paste0(destination, ".bam.bai")),
    length(Rsamtools::scanBam(bam_path)[[1L]]$qname) == 15L
)

message("Wrote ", bam_path, " with 15 designed reads and no BAM index.")
