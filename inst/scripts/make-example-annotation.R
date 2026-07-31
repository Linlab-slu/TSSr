#!/usr/bin/env Rscript

# Rebuild the compact GFF3 fixture used by annotateCluster() tests. The source
# is the same SGD R64-2-1 annotation distributed with the original TSSr
# workflow. The fixture is a contiguous chrI slice containing the first 12
# complete gene records and all accompanying features in that interval.

script_argument <- grep("^--file=", commandArgs(), value = TRUE)
if (length(script_argument) == 1L) {
    script_path <- sub("^--file=", "", script_argument)
    if (!grepl("^/", script_path)) {
        script_path <- file.path(Sys.getenv("PWD"), script_path)
    }
    script_path <- normalizePath(script_path)
    setwd(dirname(dirname(dirname(script_path))))
}

source_url <- paste0(
    "https://www.zlinlab.org/TSSr/",
    "saccharomyces_cerevisiae_R64-2-1.gff"
)
expected_md5 <- "73f6441bcfcd519ef4bf2005a25670de"
source_path <- Sys.getenv("TSSR_SOURCE_GFF", unset = "")

if (!nzchar(source_path)) {
    source_path <- tempfile(fileext = ".gff")
    on.exit(unlink(source_path), add = TRUE)
    download.file(source_url, source_path, mode = "wb", quiet = TRUE)
}
source_path <- normalizePath(source_path, mustWork = TRUE)

stopifnot(
    identical(unname(tools::md5sum(source_path)), expected_md5)
)

source_lines <- readLines(source_path, warn = FALSE)
first_record <- which(!grepl("^#", source_lines))[1L]
header <- source_lines[seq_len(first_record - 1L)]
record_lines <- source_lines[-seq_len(first_record - 1L)]
fields <- strsplit(record_lines, "\t", fixed = TRUE)
is_gff_record <- lengths(fields) == 9L
fields <- fields[is_gff_record]
record_lines <- record_lines[is_gff_record]

seqname <- vapply(fields, `[[`, character(1), 1L)
feature_type <- vapply(fields, `[[`, character(1), 3L)
feature_end <- as.integer(vapply(fields, `[[`, character(1), 5L))

# The twelfth chrI gene is YAL063C and ends at 27,968. Retain the chromosome
# declaration plus every original chrI feature wholly contained in this
# contiguous interval, so parent/child and neighboring feature context remain
# exactly as distributed by SGD.
cutoff <- 27968L
is_chrI_declaration <- seqname == "chrI" & feature_type == "chromosome"
in_slice <- seqname == "chrI" & feature_end <= cutoff
fixture_lines <- c(header, record_lines[is_chrI_declaration | in_slice])

fixture_fields <- strsplit(
    fixture_lines[!grepl("^#", fixture_lines)],
    "\t",
    fixed = TRUE
)
fixture_type <- vapply(fixture_fields, `[[`, character(1), 3L)
fixture_strand <- vapply(fixture_fields, `[[`, character(1), 7L)

stopifnot(
    sum(fixture_type == "gene") == 12L,
    sum(fixture_type == "mRNA") == 12L,
    sum(fixture_type == "CDS") == 12L,
    all(c("+", "-") %in% fixture_strand)
)

destination <- file.path("inst", "extdata", "example-annotation.gff3")
writeLines(fixture_lines, destination, useBytes = TRUE)

message(
    "Wrote ", destination, " with 12 complete genes from SGD R64-2-1 chrI."
)
