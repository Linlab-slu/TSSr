#!/usr/bin/env Rscript

# Rebuild the TSSr example data using chromosome I only.
#
# Run this script from the package root after installing TSSr. The installed
# exampleTSSr object is used as the source so that rebuilding never depends on
# partially modified objects in the source tree.

script_argument <- grep("^--file=", commandArgs(), value = TRUE)
if (length(script_argument) == 1L) {
    script_path <- sub("^--file=", "", script_argument)
    if (!grepl("^/", script_path)) {
        script_path <- file.path(Sys.getenv("PWD"), script_path)
    }
    script_path <- normalizePath(script_path)
    setwd(dirname(dirname(dirname(script_path))))
}

suppressPackageStartupMessages(library(TSSr))

data("exampleTSSr", package = "TSSr", envir = environment())

exampleTSSr@inputFiles <- character()
exampleTSSr@refSource <- character()
exampleTSSr@TSSrawMatrix <- exampleTSSr@TSSrawMatrix[
    exampleTSSr@TSSrawMatrix$chr == "chrI"
]

# Rebuild the reference table from the frozen bundled GFF instead of inheriting
# annotation state from the installed example object.
annotation_file <- "inst/extdata/example-annotation.gff3"
annotation_lines <- readLines(annotation_file, warn = FALSE)
annotation_lines <- annotation_lines[!startsWith(annotation_lines, "#")]
annotation <- data.table::fread(
    text = paste(annotation_lines, collapse = "\n"),
    sep = "\t",
    header = FALSE,
    col.names = c(
        "seqnames", "source", "type", "start", "end", "score",
        "strand", "phase", "attributes"
    )
)
reference_table <- annotation[
    type == "gene",
    .(
        seqnames,
        start,
        end,
        width = end - start + 1L,
        strand,
        gene_id = sub(";.*", "", sub("^ID=", "", attributes))
    )
]
data.table::setorder(reference_table, seqnames, start, end, strand, gene_id)
data.table::setkey(reference_table, NULL)
exampleTSSr@refTable <- reference_table

# Build a promoter-aware TSStable fixture. Six well-covered promoter
# neighborhoods provide realistic assigned clusters, while an equal number of
# background positions keeps the unassigned path represented.
fixture_columns <- c(
    "chr", "pos", "strand", exampleTSSr@sampleLabels
)
fixture_source <- data.table::copy(
    exampleTSSr@TSSrawMatrix[, .SD, .SDcols = fixture_columns]
)
data.table::setorder(fixture_source, strand, chr, pos)

fixture_promoter_genes <- c(
    "YAL067C", "YAL066W", "YAL064W-B",
    "YAL064W", "YAL063C-A", "YAL063C"
)
promoters <- data.table::copy(
    reference_table[gene_id %in% fixture_promoter_genes]
)
if (!setequal(promoters$gene_id, fixture_promoter_genes)) {
    stop("The bundled GFF does not contain every fixture promoter gene.")
}
promoters[, promoter_start := ifelse(
    strand == "+",
    pmax(1L, start - 1000L),
    end
)]
promoters[, promoter_end := ifelse(
    strand == "+",
    start,
    end + 1000L
)]

promoter_indices_by_gene <- lapply(seq_len(nrow(promoters)), function(i) {
    which(
        fixture_source$chr == promoters$seqnames[i] &
            fixture_source$strand == promoters$strand[i] &
            fixture_source$pos >= promoters$promoter_start[i] &
            fixture_source$pos <= promoters$promoter_end[i]
    )
})
names(promoter_indices_by_gene) <- promoters$gene_id
if (any(lengths(promoter_indices_by_gene) == 0L)) {
    stop("Every selected promoter must contain at least one TSS position.")
}

group_coverage <- vapply(promoter_indices_by_gene, function(indices) {
    rows <- fixture_source[indices]
    c(
        control = sum(rows$SL01 + rows$SL02),
        treat = sum(rows$SL03 + rows$SL04)
    )
}, numeric(2L))
if (any(group_coverage <= 0)) {
    stop("Every selected promoter must have nonzero control and treat coverage.")
}

promoter_indices <- sort(unique(unlist(promoter_indices_by_gene)))
promoter_rows <- fixture_source[promoter_indices]
background_source <- fixture_source[-promoter_indices]
background_source <- background_source[
    (SL01 + SL02) > 0 & (SL03 + SL04) > 0
]

select_evenly <- function(rows, n) {
    if (nrow(rows) < n) {
        stop("Not enough background TSS positions for the fixture.")
    }
    rows[as.integer(round(seq(1, nrow(rows), length.out = n)))]
}

promoter_strand_counts <- promoter_rows[, .N, by = strand]
background_rows <- data.table::rbindlist(lapply(
    seq_len(nrow(promoter_strand_counts)),
    function(i) {
        select_evenly(
            background_source[strand == promoter_strand_counts$strand[i]],
            promoter_strand_counts$N[i]
        )
    }
))
example_tss_table <- data.table::rbindlist(
    list(promoter_rows, background_rows),
    use.names = TRUE
)
data.table::setorder(example_tss_table, strand, chr, pos)
stopifnot(
    nrow(promoter_rows) == nrow(background_rows),
    nrow(example_tss_table) == 2L * nrow(promoter_rows),
    identical(unique(example_tss_table$strand), c("+", "-"))
)
data.table::fwrite(
    example_tss_table,
    file = "inst/extdata/example-tss-table.tsv",
    sep = "\t",
    quote = FALSE
)

exampleTSSr@tagClusters <- list()
exampleTSSr@consensusClusters <- list()
exampleTSSr@clusterShape <- list()
exampleTSSr@assignedClusters <- list()
exampleTSSr@unassignedClusters <- list()
exampleTSSr@filteredClusters <- list()
exampleTSSr@enhancers <- list()
exampleTSSr@DEtables <- list()
exampleTSSr@TAGtables <- list()
exampleTSSr@PromoterShift <- list()

exampleTSSr <- mergeSamples(exampleTSSr)
exampleTSSr <- normalizeTSS(exampleTSSr)
exampleTSSr <- filterTSS(exampleTSSr, method = "TPM", tpmLow = 2)
exampleTSSr <- clusterTSS(
    exampleTSSr,
    method = "peakclu",
    clusterThreshold = 1,
    useMultiCore = FALSE
)
exampleTSSr <- consensusCluster(exampleTSSr, useMultiCore = FALSE)
exampleTSSr <- annotateCluster(
    exampleTSSr,
    clusters = "consensusClusters",
    filterCluster = TRUE
)
exampleTSSr <- shapeCluster(
    exampleTSSr,
    clusters = "consensusClusters",
    method = "PSS",
    useMultiCore = FALSE
)
exampleTSSr <- deGene(
    exampleTSSr,
    comparePairs = list(c("control", "treat")),
    pval = 0.01,
    useMultiCore = FALSE
)
exampleTSSr <- callEnhancer(
    exampleTSSr,
    flanking = 400,
    dis2gene = 2000
)
exampleTSSr <- shiftPromoter(
    exampleTSSr,
    comparePairs = list(c("control", "treat")),
    pval = 0.01
)

save(exampleTSSr, file = "data/exampleTSSr.rda", compress = "xz")
