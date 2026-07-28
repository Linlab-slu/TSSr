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

exampleTSSr@TSSrawMatrix <- exampleTSSr@TSSrawMatrix[
    exampleTSSr@TSSrawMatrix$chr == "chrI"
]
exampleTSSr@refTable <- exampleTSSr@refTable[
    as.character(exampleTSSr@refTable$seqnames) == "chrI"
]

sample_columns <- exampleTSSr@sampleLabels
exampleTSSr@librarySizes <- colSums(
    as.data.frame(exampleTSSr@TSSrawMatrix)[, sample_columns, drop = FALSE],
    na.rm = TRUE
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
exampleTSSr <- filterTSS(exampleTSSr, method = "TPM", tpmLow = 0.1)
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
