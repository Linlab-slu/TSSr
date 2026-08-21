# Test the core TSSr workflow: mergeSamples -> normalizeTSS -> filterTSS ->
# clusterTSS -> consensusCluster
# Run workflow once, test each step's output

data(exampleTSSr)

make_small_workflow_object <- function() {
    object <- exampleTSSr
    processed <- object@TSSprocessedMatrix
    selected_rows <- unlist(lapply(c("+", "-"), function(strand_value) {
        head(which(processed$strand == strand_value), 500L)
    }))
    object@TSSprocessedMatrix <- data.table::copy(processed[selected_rows])
    object@tagClusters <- list()
    object@consensusClusters <- list()
    object
}

test_that("mergeSamples merges raw TSS data correctly", {
    result <- mergeSamples(exampleTSSr)

    processed <- result@TSSprocessedMatrix
    expect_true(nrow(processed) > 0)
    merged_labels <- result@sampleLabelsMerged
    expect_true(all(merged_labels %in% names(processed)))
})

test_that("normalizeTSS normalizes to TPM", {
    object <- mergeSamples(exampleTSSr)
    result <- normalizeTSS(object)

    processed <- result@TSSprocessedMatrix
    merged_labels <- result@sampleLabelsMerged
    first_col <- processed[[merged_labels[1]]]
    nonzero <- first_col[first_col > 0]
    ## TPM values should have decimals
    expect_true(any(nonzero != round(nonzero)))
})

test_that("filterTSS with TPM method reduces rows", {
    object <- normalizeTSS(mergeSamples(exampleTSSr))
    rows_before <- nrow(object@TSSprocessedMatrix)
    result <- filterTSS(object, method = "TPM", tpmLow = 2)
    rows_after <- nrow(result@TSSprocessedMatrix)

    expect_true(rows_after > 0)
    expect_lt(rows_after, rows_before)
})

test_that("clusterTSS produces tagClusters", {
    object <- make_small_workflow_object()
    before <- tssr_content(object)
    result <- clusterTSS(object,
        method = "peakclu", clusterThreshold = 1,
        useMultiCore = FALSE
    )

    expect_tssr_content_equal(object, before)
    expect_s4_class(result, "TSSr")
    tc <- result@tagClusters
    expect_type(tc, "list")
    expect_true(length(tc) > 0)

    first_tc <- tc[[1]]
    expect_true(is.data.frame(first_tc) || inherits(first_tc, "data.table"))
    expect_true(nrow(first_tc) > 0)
    expect_true(all(c("chr", "start", "end", "strand") %in% names(first_tc)))
})

test_that("consensusCluster produces consensus clusters", {
    object <- clusterTSS(
        make_small_workflow_object(),
        method = "peakclu",
        clusterThreshold = 1,
        useMultiCore = FALSE
    )
    before <- tssr_content(object)
    result <- consensusCluster(object, useMultiCore = FALSE)

    expect_tssr_content_equal(object, before)
    expect_s4_class(result, "TSSr")
    cc <- result@consensusClusters
    expect_type(cc, "list")
    expect_true(length(cc) > 0)

    first_cc <- cc[[1]]
    expect_true(is.data.frame(first_cc) || inherits(first_cc, "data.table"))
    expect_true(all(c("chr", "start", "end", "strand") %in% names(first_cc)))
})

test_that("workflow methods report missing upstream stages clearly", {
    empty <- methods::new(
        "TSSr",
        sampleLabels = c("control", "treat"),
        sampleLabelsMerged = c("control", "treat"),
        mergeIndex = c(1, 2)
    )
    clustered_input <- empty
    clustered_input@TSSprocessedMatrix <- data.table::data.table(
        chr = "chrI",
        pos = 1L,
        strand = "+",
        control = 0,
        treat = 0
    )

    expect_error(
        mergeSamples(empty),
        "TSSrawMatrix.*empty.*getTSS"
    )
    expect_error(
        normalizeTSS(empty),
        "TSSprocessedMatrix.*empty.*getTSS"
    )
    expect_error(
        clusterTSS(empty),
        "TSSprocessedMatrix.*empty.*getTSS"
    )
    expect_error(
        consensusCluster(clustered_input),
        "tagClusters.*empty.*clusterTSS"
    )
    expect_error(
        shapeCluster(clustered_input),
        "consensusClusters.*empty.*consensusCluster"
    )
    expect_error(
        deGene(empty),
        "assignedClusters.*empty.*annotateCluster"
    )
    expect_error(
        callEnhancer(empty),
        "assignedClusters.*empty.*annotateCluster"
    )
    expect_error(
        shiftPromoter(empty, comparePairs = list(c("control", "treat"))),
        "assignedClusters.*empty.*annotateCluster"
    )
})

test_that("clusterTSS explains when filtering removes every TSS", {
    input_file <- system.file(
        "extdata", "example-tss-table.tsv",
        package = "TSSr",
        mustWork = TRUE
    )
    object <- TSSr(
        genomeName = "BSgenome.Scerevisiae.UCSC.sacCer3",
        inputFiles = input_file,
        inputFilesType = "TSStable",
        sampleLabels = c("SL01", "SL02", "SL03", "SL04"),
        sampleLabelsMerged = c("control", "treat"),
        mergeIndex = c(1, 1, 2, 2)
    )
    object <- getTSS(object)
    object <- mergeSamples(object)
    object <- normalizeTSS(object)
    object <- filterTSS(object, method = "TPM", tpmLow = 1e12)

    expect_equal(nrow(TSSmatrix(object, data = "processed")), 0L)
    expect_error(
        clusterTSS(object, useMultiCore = FALSE),
        paste0(
            "TSSprocessedMatrix.*empty.*relax filterTSS\\(\\) ",
            "thresholds"
        )
    )
})
