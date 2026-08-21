make_gff_annotation_object <- function() {
    data(exampleTSSr)
    object <- exampleTSSr
    object@consensusClusters <- lapply(
        object@consensusClusters,
        function(clusters) {
            clusters <- data.table::copy(clusters)
            clusters[dominant_tss <= 29000L]
        }
    )
    object@refSource <- system.file(
        "extdata", "example-annotation.gff3",
        package = "TSSr",
        mustWork = TRUE
    )
    object@refTable <- data.table::data.table()
    object@assignedClusters <- list()
    object@unassignedClusters <- list()
    object@filteredClusters <- list()
    object
}

example_annotation_reference <- function() {
    data.table::data.table(
        seqnames = rep("chrI", 12L),
        start = c(
            335L, 538L, 1807L, 2480L, 7235L, 10091L,
            11565L, 12046L, 13363L, 21566L, 22395L, 24000L
        ),
        end = c(
            649L, 792L, 2169L, 2707L, 9016L, 10399L,
            11951L, 12426L, 13743L, 21850L, 22685L, 27968L
        ),
        width = c(
            315L, 255L, 363L, 228L, 1782L, 309L,
            387L, 381L, 381L, 285L, 291L, 3969L
        ),
        strand = c(
            "+", "+", "-", "+", "-", "+",
            "-", "+", "-", "+", "-", "-"
        ),
        gene_id = c(
            "YAL069W", "YAL068W-A", "YAL068C", "YAL067W-A",
            "YAL067C", "YAL066W", "YAL065C", "YAL064W-B",
            "YAL064C-A", "YAL064W", "YAL063C-A", "YAL063C"
        )
    )
}

annotation_result_tables <- function(object, slot_name) {
    lapply(methods::slot(object, slot_name), function(value) {
        data.frame(value, check.names = FALSE)
    })
}

minimal_consensus_clusters <- function(position, strand = "+") {
    n <- length(position)
    data.table::data.table(
        cluster = seq_len(n),
        chr = rep("chrI", n),
        start = as.integer(position),
        end = as.integer(position),
        strand = rep(strand, length.out = n),
        dominant_tss = as.integer(position),
        tags = rep(10, n),
        tags.dominant_tss = rep(10, n),
        q_0.1 = as.integer(position),
        q_0.9 = as.integer(position),
        interquantile_width = rep(1, n)
    )
}

minimal_annotation_object <- function(clusters, reference) {
    object <- methods::new(
        "TSSr",
        genomeName = "BSgenome.Scerevisiae.UCSC.sacCer3",
        sampleLabels = "sample",
        sampleLabelsMerged = "sample",
        mergeIndex = 1,
        consensusClusters = list(sample = clusters)
    )
    object@refTable <- reference
    object
}

test_that("GFF annotation matches the equivalent reference-table path", {
    gff_object <- make_gff_annotation_object()
    before <- tssr_content(gff_object)

    table_object <- gff_object
    table_object@refSource <- character()
    table_object@refTable <- example_annotation_reference()

    expect_no_warning(
        from_gff <- annotateCluster(
            gff_object,
            clusters = "consensusClusters",
            filterCluster = TRUE
        )
    )
    from_table <- annotateCluster(
        table_object,
        clusters = "consensusClusters",
        filterCluster = TRUE
    )

    expect_tssr_content_equal(gff_object, before)
    expect_equal(
        annotation_result_tables(from_gff, "assignedClusters"),
        annotation_result_tables(from_table, "assignedClusters"),
        tolerance = 0
    )
    expect_equal(
        annotation_result_tables(from_gff, "unassignedClusters"),
        annotation_result_tables(from_table, "unassignedClusters"),
        tolerance = 0
    )
    expect_equal(
        annotation_result_tables(from_gff, "filteredClusters"),
        annotation_result_tables(from_table, "filteredClusters"),
        tolerance = 0
    )

    expect_setequal(
        from_gff@refTable$gene_id,
        example_annotation_reference()$gene_id
    )
    expect_identical(nrow(from_gff@refTable), 12L)
    expect_false("subset" %in% names(from_gff@refTable))
    expect_null(data.table::key(from_gff@refTable))

    for (sample_label in from_gff@sampleLabelsMerged) {
        assigned <- from_gff@assignedClusters[[sample_label]]
        unassigned <- from_gff@unassignedClusters[[sample_label]]

        expect_gt(nrow(assigned), 0L)
        expect_gt(nrow(unassigned), 0L)
        expect_true(all(!is.na(assigned$gene)))
        expect_true(all(is.na(unassigned$gene)))
        expect_false("subset" %in% names(assigned))
        expect_false("subset" %in% names(unassigned))
    }
})

test_that("organismName is optional descriptive metadata for GFF annotation", {
    annotate_with_name <- function(organism_name) {
        object <- make_gff_annotation_object()
        object@organismName <- organism_name
        expect_no_warning(
            result <- annotateCluster(
                object,
                clusters = "consensusClusters",
                filterCluster = FALSE
            )
        )
        result
    }

    without_name <- annotate_with_name(character())
    lowercase_name <- annotate_with_name("saccharomyces cerevisiae")
    canonical_name <- annotate_with_name("Saccharomyces cerevisiae")

    for (slot_name in c("refTable", "assignedClusters", "unassignedClusters")) {
        expect_equal(
            annotation_result_tables(without_name, slot_name),
            annotation_result_tables(lowercase_name, slot_name),
            tolerance = 0
        )
        expect_equal(
            annotation_result_tables(without_name, slot_name),
            annotation_result_tables(canonical_name, slot_name),
            tolerance = 0
        )
    }
    expect_length(without_name@organismName, 0L)
    expect_identical(
        lowercase_name@organismName,
        "saccharomyces cerevisiae"
    )
    expect_identical(
        canonical_name@organismName,
        "Saccharomyces cerevisiae"
    )
})

test_that("GFF transcript annotation uses transcript names as identifiers", {
    object <- make_gff_annotation_object()

    expect_no_warning(
        result <- annotateCluster(
            object,
            annotationType = "transcripts",
            filterCluster = FALSE
        )
    )

    expect_true("gene_id" %in% names(result@refTable))
    expect_true(all(grepl("_mRNA$", result@refTable$gene_id)))
    expect_true(all(vapply(
        result@assignedClusters,
        nrow,
        integer(1L)
    ) > 0L))
})

test_that("annotation preserves clusters with no promoter overlap", {
    reference <- data.table::data.table(
        seqnames = "chrI",
        start = 1000L,
        end = 2000L,
        width = 1001L,
        strand = "+",
        gene_id = "gene_1"
    )
    object <- minimal_annotation_object(
        minimal_consensus_clusters(5000L),
        reference
    )

    result <- annotateCluster(object, filterCluster = TRUE)

    expect_equal(nrow(result@unassignedClusters$sample), 1L)
    expect_true(is.na(result@unassignedClusters$sample$gene[[1L]]))
    expect_equal(nrow(result@filteredClusters$sample), 1L)
    expect_named(
        result@filteredClusters$sample,
        c(
            "cluster", "chr", "start", "end", "strand", "dominant_tss",
            "tags", "tags.dominant_tss", "q_0.1", "q_0.9",
            "interquantile_width", "gene"
        )
    )
})

test_that("annotation reports a missing consensus-cluster stage clearly", {
    object <- methods::new(
        "TSSr",
        genomeName = "BSgenome.DoesNotExist",
        sampleLabels = "sample",
        sampleLabelsMerged = "sample"
    )

    expect_error(
        annotateCluster(object),
        regexp = paste0(
            "consensusClusters.*empty.*run clusterTSS\\(\\) and ",
            "consensusCluster\\(\\)"
        )
    )
})

test_that("annotation rejects incompatible chromosome names", {
    reference <- data.table::data.table(
        seqnames = "chrII",
        start = 1000L,
        end = 2000L,
        width = 1001L,
        strand = "+",
        gene_id = "gene_1"
    )
    object <- minimal_annotation_object(
        minimal_consensus_clusters(500L),
        reference
    )

    expect_error(
        annotateCluster(object, filterCluster = FALSE),
        regexp = "No chromosome names are shared.*clusters.*reference"
    )
})

test_that("filtering handles clusters with no coding overlap", {
    reference <- data.table::data.table(
        seqnames = "chrI",
        start = 1000L,
        end = 2000L,
        width = 1001L,
        strand = "+",
        gene_id = "gene_1"
    )
    object <- minimal_annotation_object(
        minimal_consensus_clusters(500L),
        reference
    )

    expect_no_error(result <- annotateCluster(object, filterCluster = TRUE))
    expect_equal(nrow(result@assignedClusters$sample), 1L)
    expect_equal(nrow(result@filteredClusters$sample), 1L)
    expect_named(
        result@filteredClusters$sample,
        c(
            "cluster", "chr", "start", "end", "strand", "dominant_tss",
            "tags", "tags.dominant_tss", "q_0.1", "q_0.9",
            "interquantile_width", "gene"
        )
    )
})

test_that("negative-strand coding overlaps use coding-reference gene names", {
    reference <- data.table::data.table(
        seqnames = c("chrI", "chrI"),
        start = c(100L, 500L),
        end = c(1000L, 600L),
        width = c(901L, 101L),
        strand = c("-", "-"),
        gene_id = c("outer_gene", "inner_gene")
    )
    object <- minimal_annotation_object(
        minimal_consensus_clusters(200L, strand = "-"),
        reference
    )

    result <- annotateCluster(object, filterCluster = TRUE)

    expect_identical(
        result@unassignedClusters$sample$inCoding,
        "outer_gene"
    )
})

test_that("a TSStable import reaches GFF annotation through the full workflow", {
    input_file <- system.file(
        "extdata", "example-tss-table.tsv",
        package = "TSSr",
        mustWork = TRUE
    )
    annotation_file <- system.file(
        "extdata", "example-annotation.gff3",
        package = "TSSr",
        mustWork = TRUE
    )
    object <- TSSr(
        genomeName = "BSgenome.Scerevisiae.UCSC.sacCer3",
        inputFiles = input_file,
        inputFilesType = "TSStable",
        sampleLabels = c("SL01", "SL02", "SL03", "SL04"),
        sampleLabelsMerged = c("control", "treat"),
        mergeIndex = c(1, 1, 2, 2),
        refSource = annotation_file
    )
    expect_equal(nrow(object@refTable), 0L)

    object <- getTSS(object)
    object <- mergeSamples(object)
    object <- normalizeTSS(object)
    object <- filterTSS(object, method = "TPM", tpmLow = 2)
    object <- clusterTSS(
        object,
        method = "peakclu",
        clusterThreshold = 1,
        useMultiCore = FALSE
    )
    object <- consensusCluster(object, useMultiCore = FALSE)
    expect_no_warning(
        object <- annotateCluster(object, filterCluster = TRUE)
    )

    expect_true(all(vapply(
        object@assignedClusters,
        nrow,
        integer(1L)
    ) > 0L))
    expect_true(all(vapply(
        object@unassignedClusters,
        nrow,
        integer(1L)
    ) > 0L))
    expect_equal(nrow(object@refTable), 12L)
})

# The full-workflow test above is an integration smoke test. The deliberate
# zero-overlap and negative-strand cases remain the acceptance evidence for the
# annotation bugs because this fixture contains many ordinary overlaps.
