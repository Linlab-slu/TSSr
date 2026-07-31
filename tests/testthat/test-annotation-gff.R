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

test_that("GFF annotation matches the equivalent reference-table path", {
    gff_object <- make_gff_annotation_object()
    before <- tssr_content(gff_object)

    table_object <- gff_object
    table_object@refSource <- character()
    table_object@refTable <- example_annotation_reference()

    expect_warning(
        from_gff <- annotateCluster(
            gff_object,
            clusters = "consensusClusters",
            filterCluster = TRUE
        ),
        "genome version information is not available"
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
        expect_warning(
            result <- annotateCluster(
                object,
                clusters = "consensusClusters",
                filterCluster = FALSE
            ),
            "genome version information is not available"
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
