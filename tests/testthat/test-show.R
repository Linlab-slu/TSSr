format_expected_count <- function(value) {
    format(value, big.mark = ",", scientific = FALSE, trim = TRUE)
}

format_expected_table_summary <- function(tables) {
    paste(
        paste(
            names(tables),
            vapply(tables, function(table) {
                format_expected_count(nrow(table))
            }, character(1))
        ),
        collapse = "; "
    )
}

test_that("show keeps a populated TSSr object within one screen", {
    data(exampleTSSr)

    output <- capture.output(show(exampleTSSr))

    expect_lte(length(output), 20L)
})

test_that("show reports the core identity and data volume", {
    data(exampleTSSr)

    output <- capture.output(show(exampleTSSr))

    expect_true(any(output == sprintf(
        "  Genome: %s", exampleTSSr@genomeName
    )))
    expect_true(any(output == sprintf(
        "  Samples: %d (%s)",
        length(exampleTSSr@sampleLabels),
        paste(exampleTSSr@sampleLabels, collapse = ", ")
    )))
    expect_true(any(output == sprintf(
        "  Merged samples: %d (%s)",
        length(exampleTSSr@sampleLabelsMerged),
        paste(exampleTSSr@sampleLabelsMerged, collapse = ", ")
    )))
    expect_true(any(output == sprintf(
        "  TSSs: raw %s; processed %s",
        format_expected_count(nrow(exampleTSSr@TSSrawMatrix)),
        format_expected_count(nrow(exampleTSSr@TSSprocessedMatrix))
    )))
})

test_that("show summarizes completed analyses without expanding result tables", {
    data(exampleTSSr)

    output <- capture.output(show(exampleTSSr))
    expected <- c(
        "  Analyses:",
        paste0(
            "    Tag clusters: ",
            format_expected_table_summary(exampleTSSr@tagClusters)
        ),
        paste0(
            "    Consensus clusters: ",
            format_expected_table_summary(exampleTSSr@consensusClusters)
        ),
        paste0(
            "    Cluster shapes: ",
            format_expected_table_summary(exampleTSSr@clusterShape)
        ),
        paste0(
            "    Assigned clusters: ",
            format_expected_table_summary(exampleTSSr@assignedClusters)
        ),
        paste0(
            "    Unassigned clusters: ",
            format_expected_table_summary(exampleTSSr@unassignedClusters)
        ),
        paste0(
            "    Filtered clusters: ",
            format_expected_table_summary(exampleTSSr@filteredClusters)
        ),
        paste0(
            "    Enhancers: ",
            format_expected_table_summary(exampleTSSr@enhancers)
        ),
        sprintf(
            "    DE comparisons: %d (%s)",
            length(exampleTSSr@DEtables),
            paste(names(exampleTSSr@DEtables), collapse = ", ")
        ),
        paste0(
            "    TAG tables: ",
            format_expected_table_summary(exampleTSSr@TAGtables)
        ),
        paste0(
            "    Promoter shifts: ",
            format_expected_table_summary(exampleTSSr@PromoterShift)
        )
    )

    expect_true(all(expected %in% output))
})

test_that("show handles an empty TSSr object with explicit placeholders", {
    output <- capture.output(show(methods::new("TSSr")))

    expect_lte(length(output), 20L)
    expect_true(any(output == "  Genome: <unset>"))
    expect_true(any(output == "  Samples: 0 (<none>)"))
    expect_true(any(output == "  Merged samples: 0 (<none>)"))
    expect_true(any(output == "  TSSs: raw 0; processed 0"))
    expect_equal(sum(grepl("<not run>", output, fixed = TRUE)), 10L)
})

test_that("show truncates long sample and result label lists", {
    labels <- paste0("sample", seq_len(7L))
    groups <- paste0("group", seq_len(7L))
    tables <- stats::setNames(
        replicate(7L, data.frame(value = 1), simplify = FALSE),
        groups
    )
    object <- new(
        "TSSr",
        sampleLabels = labels,
        sampleLabelsMerged = groups,
        mergeIndex = seq_len(7L),
        tagClusters = tables
    )

    output <- capture.output(show(object))

    expect_true(any(output == paste0(
        "  Samples: 7 (sample1, sample2, sample3, sample4, ... +3)"
    )))
    expect_true(any(output == paste0(
        "    Tag clusters: group1 1; group2 1; group3 1; group4 1; ... +3"
    )))
})

test_that("displaying a TSSr object is read-only and print uses its show method", {
    data(exampleTSSr)
    before <- tssr_content(exampleTSSr)

    show_output <- capture.output(visibility <- withVisible(show(exampleTSSr)))
    print_output <- capture.output(print(exampleTSSr))

    expect_false(visibility$visible)
    expect_s4_class(visibility$value, "TSSr")
    expect_tssr_content_equal(exampleTSSr, before)
    expect_identical(print_output, show_output)
})

test_that("the TSSr show method is part of the documented package interface", {
    expect_true("show" %in% getNamespaceExports("TSSr"))
})
