test_that("show keeps a populated TSSr object within one screen", {
    data(exampleTSSr)

    output <- capture.output(show(exampleTSSr))

    expect_lte(length(output), 20L)
})

test_that("show reports the core identity and data volume", {
    data(exampleTSSr)

    output <- capture.output(show(exampleTSSr))

    expect_true(any(output == "  Genome: BSgenome.Scerevisiae.UCSC.sacCer3"))
    expect_true(any(output == "  Samples: 4 (SL01, SL02, SL03, SL04)"))
    expect_true(any(output == "  Merged samples: 2 (control, treat)"))
    expect_true(any(output == "  TSSs: raw 29,456; processed 29,456"))
})

test_that("show summarizes completed analyses without expanding result tables", {
    data(exampleTSSr)

    output <- capture.output(show(exampleTSSr))
    expected <- c(
        "  Analyses:",
        "    Tag clusters: control 1,066; treat 1,098",
        "    Consensus clusters: control 1,065; treat 1,097",
        "    Cluster shapes: control 1,065; treat 1,097",
        "    Assigned clusters: control 267; treat 264",
        "    Unassigned clusters: control 798; treat 833",
        "    Filtered clusters: control 781; treat 868",
        "    Enhancers: control 21; treat 14",
        "    DE comparisons: 1 (control_VS_treat)",
        "    TAG tables: control 267; treat 264",
        "    Promoter shifts: control_VS_treat 30"
    )

    expect_true(all(expected %in% output))
})

test_that("show handles an empty TSSr object with explicit placeholders", {
    output <- capture.output(show(TSSr()))

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
