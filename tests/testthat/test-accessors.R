test_that("TSSmatrix returns an independent base data.frame", {
    data("exampleTSSr")

    raw <- TSSmatrix(exampleTSSr, data = "raw")

    expect_s3_class(raw, "data.frame")
    expect_false(data.table::is.data.table(raw))
    expect_equal(raw, as.data.frame(slot(exampleTSSr, "TSSrawMatrix")))

    original_chr <- slot(exampleTSSr, "TSSrawMatrix")$chr[[1]]
    raw$chr[[1]] <- "changed"
    expect_identical(
        slot(exampleTSSr, "TSSrawMatrix")$chr[[1]],
        original_chr
    )

    expect_error(
        TSSmatrix(exampleTSSr, data = "unknown"),
        "arg.*raw.*processed"
    )
})

test_that("tagClusters returns one or all results without exposing data.table", {
    data("exampleTSSr")

    control <- tagClusters(exampleTSSr, sample = "control")
    all_samples <- tagClusters(exampleTSSr)

    expect_s3_class(control, "data.frame")
    expect_false(data.table::is.data.table(control))
    expect_equal(
        control,
        as.data.frame(slot(exampleTSSr, "tagClusters")[["control"]])
    )
    expect_identical(names(all_samples), c("control", "treat"))
    expect_true(all(vapply(all_samples, is.data.frame, logical(1))))
    expect_false(any(vapply(
        all_samples,
        data.table::is.data.table,
        logical(1)
    )))

    original_chr <- slot(exampleTSSr, "tagClusters")[["control"]]$chr[[1]]
    control$chr[[1]] <- "changed"
    all_samples$control$chr[[1]] <- "also changed"
    expect_identical(
        slot(exampleTSSr, "tagClusters")[["control"]]$chr[[1]],
        original_chr
    )

    expect_error(
        tagClusters(exampleTSSr, sample = "missing"),
        "Unknown sample.*control.*treat"
    )
})

test_that("all result accessors return independent base data.frames", {
    data("exampleTSSr")
    accessors <- c(
        "tagClusters", "consensusClusters", "clusterShape",
        "assignedClusters", "unassignedClusters", "filteredClusters",
        "enhancers", "DEtables", "TAGtables", "PromoterShift"
    )

    for (accessor_name in accessors) {
        accessor <- get(accessor_name, mode = "function")
        stored <- slot(exampleTSSr, accessor_name)
        selected_name <- names(stored)[[1]]
        argument_name <- if (accessor_name %in%
            c("DEtables", "PromoterShift")) {
            "comparison"
        } else {
            "sample"
        }
        arguments <- setNames(list(selected_name), argument_name)
        expected <- stored[[selected_name]]
        if (identical(accessor_name, "DEtables")) {
            arguments$result <- "all"
            expected <- expected[["DEtable"]]
        }
        selected <- do.call(accessor, c(list(exampleTSSr), arguments))
        all_results <- if (identical(accessor_name, "DEtables")) {
            accessor(exampleTSSr, result = "all")
        } else {
            accessor(exampleTSSr)
        }

        expect_s3_class(selected, "data.frame")
        expect_false(
            data.table::is.data.table(selected),
            info = accessor_name
        )
        expect_equal(
            selected,
            as.data.frame(expected),
            info = accessor_name
        )
        expect_identical(
            names(all_results),
            names(stored),
            info = accessor_name
        )
        expect_true(
            all(vapply(all_results, is.data.frame, logical(1))),
            info = accessor_name
        )
        expect_false(
            any(vapply(all_results, data.table::is.data.table, logical(1))),
            info = accessor_name
        )

        original_names <- names(stored[[selected_name]])
        names(selected)[[1]] <- "changed"
        names(all_results[[selected_name]])[[1]] <- "also_changed"
        expect_identical(
            names(slot(exampleTSSr, accessor_name)[[selected_name]]),
            original_names,
            info = accessor_name
        )
    }
})

test_that("metadata accessors return independent user-facing values", {
    data("exampleTSSr")

    reference <- refTable(exampleTSSr)
    sizes <- librarySizes(exampleTSSr)

    expect_s3_class(reference, "data.frame")
    expect_false(data.table::is.data.table(reference))
    expect_equal(reference, as.data.frame(slot(exampleTSSr, "refTable")))
    expect_identical(sizes, slot(exampleTSSr, "librarySizes"))

    original_names <- names(slot(exampleTSSr, "refTable"))
    names(reference)[[1]] <- "changed"
    sizes[[1]] <- -1
    expect_identical(
        names(slot(exampleTSSr, "refTable")),
        original_names
    )
    expect_false(any(slot(exampleTSSr, "librarySizes") < 0))
})

test_that("DEtables selects all or significant results explicitly", {
    data("exampleTSSr")

    all_genes <- DEtables(
        exampleTSSr,
        comparison = "control_VS_treat",
        result = "all"
    )
    significant <- DEtables(
        exampleTSSr,
        comparison = "control_VS_treat",
        result = "significant"
    )

    expect_equal(
        all_genes,
        as.data.frame(
            slot(exampleTSSr, "DEtables")$control_VS_treat$DEtable
        )
    )
    expect_equal(
        significant,
        as.data.frame(
            slot(exampleTSSr, "DEtables")$control_VS_treat$DEsig
        )
    )
    expect_lte(nrow(significant), nrow(all_genes))
    expect_error(
        DEtables(exampleTSSr, comparison = "missing"),
        "Unknown comparison.*control_VS_treat"
    )
    expect_error(
        DEtables(exampleTSSr, result = "unknown"),
        "arg.*all.*significant"
    )
})

test_that("accessors handle an empty TSSr object", {
    object <- methods::new("TSSr")

    expect_s3_class(TSSmatrix(object, data = "raw"), "data.frame")
    expect_s3_class(TSSmatrix(object, data = "processed"), "data.frame")
    expect_s3_class(refTable(object), "data.frame")
    expect_identical(librarySizes(object), numeric())
    expect_identical(tagClusters(object), list())
    expect_identical(DEtables(object), list())
    expect_error(
        tagClusters(object, sample = "missing"),
        "Available values: none"
    )
})
