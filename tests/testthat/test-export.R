# Test export functions — all write to tempdir to avoid side effects

expect_exported_table <- function(path, expected) {
    observed <- as.data.frame(
        data.table::fread(path),
        check.names = FALSE
    )
    expected <- as.data.frame(
        TSSr:::.prepareExportTable(expected),
        check.names = FALSE
    )

    expect_identical(names(observed), names(expected))
    expect_identical(nrow(observed), nrow(expected))
    for (column in names(expected)) {
        if (is.integer(expected[[column]])) {
            expect_identical(
                as.integer(observed[[column]]),
                expected[[column]],
                info = column
            )
        } else if (is.numeric(expected[[column]])) {
            expect_equal(
                as.numeric(observed[[column]]),
                expected[[column]],
                tolerance = 1e-12,
                info = column
            )
        } else {
            expect_identical(
                as.character(observed[[column]]),
                as.character(expected[[column]]),
                info = column
            )
        }
    }
}

test_that("exportTSStable keeps raw and processed data semantics distinct", {
    data(exampleTSSr)
    tmpdir <- tempdir()
    withr::with_dir(tmpdir, {
        exportTSStable(exampleTSSr, data = "raw", merged = TRUE)
        exportTSStable(exampleTSSr, data = "raw", merged = FALSE)
        exportTSStable(exampleTSSr, data = "processed")

        merged_file <- "ALL.samples.TSS.raw.merged.txt"
        raw_file <- "ALL.samples.TSS.raw.txt"
        processed_file <- "ALL.samples.TSS.processed.txt"
        expect_true(all(file.exists(c(
            merged_file, raw_file, processed_file
        ))))

        merged <- data.table::fread(merged_file)
        raw <- data.table::fread(raw_file)
        processed <- data.table::fread(processed_file)
        raw_matrix <- exampleTSSr@TSSrawMatrix
        expected_merged <- data.table::data.table(
            control = rowSums(raw_matrix[, c("SL01", "SL02")]),
            treat = rowSums(raw_matrix[, c("SL03", "SL04")])
        )

        expect_identical(names(raw), c(
            "chr", "pos", "strand", exampleTSSr@sampleLabels
        ))
        expect_identical(names(merged), c(
            "chr", "pos", "strand", exampleTSSr@sampleLabelsMerged
        ))
        expect_identical(names(processed), names(
            exampleTSSr@TSSprocessedMatrix
        ))
        expect_equal(merged[, c("control", "treat")], expected_merged)
        expect_equal(raw, exampleTSSr@TSSrawMatrix)
        expect_equal(processed, exampleTSSr@TSSprocessedMatrix)
    })
})

test_that("exportTSStable uses a custom prefix with the established suffix", {
    data(exampleTSSr)

    withr::with_tempdir({
        exportTSStable(
            exampleTSSr,
            data = "raw",
            merged = FALSE,
            outputPrefix = "project1"
        )
        exportTSStable(
            exampleTSSr,
            data = "raw",
            merged = TRUE,
            outputPrefix = "project2"
        )
        exportTSStable(
            exampleTSSr,
            data = "processed",
            outputPrefix = "project3"
        )

        expect_true(file.exists("project1.TSS.raw.txt"))
        expect_true(file.exists("project2.TSS.raw.merged.txt"))
        expect_true(file.exists("project3.TSS.processed.txt"))
        expect_false(file.exists("ALL.samples.TSS.raw.txt"))
        expect_exported_table(
            "project1.TSS.raw.txt",
            exampleTSSr@TSSrawMatrix
        )
        expect_exported_table(
            "project3.TSS.processed.txt",
            exampleTSSr@TSSprocessedMatrix
        )
        header <- readLines("project1.TSS.raw.txt", n = 5L)
        expect_true(any(grepl("^# dataType: raw$", header)))
    })
})

test_that("exportTSStable accepts an explicit dot-separated suffix", {
    data(exampleTSSr)

    withr::with_tempdir({
        exportTSStable(
            exampleTSSr,
            data = "raw",
            merged = FALSE,
            outputPrefix = "project1",
            outputSuffix = "TSS.raw.txt"
        )

        expect_true(file.exists("project1.TSS.raw.txt"))
        expect_false(file.exists("ALL.samples.TSS.raw.txt"))
    })
})

test_that("exportTSStable refuses to overwrite an existing file by default", {
    data(exampleTSSr)

    withr::with_tempdir({
        target <- "project1.TSS.raw.txt"
        writeLines("keep this file", target)

        expect_error(
            exportTSStable(
                exampleTSSr,
                data = "raw",
                merged = FALSE,
                outputPrefix = "project1"
            ),
            "already exists.*outputPrefix.*overwrite = TRUE"
        )
        expect_identical(readLines(target), "keep this file")

        expect_no_error(exportTSStable(
            exampleTSSr,
            data = "raw",
            merged = FALSE,
            outputPrefix = "project1",
            overwrite = TRUE
        ))
        expect_exported_table(target, exampleTSSr@TSSrawMatrix)
    })
})

test_that("exportTSStable rejects a suffix that contradicts the data type", {
    data(exampleTSSr)

    withr::with_tempdir({
        expect_error(
            exportTSStable(
                exampleTSSr,
                data = "raw",
                merged = FALSE,
                outputPrefix = "project1",
                outputSuffix = "TSS.processed.txt"
            ),
            "raw data.*processed.*outputSuffix"
        )
        expect_error(
            exportTSStable(
                exampleTSSr,
                data = "processed",
                outputPrefix = "project1",
                outputSuffix = "TSS.raw.txt"
            ),
            "processed data.*raw.*outputSuffix"
        )
    })
})

test_that("exportTSStable keeps directory paths out of outputSuffix", {
    data(exampleTSSr)

    withr::with_tempdir({
        expect_error(
            exportTSStable(
                exampleTSSr,
                data = "raw",
                merged = FALSE,
                outputPrefix = "project1",
                outputSuffix = "nested/TSS.raw.txt"
            ),
            "outputSuffix.*path separators"
        )
    })
})

test_that("exportTSStable requires one non-empty outputPrefix", {
    data(exampleTSSr)

    withr::with_tempdir({
        expect_error(
            exportTSStable(
                exampleTSSr,
                data = "raw",
                merged = FALSE,
                outputPrefix = ""
            ),
            "outputPrefix.*one non-empty string"
        )
        expect_error(
            exportTSStable(
                exampleTSSr,
                data = "raw",
                merged = FALSE,
                outputPrefix = c("project1", "project2")
            ),
            "outputPrefix.*one non-empty string"
        )
    })
})

test_that("exportTSStable requires a single logical overwrite value", {
    data(exampleTSSr)

    withr::with_tempdir({
        expect_error(
            exportTSStable(
                exampleTSSr,
                data = "raw",
                merged = FALSE,
                outputPrefix = "project1",
                overwrite = NA
            ),
            "overwrite must be TRUE or FALSE"
        )
    })
})

test_that("exportTSStable accepts an existing directory in outputPrefix", {
    data(exampleTSSr)

    withr::with_tempdir({
        dir.create("results")
        exportTSStable(
            exampleTSSr,
            data = "raw",
            merged = FALSE,
            outputPrefix = file.path("results", "project1")
        )
        expect_true(file.exists(file.path(
            "results", "project1.TSS.raw.txt"
        )))

        expect_error(
            exportTSStable(
                exampleTSSr,
                data = "raw",
                merged = FALSE,
                outputPrefix = file.path("missing", "project1")
            ),
            "Output directory does not exist.*missing"
        )
    })
})

test_that("exportClustersTable writes tag cluster tables", {
    data(exampleTSSr)
    tmpdir <- tempdir()
    withr::with_dir(tmpdir, {
        exportClustersTable(exampleTSSr, data = "tagClusters")
        samples <- exampleTSSr@sampleLabelsMerged
        for (s in samples) {
            fname <- paste0(s, ".tagClusters.txt")
            expect_exported_table(fname, tagClusters(exampleTSSr, s))
            unlink(fname)
        }
    })
})

test_that("exportClustersTable writes consensus cluster tables", {
    data(exampleTSSr)
    tmpdir <- tempdir()
    withr::with_dir(tmpdir, {
        exportClustersTable(exampleTSSr, data = "consensusClusters")
        samples <- exampleTSSr@sampleLabelsMerged
        for (s in samples) {
            fname <- paste0(s, ".consensusClusters.txt")
            expect_exported_table(
                fname,
                consensusClusters(exampleTSSr, s)
            )
            unlink(fname)
        }
    })
})

test_that("exportClustersTable writes assigned cluster tables", {
    data(exampleTSSr)
    tmpdir <- tempdir()
    expect_gt(length(exampleTSSr@assignedClusters), 0L)
    withr::with_dir(tmpdir, {
        exportClustersTable(exampleTSSr, data = "assigned")
        samples <- exampleTSSr@sampleLabelsMerged
        for (s in samples) {
            fname <- paste0(s, ".assignedClusters.txt")
            expect_exported_table(fname, assignedClusters(exampleTSSr, s))
            unlink(fname)
        }
    })
})

test_that("exportClustersTable writes unassigned cluster table contents", {
    data(exampleTSSr)
    tmpdir <- tempdir()
    expect_gt(length(exampleTSSr@unassignedClusters), 0L)
    withr::with_dir(tmpdir, {
        exportClustersTable(exampleTSSr, data = "unassigned")
        for (sample in exampleTSSr@sampleLabelsMerged) {
            path <- paste0(sample, ".unassignedClusters.txt")
            expect_exported_table(
                path,
                unassignedClusters(exampleTSSr, sample)
            )
            unlink(path)
        }
    })
})

test_that("exportShapeTable writes shape tables", {
    data(exampleTSSr)
    tmpdir <- tempdir()
    expect_gt(length(exampleTSSr@clusterShape), 0L)
    withr::with_dir(tmpdir, {
        exportShapeTable(exampleTSSr)
        samples <- exampleTSSr@sampleLabelsMerged
        for (s in samples) {
            fname <- paste0(s, ".promoter.shape.txt")
            expect_exported_table(fname, clusterShape(exampleTSSr, s))
            unlink(fname)
        }
    })
})

test_that("exportEnhancerTable writes enhancer tables", {
    data(exampleTSSr)
    tmpdir <- tempdir()
    expect_gt(length(exampleTSSr@enhancers), 0L)
    withr::with_dir(tmpdir, {
        exportEnhancerTable(exampleTSSr)
        samples <- exampleTSSr@sampleLabelsMerged
        for (s in samples) {
            fname <- paste0(s, ".enhancers.txt")
            expect_exported_table(fname, enhancers(exampleTSSr, s))
            unlink(fname)
        }
    })
})

test_that("exportDETable writes DE tables", {
    data(exampleTSSr)
    tmpdir <- tempdir()
    expect_gt(length(exampleTSSr@DEtables), 0L)
    withr::with_dir(tmpdir, {
        exportDETable(exampleTSSr, data = "sig")
        d_names <- names(exampleTSSr@DEtables)
        for (d in d_names) {
            fname <- paste0(d, ".DE.table.sig.txt")
            expect_exported_table(
                fname,
                DEtables(exampleTSSr, d, result = "significant")
            )
            unlink(fname)
        }
    })
})

test_that("exportShiftTable writes shift tables", {
    data(exampleTSSr)
    tmpdir <- tempdir()
    expect_gt(length(exampleTSSr@PromoterShift), 0L)
    withr::with_dir(tmpdir, {
        exportShiftTable(exampleTSSr)
        d_names <- names(exampleTSSr@PromoterShift)
        for (d in d_names) {
            fname <- paste0(d, ".promoter.shift.table.txt")
            expect_exported_table(fname, PromoterShift(exampleTSSr, d))
            unlink(fname)
        }
    })
})
