# Test export functions — all write to tempdir to avoid side effects

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

test_that("exportClustersTable writes tag cluster tables", {
    data(exampleTSSr)
    tmpdir <- tempdir()
    withr::with_dir(tmpdir, {
        exportClustersTable(exampleTSSr, data = "tagClusters")
        samples <- exampleTSSr@sampleLabelsMerged
        for (s in samples) {
            fname <- paste0(s, ".tagClusters.txt")
            expect_true(file.exists(fname))
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
            expect_true(file.exists(fname))
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
            expect_true(file.exists(fname))
            unlink(fname)
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
            expect_true(file.exists(fname))
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
            expect_true(file.exists(fname))
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
            expect_true(file.exists(fname))
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
            expect_true(file.exists(fname))
            unlink(fname)
        }
    })
})
