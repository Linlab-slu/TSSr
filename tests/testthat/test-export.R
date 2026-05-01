# Test export functions — all write to tempdir to avoid side effects

test_that("exportTSStable writes TSS table to file", {
    data(exampleTSSr)
    tmpdir <- tempdir()
    withr::with_dir(tmpdir, {
        exportTSStable(exampleTSSr, data = "raw", merged = "TRUE")
        expect_true(file.exists("ALL.samples.TSS.raw.txt"))
        content <- read.table("ALL.samples.TSS.raw.txt", header = TRUE, sep = "\t")
        expect_true(nrow(content) > 0)
        unlink("ALL.samples.TSS.raw.txt")
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
    if (length(exampleTSSr@assignedClusters) > 0) {
        withr::with_dir(tmpdir, {
            exportClustersTable(exampleTSSr, data = "assigned")
            samples <- exampleTSSr@sampleLabelsMerged
            for (s in samples) {
                fname <- paste0(s, ".assignedClusters.txt")
                expect_true(file.exists(fname))
                unlink(fname)
            }
        })
    }
})

test_that("exportShapeTable writes shape tables", {
    data(exampleTSSr)
    tmpdir <- tempdir()
    if (length(exampleTSSr@clusterShape) > 0) {
        withr::with_dir(tmpdir, {
            exportShapeTable(exampleTSSr)
            samples <- exampleTSSr@sampleLabelsMerged
            for (s in samples) {
                fname <- paste0(s, ".promoter.shape.txt")
                expect_true(file.exists(fname))
                unlink(fname)
            }
        })
    }
})

test_that("exportEnhancerTable writes enhancer tables", {
    data(exampleTSSr)
    tmpdir <- tempdir()
    if (length(exampleTSSr@enhancers) > 0) {
        withr::with_dir(tmpdir, {
            exportEnhancerTable(exampleTSSr)
            samples <- exampleTSSr@sampleLabelsMerged
            for (s in samples) {
                fname <- paste0(s, ".enhancers.txt")
                expect_true(file.exists(fname))
                unlink(fname)
            }
        })
    }
})

test_that("exportDETable writes DE tables", {
    data(exampleTSSr)
    tmpdir <- tempdir()
    if (length(exampleTSSr@DEtables) > 0) {
        withr::with_dir(tmpdir, {
            exportDETable(exampleTSSr, data = "sig")
            d_names <- names(exampleTSSr@DEtables)
            for (d in d_names) {
                fname <- paste0(d, ".DE.table.sig.txt")
                expect_true(file.exists(fname))
                unlink(fname)
            }
        })
    }
})

test_that("exportShiftTable writes shift tables", {
    data(exampleTSSr)
    tmpdir <- tempdir()
    if (length(exampleTSSr@PromoterShift) > 0) {
        withr::with_dir(tmpdir, {
            exportShiftTable(exampleTSSr)
            d_names <- names(exampleTSSr@PromoterShift)
            for (d in d_names) {
                fname <- paste0(d, ".promoter.shift.table.txt")
                expect_true(file.exists(fname))
                unlink(fname)
            }
        })
    }
})
