test_that("sample TSS tables combine without changing their inputs", {
    first <- data.table::data.table(
        chr = c("chrI", "chrI"),
        pos = c(100L, 300L),
        strand = c("+", "-"),
        sample_1 = c(2L, 3L)
    )
    second <- data.table::data.table(
        chr = c("chrI", "chrI"),
        pos = c(100L, 200L),
        strand = c("+", "+"),
        sample_2 = c(5L, 7L)
    )
    first_before <- as.data.frame(data.table::copy(first))
    second_before <- as.data.frame(data.table::copy(second))

    observed <- TSSr:::.combineSampleTSSTables(
        list(first, second),
        c("sample_1", "sample_2")
    )

    expect_identical(
        as.data.frame(observed),
        data.frame(
            chr = rep("chrI", 3L),
            pos = c(100L, 200L, 300L),
            strand = c("+", "+", "-"),
            sample_1 = c(2, 0, 3),
            sample_2 = c(5, 7, 0)
        )
    )
    expect_identical(as.data.frame(first), first_before)
    expect_identical(as.data.frame(second), second_before)
})

test_that("sample TSS tables preserve count types for empty samples", {
    integer_counts <- data.table::data.table(
        chr = "chrI", pos = 100L, strand = "+", integer_sample = 2L
    )
    empty_double_counts <- data.table::data.table(
        chr = character(), pos = integer(), strand = character(),
        double_sample = numeric()
    )

    observed <- TSSr:::.combineSampleTSSTables(
        list(integer_counts, empty_double_counts),
        c("integer_sample", "double_sample")
    )

    expect_identical(observed$integer_sample, 2L)
    expect_identical(observed$double_sample, 0)
    expect_identical(typeof(observed$pos), "integer")
})

test_that("sample TSS tables reject duplicate coordinates explicitly", {
    duplicated <- data.table::data.table(
        chr = c("chrI", "chrI"),
        pos = c(100L, 100L),
        strand = c("+", "+"),
        sample = c(2L, 3L)
    )

    expect_error(
        TSSr:::.combineSampleTSSTables(list(duplicated), "sample"),
        "Duplicate chr/pos/strand coordinates.*sample"
    )
})

test_that("sample TSS tables replace missing counts with typed zeros", {
    missing_integer <- data.table::data.table(
        chr = "chrI", pos = 100L, strand = "+", sample_1 = NA_integer_
    )
    missing_double <- data.table::data.table(
        chr = "chrI", pos = 200L, strand = "-", sample_2 = NA_real_
    )

    observed <- TSSr:::.combineSampleTSSTables(
        list(missing_integer, missing_double),
        c("sample_1", "sample_2")
    )

    expect_identical(observed$sample_1, c(0, 0))
    expect_identical(observed$sample_2, c(0, 0))
})

test_that("sample TSS tables retain every sample when all inputs are empty", {
    first <- data.table::data.table(
        chr = character(), pos = integer(), strand = character(),
        sample_1 = integer()
    )
    second <- data.table::data.table(
        chr = character(), pos = integer(), strand = character(),
        sample_2 = integer()
    )

    observed <- TSSr:::.combineSampleTSSTables(
        list(first, second),
        c("sample_1", "sample_2")
    )

    expect_identical(
        names(observed),
        c("chr", "pos", "strand", "sample_1", "sample_2")
    )
    expect_identical(nrow(observed), 0L)
    expect_identical(
        vapply(observed, typeof, character(1)),
        c(
            chr = "character", pos = "integer", strand = "character",
            sample_1 = "integer", sample_2 = "integer"
        )
    )
})
