tsstable_fixture <- function() {
    data.table::fread(system.file(
        "extdata", "example-tss-table.tsv",
        package = "TSSr"
    ))
}

canonical_tsstable <- function(x) {
    x <- data.table::as.data.table(x)
    data.table::setorderv(x, c("chr", "pos", "strand"))
    as.data.frame(x)
}

write_tsstable_with_metadata <- function(table, path, metadata) {
    metadata_lines <- paste0(
        "# ", names(metadata), ": ", unname(metadata)
    )
    table_path <- tempfile()
    on.exit(unlink(table_path))
    data.table::fwrite(table, table_path, sep = "\t")
    writeLines(c(metadata_lines, readLines(table_path)), path)
}

test_that("multiple single-sample TSS tables combine into the original table", {
    fixture <- tsstable_fixture()
    sample_names <- names(fixture)[-(1:3)]

    withr::with_tempdir({
        input_files <- vapply(sample_names, function(sample_name) {
            path <- paste0(sample_name, ".tsv")
            data.table::fwrite(
                fixture[, c("chr", "pos", "strand", sample_name),
                    with = FALSE
                ],
                path, sep = "\t"
            )
            path
        }, character(1))

        combined <- combineTSSTables(input_files)

        expect_s3_class(combined, "data.frame")
        expect_false(data.table::is.data.table(combined))
        expect_identical(
            canonical_tsstable(combined),
            canonical_tsstable(fixture)
        )
    })
})

test_that(
    "sparse single- and multi-sample tables combine by coordinate union",
    {
    fixture <- tsstable_fixture()
    sample_groups <- list("SL01", c("SL02", "SL03"), "SL04")

    withr::with_tempdir({
        input_files <- vapply(seq_along(sample_groups), function(i) {
            samples <- sample_groups[[i]]
            table <- fixture[
                rowSums(fixture[, ..samples]) > 0,
                c("chr", "pos", "strand", samples),
                with = FALSE
            ]
            path <- paste0("part", i, ".tsv")
            data.table::fwrite(table, path, sep = "\t")
            path
        }, character(1))

        combined <- combineTSSTables(input_files)

        expect_identical(
            canonical_tsstable(combined),
            canonical_tsstable(fixture)
        )
    })
    }
)

test_that("selected samples split into one table without all-zero rows", {
    fixture_path <- system.file(
        "extdata", "example-tss-table.tsv",
        package = "TSSr"
    )
    fixture <- tsstable_fixture()
    samples <- c("SL01", "SL02")
    expected <- fixture[
        rowSums(fixture[, ..samples]) > 0,
        c("chr", "pos", "strand", samples),
        with = FALSE
    ]

    result <- splitTSSTable(fixture_path, samples = samples)

    expect_s3_class(result, "data.frame")
    expect_false(data.table::is.data.table(result))
    expect_identical(result, as.data.frame(expected))
    expect_identical(
        nrow(splitTSSTable(
            fixture_path, samples = samples, dropZero = FALSE
        )),
        nrow(fixture)
    )
})

test_that("sample-wise split and combine round trip preserves the table", {
    fixture_path <- system.file(
        "extdata", "example-tss-table.tsv",
        package = "TSSr"
    )
    fixture <- tsstable_fixture()
    sample_names <- names(fixture)[-(1:3)]

    withr::with_tempdir({
        split_files <- vapply(sample_names, function(sample_name) {
            path <- paste0(sample_name, ".tsv")
            splitTSSTable(
                fixture_path,
                samples = sample_name,
                outputFile = path
            )
            path
        }, character(1))

        combined <- combineTSSTables(
            split_files,
            outputFile = "combined.tsv"
        )

        expect_true(file.exists("combined.tsv"))
        expect_identical(
            canonical_tsstable(combined),
            canonical_tsstable(fixture)
        )
        expect_identical(
            canonical_tsstable(data.table::fread("combined.tsv")),
            canonical_tsstable(fixture)
        )
    })
})

test_that("TSS table coordinate columns are validated before combining", {
    fixture <- tsstable_fixture()

    withr::with_tempdir({
        data.table::fwrite(
            fixture[, c("chr", "pos", "strand", "SL01"), with = FALSE],
            "valid.tsv", sep = "\t"
        )
        invalid <- fixture[, c("chr", "pos", "strand", "SL02"),
            with = FALSE
        ]
        data.table::setnames(invalid, "chr", "chromosome")
        data.table::fwrite(invalid, "invalid.tsv", sep = "\t")

        expect_error(
            combineTSSTables(c("valid.tsv", "invalid.tsv")),
            "invalid.tsv.*first three columns.*chr.*pos.*strand"
        )
    })
})

test_that("non-integer sample values are rejected as non-raw counts", {
    fixture <- tsstable_fixture()

    withr::with_tempdir({
        data.table::fwrite(
            fixture[, c("chr", "pos", "strand", "SL01"), with = FALSE],
            "raw.tsv", sep = "\t"
        )
        normalized <- fixture[, c("chr", "pos", "strand", "SL02"),
            with = FALSE
        ]
        normalized[, SL02 := as.numeric(SL02)]
        normalized[1L, SL02 := 0.5]
        data.table::fwrite(normalized, "normalized.tsv", sep = "\t")

        expect_error(
            combineTSSTables(c("raw.tsv", "normalized.tsv")),
            "normalized.tsv.*non-negative integer raw counts"
        )
    })
})

test_that("duplicate genomic coordinates are rejected", {
    fixture <- tsstable_fixture()

    withr::with_tempdir({
        first <- fixture[, c("chr", "pos", "strand", "SL01"),
            with = FALSE
        ]
        duplicated <- fixture[, c("chr", "pos", "strand", "SL02"),
            with = FALSE
        ]
        duplicated <- data.table::rbindlist(
            list(duplicated, duplicated[1L])
        )
        data.table::fwrite(first, "first.tsv", sep = "\t")
        data.table::fwrite(duplicated, "duplicated.tsv", sep = "\t")

        expect_error(
            combineTSSTables(c("first.tsv", "duplicated.tsv")),
            "duplicated.tsv.*duplicate chr.*pos.*strand"
        )
    })
})

test_that("TSS table genomic coordinates are validated", {
    fixture <- tsstable_fixture()[, c("chr", "pos", "strand", "SL01"),
        with = FALSE
    ]

    withr::with_tempdir({
        invalid_pos <- data.table::copy(fixture)
        invalid_pos[, pos := as.numeric(pos)]
        invalid_pos[1L, pos := 1.5]
        data.table::fwrite(invalid_pos, "invalid-pos.tsv", sep = "\t")

        invalid_chr <- data.table::copy(fixture)
        invalid_chr[1L, chr := ""]
        data.table::fwrite(invalid_chr, "invalid-chr.tsv", sep = "\t")

        invalid_strand <- data.table::copy(fixture)
        invalid_strand[1L, strand := "?"]
        data.table::fwrite(
            invalid_strand, "invalid-strand.tsv", sep = "\t"
        )

        expect_error(
            splitTSSTable("invalid-pos.tsv", "SL01"),
            "invalid-pos.tsv.*pos.*positive integers"
        )
        expect_error(
            splitTSSTable("invalid-chr.tsv", "SL01"),
            "invalid-chr.tsv.*chr.*non-empty"
        )
        expect_error(
            splitTSSTable("invalid-strand.tsv", "SL01"),
            "invalid-strand.tsv.*strand.*\\+.*-"
        )
    })
})

test_that("strict chromosome checking rejects incompatible names", {
    fixture <- tsstable_fixture()

    withr::with_tempdir({
        first <- fixture[, c("chr", "pos", "strand", "SL01"),
            with = FALSE
        ]
        second <- fixture[, c("chr", "pos", "strand", "SL02"),
            with = FALSE
        ]
        second[, chr := "I"]
        data.table::fwrite(first, "chr-prefixed.tsv", sep = "\t")
        data.table::fwrite(second, "unprefixed.tsv", sep = "\t")

        expect_error(
            combineTSSTables(c("chr-prefixed.tsv", "unprefixed.tsv")),
            "Chromosome names differ.*chr-prefixed.tsv.*unprefixed.tsv"
        )
    })
})

test_that("allow_missing accepts subsets but not naming-style conflicts", {
    fixture <- tsstable_fixture()

    withr::with_tempdir({
        broader <- fixture[, c("chr", "pos", "strand", "SL01"),
            with = FALSE
        ]
        broader[seq.int(51L, nrow(broader)), chr := "chrII"]
        subset <- fixture[seq_len(50L),
            c("chr", "pos", "strand", "SL02"),
            with = FALSE
        ]
        conflicting <- data.table::copy(subset)
        conflicting[, chr := "I"]
        data.table::fwrite(broader, "broader.tsv", sep = "\t")
        data.table::fwrite(subset, "subset.tsv", sep = "\t")
        data.table::fwrite(conflicting, "conflicting.tsv", sep = "\t")

        expect_no_error(combineTSSTables(
            c("broader.tsv", "subset.tsv"),
            chromosomePolicy = "allow_missing"
        ))
        expect_error(
            combineTSSTables(
                c("broader.tsv", "conflicting.tsv"),
                chromosomePolicy = "allow_missing"
            ),
            "Chromosome naming conflict.*chrI.*I"
        )
    })
})

test_that("all missing input files are reported together", {
    expect_error(
        combineTSSTables(c("missing-one.tsv", "missing-two.tsv")),
        "missing-one.tsv.*missing-two.tsv"
    )
})

test_that("duplicate and missing sample selections are reported clearly", {
    fixture <- tsstable_fixture()
    fixture_path <- system.file(
        "extdata", "example-tss-table.tsv",
        package = "TSSr"
    )

    withr::with_tempdir({
        duplicate_one <- fixture[, c("chr", "pos", "strand", "SL01"),
            with = FALSE
        ]
        duplicate_two <- data.table::copy(duplicate_one)
        data.table::fwrite(duplicate_one, "one.tsv", sep = "\t")
        data.table::fwrite(duplicate_two, "two.tsv", sep = "\t")

        expect_error(
            combineTSSTables(c("one.tsv", "two.tsv")),
            "Duplicate sample names.*SL01"
        )
        expect_error(
            splitTSSTable(fixture_path, c("SL01", "SL01")),
            "samples.*must not contain duplicates.*SL01"
        )
        expect_error(
            splitTSSTable(fixture_path, c("SL01", "unknown")),
            "Samples not found.*unknown"
        )
        expect_error(
            splitTSSTable(fixture_path, "chr"),
            "samples must name sample columns.*chr"
        )
    })
})

test_that("logical options require one non-missing logical value", {
    fixture_path <- system.file(
        "extdata", "example-tss-table.tsv",
        package = "TSSr"
    )

    expect_error(
        splitTSSTable(fixture_path, "SL01", dropZero = NA),
        "dropZero must be TRUE or FALSE"
    )
    expect_error(
        splitTSSTable(fixture_path, "SL01", overwrite = c(TRUE, FALSE)),
        "overwrite must be TRUE or FALSE"
    )
    expect_error(
        combineTSSTables(
            c(fixture_path, fixture_path), overwrite = "yes"
        ),
        "overwrite must be TRUE or FALSE"
    )
})

test_that("metadata-marked processed tables cannot be combined as raw counts", {
    fixture <- tsstable_fixture()

    withr::with_tempdir({
        raw <- fixture[, c("chr", "pos", "strand", "SL01"),
            with = FALSE
        ]
        processed <- fixture[, c("chr", "pos", "strand", "SL02"),
            with = FALSE
        ]
        write_tsstable_with_metadata(
            raw, "raw.tsv",
            c("TSSr-TSStable-Version" = "1", dataType = "raw")
        )
        write_tsstable_with_metadata(
            processed, "processed.tsv",
            c("TSSr-TSStable-Version" = "1", dataType = "processed")
        )

        expect_error(
            combineTSSTables(c("raw.tsv", "processed.tsv")),
            "processed.tsv.*marked as processed.*raw counts"
        )
    })
})

test_that("conflicting genome metadata prevents table combination", {
    fixture <- tsstable_fixture()

    withr::with_tempdir({
        first <- fixture[, c("chr", "pos", "strand", "SL01"),
            with = FALSE
        ]
        second <- fixture[, c("chr", "pos", "strand", "SL02"),
            with = FALSE
        ]
        write_tsstable_with_metadata(
            first, "first.tsv",
            c(dataType = "raw", genomeName = "genome-A")
        )
        write_tsstable_with_metadata(
            second, "second.tsv",
            c(dataType = "raw", genomeName = "genome-B")
        )

        expect_error(
            combineTSSTables(c("first.tsv", "second.tsv")),
            "Conflicting genomeName metadata.*genome-A.*genome-B"
        )
    })
})

test_that("split and combined output files preserve raw-table metadata", {
    fixture <- tsstable_fixture()

    withr::with_tempdir({
        write_tsstable_with_metadata(
            fixture, "source.tsv",
            c(
                "TSSr-TSStable-Version" = "1",
                dataType = "raw",
                genomeName = "sacCer3",
                chromosomes = "chrI"
            )
        )
        splitTSSTable(
            "source.tsv", "SL01", outputFile = "SL01.tsv"
        )
        splitTSSTable(
            "source.tsv", c("SL02", "SL03", "SL04"),
            outputFile = "others.tsv"
        )
        combineTSSTables(
            c("SL01.tsv", "others.tsv"),
            outputFile = "combined.tsv"
        )

        split_header <- readLines("SL01.tsv", n = 5L)
        combined_header <- readLines("combined.tsv", n = 5L)
        for (header in list(split_header, combined_header)) {
            expect_true(any(grepl(
                "^# TSSr-TSStable-Version: 1$", header
            )))
            expect_true(any(grepl("^# dataType: raw$", header)))
            expect_true(any(grepl("^# genomeName: sacCer3$", header)))
            expect_true(any(grepl("^# chromosomes: chrI$", header)))
        }
        expect_identical(
            canonical_tsstable(data.table::fread("combined.tsv")),
            canonical_tsstable(fixture)
        )
    })
})

test_that("exportTSStable identifies raw and processed output metadata", {
    data(exampleTSSr)

    withr::with_tempdir({
        exportTSStable(exampleTSSr, data = "raw", merged = FALSE)
        exportTSStable(exampleTSSr, data = "processed")
        raw_file <- "ALL.samples.TSS.raw.txt"
        processed_file <- "ALL.samples.TSS.processed.txt"
        raw_header <- readLines(raw_file, n = 5L)
        processed_header <- readLines(processed_file, n = 5L)

        expect_true(any(grepl("^# dataType: raw$", raw_header)))
        expect_true(any(grepl(
            "^# dataType: processed$", processed_header
        )))
        expect_true(any(grepl(
            "^# genomeName: BSgenome.Scerevisiae.UCSC.sacCer3$",
            raw_header
        )))
        expect_error(
            combineTSSTables(c(raw_file, processed_file)),
            "processed.*only raw counts"
        )
    })
})

test_that("full chromosome-I raw matrix survives sample split and combine", {
    data(exampleTSSr)
    sample_names <- exampleTSSr@sampleLabels

    withr::with_tempdir({
        exportTSStable(exampleTSSr, data = "raw", merged = FALSE)
        source_file <- "ALL.samples.TSS.raw.txt"
        split_files <- vapply(sample_names, function(sample_name) {
            output_file <- paste0(sample_name, ".tsv")
            splitTSSTable(
                source_file,
                samples = sample_name,
                outputFile = output_file
            )
            output_file
        }, character(1))

        combined <- combineTSSTables(split_files)

        expect_equal(
            canonical_tsstable(combined),
            canonical_tsstable(exampleTSSr@TSSrawMatrix)
        )
    })
})

test_that("duplicate sample columns within one table are rejected", {
    withr::with_tempdir({
        writeLines(
            c(
                "chr\tpos\tstrand\tS1\tS1",
                "chrI\t1\t+\t1\t2"
            ),
            "duplicated-header.tsv"
        )

        expect_error(
            splitTSSTable("duplicated-header.tsv", "S1"),
            "duplicated-header.tsv.*duplicate sample columns.*S1"
        )
    })
})

test_that("invalid raw-count values are rejected", {
    invalid_rows <- c(
        negative = "chrI\t1\t+\t-1",
        missing = "chrI\t1\t+\tNA",
        infinite = "chrI\t1\t+\tInf"
    )

    withr::with_tempdir({
        for (case in names(invalid_rows)) {
            path <- paste0(case, ".tsv")
            writeLines(
                c("chr\tpos\tstrand\tS1", invalid_rows[[case]]),
                path
            )
            expect_error(
                splitTSSTable(path, "S1"),
                "non-negative integer raw counts",
                info = case
            )
        }
    })
})

test_that("table count and output overwrite rules are explicit", {
    fixture_path <- system.file(
        "extdata", "example-tss-table.tsv",
        package = "TSSr"
    )

    expect_error(
        combineTSSTables(fixture_path),
        "At least two TSS table files"
    )
    withr::with_tempdir({
        splitTSSTable(fixture_path, "SL01", outputFile = "selected.tsv")
        expect_error(
            splitTSSTable(
                fixture_path, "SL02", outputFile = "selected.tsv"
            ),
            "Output file already exists"
        )
        expect_no_error(splitTSSTable(
            fixture_path, "SL02", outputFile = "selected.tsv",
            overwrite = TRUE
        ))
        expect_identical(
            names(data.table::fread("selected.tsv")),
            c("chr", "pos", "strand", "SL02")
        )
    })
})
