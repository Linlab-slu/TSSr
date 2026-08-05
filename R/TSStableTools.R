###############################################################################
.readTSStableMetadata <- function(file) {
    lines <- readLines(file, n = 100L, warn = FALSE)
    if (length(lines) == 0L) {
        return(character())
    }
    first_table_line <- which(!grepl("^#", lines))[1L]
    if (is.na(first_table_line) || first_table_line == 1L) {
        return(character())
    }

    metadata_lines <- lines[seq_len(first_table_line - 1L)]
    matched <- regexec("^#\\s*([^:]+):\\s*(.*)$", metadata_lines)
    fields <- regmatches(metadata_lines, matched)
    fields <- fields[lengths(fields) == 3L]
    if (length(fields) == 0L) {
        return(character())
    }
    keys <- tolower(trimws(vapply(fields, `[[`, character(1), 2L)))
    values <- trimws(vapply(fields, `[[`, character(1), 3L))
    if (anyDuplicated(keys)) {
        stop(file, ": duplicate TSStable metadata keys are not allowed.")
    }
    stats::setNames(values, keys)
}

###############################################################################
.metadataChromosomes <- function(metadata, table) {
    declared <- unname(metadata["chromosomes"])
    if (length(declared) == 1L && !is.na(declared) && nzchar(declared)) {
        chromosomes <- trimws(strsplit(declared, ",", fixed = TRUE)[[1L]])
        return(chromosomes[nzchar(chromosomes)])
    }
    unique(as.character(table[["chr"]]))
}

.TSStableMetadata <- function(dataType, genomeName, chromosomes) {
    metadata <- c(
        "TSSr-TSStable-Version" = "1",
        dataType = dataType
    )
    if (length(genomeName) == 1L && !is.na(genomeName) &&
        nzchar(genomeName)) {
        metadata <- c(metadata, genomeName = genomeName)
    }
    c(
        metadata,
        chromosomes = paste(unique(as.character(chromosomes)), collapse = ",")
    )
}

.rawTSStableMetadata <- function(genomeName, chromosomes) {
    .TSStableMetadata("raw", genomeName, chromosomes)
}

.validateTSStableFlag <- function(value, name) {
    if (!is.logical(value) || length(value) != 1L || is.na(value)) {
        stop(name, " must be TRUE or FALSE.")
    }
    invisible(value)
}

.writeTSStableWithMetadata <- function(x, file, metadata, overwrite) {
    if (length(file) != 1L || is.na(file) || !nzchar(file)) {
        stop("outputFile must be one non-empty file path.")
    }
    if (!dir.exists(dirname(file))) {
        stop("Output directory does not exist: ", dirname(file))
    }
    if (file.exists(file) && !isTRUE(overwrite)) {
        stop("Output file already exists: ", file)
    }

    temporary_table <- tempfile("tsstable-table-")
    temporary_output <- tempfile(
        "tsstable-output-", tmpdir = dirname(file)
    )
    on.exit(unlink(c(temporary_table, temporary_output)), add = TRUE)
    .writeTSStable(x, temporary_table)
    metadata_lines <- paste0(
        "# ", names(metadata), ": ", unname(metadata)
    )
    writeLines(metadata_lines, temporary_output, useBytes = TRUE)
    if (!file.append(temporary_output, temporary_table)) {
        stop("Failed to assemble output TSS table: ", file)
    }
    backup <- character()
    if (file.exists(file)) {
        backup <- tempfile("tsstable-backup-", tmpdir = dirname(file))
        if (!file.rename(file, backup)) {
            stop("Could not prepare existing output for replacement: ", file)
        }
    }
    if (!file.rename(temporary_output, file)) {
        if (length(backup) == 1L) {
            file.rename(backup, file)
        }
        stop("Could not move completed TSS table to: ", file)
    }
    if (length(backup) == 1L) {
        unlink(backup)
    }
    invisible(file)
}

###############################################################################
.readTSStableForManipulation <- function(file) {
    if (length(file) != 1L || is.na(file) || !nzchar(file)) {
        stop("Each TSS table path must be one non-empty string.")
    }
    if (!file.exists(file)) {
        stop("Could not locate TSS table: ", file)
    }

    metadata <- .readTSStableMetadata(file)
    data_type <- unname(metadata["datatype"])
    if (length(data_type) == 1L && !is.na(data_type) &&
        !identical(tolower(data_type), "raw")) {
        stop(
            file, " is marked as ", data_type,
            "; only raw counts can be combined or split."
        )
    }
    table <- data.table::fread(file)
    coordinate_columns <- c("chr", "pos", "strand")
    if (ncol(table) < 4L ||
        !identical(names(table)[seq_len(3L)], coordinate_columns)) {
        stop(
            file, ": the first three columns must be chr, pos, and strand, ",
            "followed by at least one sample column."
        )
    }

    chr <- table[["chr"]]
    if (anyNA(chr) || any(!nzchar(trimws(as.character(chr))))) {
        stop(file, ": chr values must be non-empty.")
    }
    table[, chr := as.character(chr)]
    pos <- table[["pos"]]
    valid_positions <- is.numeric(pos) &&
        all(is.finite(pos)) &&
        all(pos >= 1) &&
        all(pos <= .Machine$integer.max) &&
        all(pos == floor(pos))
    if (!valid_positions) {
        stop(file, ": pos values must be positive integers.")
    }
    table[, pos := as.integer(pos)]
    strand <- table[["strand"]]
    if (anyNA(strand) || any(!strand %in% c("+", "-"))) {
        stop(file, ": strand values must be '+' or '-'.")
    }
    table[, strand := as.character(strand)]

    sample_columns <- names(table)[-(seq_len(3L))]
    if (anyDuplicated(sample_columns)) {
        duplicate_samples <- unique(
            sample_columns[duplicated(sample_columns)]
        )
        stop(
            file, ": duplicate sample columns are not allowed: ",
            paste(duplicate_samples, collapse = ", "), "."
        )
    }
    if (any(duplicated(table, by = coordinate_columns))) {
        stop(
            file, ": duplicate chr, pos, and strand coordinates are not ",
            "allowed."
        )
    }
    for (sample_name in sample_columns) {
        values <- table[[sample_name]]
        valid_counts <- is.numeric(values) &&
            all(is.finite(values)) &&
            all(values >= 0) &&
            all(values == floor(values))
        if (!valid_counts) {
            stop(
                file, ": sample column ", sQuote(sample_name),
                " must contain non-negative integer raw counts."
            )
        }
    }
    attr(table, "tsstable_metadata") <- metadata
    table
}

###############################################################################
.canonicalChromosomeName <- function(x) {
    canonical <- sub("^chr", "", tolower(as.character(x)))
    canonical[canonical %in% c("m", "mt")] <- "mt"
    canonical
}

.validateMissingChromosomes <- function(tables, files) {
    pairs <- utils::combn(seq_along(tables), 2L)
    for (column in seq_len(ncol(pairs))) {
        i <- pairs[1L, column]
        j <- pairs[2L, column]
        left <- unique(as.character(tables[[i]][["chr"]]))
        right <- unique(as.character(tables[[j]][["chr"]]))
        left_canonical <- .canonicalChromosomeName(left)
        right_canonical <- .canonicalChromosomeName(right)
        common <- intersect(left_canonical, right_canonical)
        if (length(common) == 0L) {
            stop(
                "No shared chromosome names between ", files[[i]],
                " and ", files[[j]], "; compatibility cannot be verified."
            )
        }
        for (name in common) {
            left_name <- left[left_canonical == name]
            right_name <- right[right_canonical == name]
            if (!identical(sort(left_name), sort(right_name))) {
                stop(
                    "Chromosome naming conflict between ", files[[i]],
                    " and ", files[[j]], ": ",
                    paste(left_name, collapse = ", "), " versus ",
                    paste(right_name, collapse = ", "), "."
                )
            }
        }
    }
}

###############################################################################
#' Combine raw-count TSS tables
#'
#' @description Combine two or more TSS tables by genomic coordinate. Each
#'   input table may contain one or more sample columns.
#'
#' @param inputFiles Character vector containing at least two TSS table paths.
#' @param outputFile Optional path for the combined table. If \code{NULL}, no
#'   file is written.
#' @param chromosomePolicy Chromosome compatibility policy. The default
#'   \code{"strict"} requires identical declared or observed chromosome sets.
#'   \code{"allow_missing"} permits missing chromosomes but still rejects
#'   inconsistent naming conventions.
#' @param overwrite Logical indicating whether an existing output file may be
#'   replaced.
#'
#' @details Input tables must begin with \code{chr}, \code{pos}, and
#'   \code{strand}, followed by one or more uniquely named sample columns.
#'   Sample values must be finite, non-negative integer raw counts. Each table
#'   may contain one or multiple samples, and any number of tables greater than
#'   one can be combined. Missing coordinate/sample combinations are filled
#'   with zero. Files produced by TSSr carry metadata identifying raw versus
#'   processed values, the genome when known, and the chromosome set. Legacy
#'   files without metadata remain supported after value validation.
#'
#' @return A base \code{data.frame} containing the coordinate union and all
#'   sample columns.
#'
#' @export
#'
#' @examples
#' fixture <- system.file(
#'     "extdata", "example-tss-table.tsv", package = "TSSr"
#' )
#' part1 <- tempfile(fileext = ".tsv")
#' part2 <- tempfile(fileext = ".tsv")
#' splitTSSTable(fixture, c("SL01", "SL02"), part1)
#' splitTSSTable(fixture, c("SL03", "SL04"), part2)
#' combined <- combineTSSTables(c(part1, part2))
#' head(combined)
combineTSSTables <- function(inputFiles, outputFile = NULL,
                             chromosomePolicy = "strict",
                             overwrite = FALSE) {
    .validateTSStableFlag(overwrite, "overwrite")
    chromosomePolicy <- match.arg(
        chromosomePolicy, c("strict", "allow_missing")
    )
    if (length(inputFiles) < 2L) {
        stop("At least two TSS table files are required.")
    }
    if (!is.character(inputFiles) || anyNA(inputFiles) ||
        any(!nzchar(inputFiles))) {
        stop("inputFiles must contain non-empty file paths.")
    }
    missing_files <- inputFiles[!file.exists(inputFiles)]
    if (length(missing_files) > 0L) {
        stop(
            "Could not locate TSS tables: ",
            paste(missing_files, collapse = ", "), "."
        )
    }

    tables <- lapply(inputFiles, .readTSStableForManipulation)
    metadata <- lapply(tables, attr, which = "tsstable_metadata")
    chromosome_sets <- Map(.metadataChromosomes, metadata, tables)
    tables <- lapply(tables, function(x) {
        attr(x, "tsstable_metadata") <- NULL
        x
    })
    genome_names <- vapply(metadata, function(x) {
        value <- unname(x["genomename"])
        if (length(value) == 1L && !is.na(value) && nzchar(value)) {
            value
        } else {
            NA_character_
        }
    }, character(1))
    known_genomes <- unique(stats::na.omit(genome_names))
    if (length(known_genomes) > 1L) {
        stop(
            "Conflicting genomeName metadata: ",
            paste(known_genomes, collapse = ", "), "."
        )
    }
    coordinate_columns <- c("chr", "pos", "strand")
    if (identical(chromosomePolicy, "strict")) {
        chromosome_sets <- lapply(chromosome_sets, function(x) {
            sort(unique(as.character(x)))
        })
        incompatible <- which(!vapply(
            chromosome_sets,
            identical,
            logical(1),
            chromosome_sets[[1L]]
        ))
        if (length(incompatible) > 0L) {
            stop(
                "Chromosome names differ between ", inputFiles[[1L]],
                " and ", inputFiles[[incompatible[[1L]]]], "."
            )
        }
    } else {
        .validateMissingChromosomes(tables, inputFiles)
    }
    sample_names <- unlist(
        lapply(tables, function(x) setdiff(names(x), coordinate_columns)),
        use.names = FALSE
    )
    if (anyDuplicated(sample_names)) {
        duplicate_names <- unique(sample_names[duplicated(sample_names)])
        stop(
            "Duplicate sample names across input TSS tables: ",
            paste(duplicate_names, collapse = ", "), "."
        )
    }

    combined <- Reduce(function(x, y) {
        merge(x, y, by = coordinate_columns, all = TRUE, sort = FALSE)
    }, tables)
    for (sample_name in sample_names) {
        data.table::set(
            combined,
            i = which(is.na(combined[[sample_name]])),
            j = sample_name,
            value = 0
        )
    }
    combined <- combined[
        rowSums(combined[, .SD, .SDcols = sample_names]) > 0
    ]
    data.table::setorderv(combined, coordinate_columns)

    result <- as.data.frame(combined)
    if (!is.null(outputFile)) {
        output_metadata <- .rawTSStableMetadata(
            known_genomes,
            unique(unlist(chromosome_sets, use.names = FALSE))
        )
        .writeTSStableWithMetadata(
            result, outputFile, output_metadata, overwrite
        )
    }
    result
}

###############################################################################
#' Extract samples from a raw-count TSS table
#'
#' @description Select sample columns from a TSS table and optionally remove
#'   positions at which every selected sample has a zero count.
#'
#' @param inputFile Path to one TSS table.
#' @param samples Character vector of sample columns to extract.
#' @param outputFile Optional path for the extracted table. If \code{NULL}, no
#'   file is written.
#' @param dropZero Logical indicating whether rows containing zero for every
#'   selected sample should be removed.
#' @param overwrite Logical indicating whether an existing output file may be
#'   replaced.
#'
#' @details This function selects columns; it does not rename or sum samples.
#'   With \code{dropZero = TRUE}, a row is retained when at least one selected
#'   sample has a non-zero raw count. Input validation and metadata handling are
#'   the same as for \code{combineTSSTables()}.
#'
#' @return A base \code{data.frame} containing coordinates and the selected
#'   samples.
#'
#' @export
#'
#' @examples
#' fixture <- system.file(
#'     "extdata", "example-tss-table.tsv", package = "TSSr"
#' )
#' control <- splitTSSTable(fixture, samples = c("SL01", "SL02"))
#' head(control)
splitTSSTable <- function(inputFile, samples, outputFile = NULL,
                          dropZero = TRUE, overwrite = FALSE) {
    .validateTSStableFlag(dropZero, "dropZero")
    .validateTSStableFlag(overwrite, "overwrite")
    table <- .readTSStableForManipulation(inputFile)
    metadata <- attr(table, "tsstable_metadata")
    chromosomes <- .metadataChromosomes(metadata, table)
    attr(table, "tsstable_metadata") <- NULL
    coordinate_columns <- c("chr", "pos", "strand")
    if (!is.character(samples) || length(samples) == 0L || anyNA(samples) ||
        any(!nzchar(samples))) {
        stop("samples must contain one or more non-empty sample names.")
    }
    if (anyDuplicated(samples)) {
        duplicate_samples <- unique(samples[duplicated(samples)])
        stop(
            "samples must not contain duplicates: ",
            paste(duplicate_samples, collapse = ", "), "."
        )
    }
    reserved_samples <- intersect(samples, coordinate_columns)
    if (length(reserved_samples) > 0L) {
        stop(
            "samples must name sample columns, not coordinate columns: ",
            paste(reserved_samples, collapse = ", "), "."
        )
    }
    missing_samples <- setdiff(samples, names(table))
    if (length(missing_samples) > 0L) {
        stop(
            "Samples not found in TSS table: ",
            paste(missing_samples, collapse = ", "), "."
        )
    }

    result <- table[, c(coordinate_columns, samples), with = FALSE]
    if (isTRUE(dropZero)) {
        result <- result[rowSums(result[, .SD, .SDcols = samples]) > 0]
    }
    result <- as.data.frame(result)

    if (!is.null(outputFile)) {
        genome_name <- unname(metadata["genomename"])
        .writeTSStableWithMetadata(
            result,
            outputFile,
            .rawTSStableMetadata(genome_name, chromosomes),
            overwrite
        )
    }
    result
}
