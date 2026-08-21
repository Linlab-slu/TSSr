#' TSSr: A package for transcription start site sequencing data analyses.
#'
#' TSSr is designed to analyze transcription start sites (TSSs) and core promoters
#' with most types of 5’end sequencing data, such as cap analysis of gene expression (CAGE)
#' (Takahashi, Lassmann et al. 2012), no-amplification non-tagging CAGE libraries for
#' Illumina next-generation sequencers (nAnT-iCAGE) (Murata, Nishiyori-Sueki et al. 2014),
#' a Super-Low Input Carrier-CAGE (SLIC-CAGE) (Cvetesic, Leitch et al. 2018),
#' NanoCAGE (Cumbie, Ivanchenko et al. 2015), TSS-seq (Malabat, Feuerbach et al. 2015),
#' transcript isoform sequencing (TIF-seq) (Pelechano, Wei et al. 2013),
#' transcript-leaders sequencing (TL-seq) (Arribere and Gilbert 2013),
#' precision nuclear run-on sequencing (PRO-Cap) (Mahat, Kwak et al. 2016),
#' and GRO-Cap/5’GRO-seq (Kruesi, Core et al. 2013).
#'
#' TSSr package provides a comprehensive workflow on TSS data starts from identification
#' of accurate TSS locations, clustering TSSs within small genomic regions corresponding
#' to core promoters, and transcriptional activity quantifications, as well as specialized
#' downstream analyses including core promoter shape, cluster annotation, gene differential
#' expression, core promoter shift. TSSr can take multiple formats of files as input, such
#' as Binary Sequence Alignment Map (BAM) files (single-ended or paired-ended), Browser
#' Extension Data (bed) files, BigWig files, ctss files or tss tables. TSSr also generates
#' various types of TSS or core promoter track files which can be visualized in the UCSC
#' Genome Browser or Integrative Genomics Viewer (IGV). TSSr also exports downstream analyses
#' result tables and plots. Multiple cores are supported on Linux or Mac platforms.
#'
#' @keywords internal
"_PACKAGE"

utils::globalVariables(c("gene", "tags"))

.requireSuggestedPackage <- function(package, feature) {
    if (!requireNamespace(package, quietly = TRUE)) {
        stop(
            feature, " requires the suggested package '", package,
            "'. Install it with BiocManager::install(\"", package, "\").",
            call. = FALSE
        )
    }
    invisible(TRUE)
}

.validateFilePaths <- function(paths, argument, allowEmpty = FALSE,
                               requireSingle = FALSE) {
    if (!is.character(paths)) {
        stop("'", argument, "' must be a character vector of file paths.",
             call. = FALSE)
    }
    if (length(paths) == 0L) {
        if (allowEmpty) {
            return(invisible(paths))
        }
        stop("'", argument, "' must contain at least one existing regular file.",
             call. = FALSE)
    }
    if (requireSingle && length(paths) != 1L) {
        stop("'", argument, "' must identify exactly one existing regular file.",
             call. = FALSE)
    }

    has_path <- !is.na(paths) & nzchar(trimws(paths))
    exists <- rep(FALSE, length(paths))
    exists[has_path] <- file.exists(paths[has_path])
    is_regular <- exists
    is_regular[exists] <- !dir.exists(paths[exists])
    invalid <- !has_path | !is_regular

    if (any(invalid)) {
        display <- paths[invalid]
        display[is.na(display)] <- "<NA>"
        display[display == ""] <- "<empty>"
        stop(
            "'", argument, "' must contain existing regular file path(s). ",
            "Invalid path(s): ", paste(sQuote(display), collapse = ", "), ".",
            call. = FALSE
        )
    }
    invisible(paths)
}

.validateNonEmptyCharacter <- function(value, argument, requireSingle = FALSE) {
    descriptor <- if (requireSingle) {
        "a non-empty character string"
    } else {
        "a non-empty character vector without missing or blank values"
    }
    invalid_value <- !is.character(value) || length(value) == 0L ||
        anyNA(value) || any(!nzchar(trimws(value)))
    invalid_length <- requireSingle && length(value) != 1L
    if (invalid_value || invalid_length) {
        stop("'", argument, "' must be ", descriptor, ".", call. = FALSE)
    }
    invisible(value)
}

.requireWorkflowArtifact <- function(object, slotName, nextStep) {
    value <- methods::slot(object, slotName)
    is_missing <- if (is.data.frame(value)) {
        nrow(value) == 0L || ncol(value) == 0L
    } else {
        length(value) == 0L
    }
    if (is_missing) {
        stop(
            "'", slotName, "' is empty; ", nextStep, ".",
            call. = FALSE
        )
    }
    value
}
