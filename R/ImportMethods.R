###############################################################################
#' Precisely identify TSSs from bam files, paired end bam files, bed files,
#' BigWig files, tss files, or tss tables.
#'
#' @description getTSS function is used to precisely identify TSSs from multiple
#' input file formats. The files include users' home-made alignment files (bam format)
#' or downloaded files from public databases. See inputFilesType for details on
#' the supported input file formats.
#'
#' @usage getTSS(object, sequencingQualityThreshold = 10,
#' mappingQualityThreshold = 20, softclippingAllowed = FALSE)
#'
#' @param object A TSSr object.
#' @param sequencingQualityThreshold Used only if inputFilesType == "bam" or
#' "bamPairedEnd", otherwise ignored.
#' @param mappingQualityThreshold Used only if inputFilesType == "bam" or
#' "bamPairedEnd", otherwise ignored.
#' @param softclippingAllowed Used only if inputFilesType == "bam" or
#' "bamPairedEnd". Default is FALSE. When FALSE, TSSr applies G-only
#' uncoded G correction to leading G bases on plus-strand reads and
#' trailing C bases on minus-strand reads as transcript-sense G bases.
#' When TRUE, TSSr uses the aligner's aligned 5' boundary directly and
#' skips uncoded G correction.
#' @return A modified TSSr object with updated \code{TSSrawMatrix},
#'   \code{TSSprocessedMatrix}, and \code{librarySizes} slots. The input
#'   object is not modified; assign the returned object to retain the changes.
#'
#' @export
#'
#' @examples
#' exampleInput <- system.file(
#'     "extdata", "example-tss-table.tsv", package = "TSSr", mustWork = TRUE
#' )
#' importedTSSr <- TSSr(
#'     genomeName = "BSgenome.Scerevisiae.UCSC.sacCer3",
#'     inputFiles = exampleInput,
#'     inputFilesType = "TSStable",
#'     sampleLabels = "example",
#'     sampleLabelsMerged = "example",
#'     mergeIndex = 1
#' )
#' importedTSSr <- getTSS(importedTSSr)
#' head(TSSmatrix(importedTSSr, data = "raw"))
setGeneric("getTSS", function(
  object,
  sequencingQualityThreshold = 10,
  mappingQualityThreshold = 20,
  softclippingAllowed = FALSE
) standardGeneric("getTSS"), signature = "object")
#' @rdname getTSS
#' @export
setMethod("getTSS", signature(object = "TSSr"), function(
  object,
  sequencingQualityThreshold,
  mappingQualityThreshold,
  softclippingAllowed
) {
    ## initialize values
    pos <- NULL
    Genome <- .getGenome(object@genomeName)
    sampleLabels <- object@sampleLabels
    inputFilesType <- object@inputFilesType
    inputFiles <- object@inputFiles

    ## Check if input files exist
    missingFiles <- inputFiles[!file.exists(inputFiles)]
    if (length(missingFiles) > 0) {
        stop("Input file(s) not found: ",
             paste(sQuote(missingFiles), collapse = ", "),
             ". Please check the 'inputFiles' slot of your TSSr object ",
             "and ensure the files exist at the specified paths.")
    }

    if (length(object@sampleLabelsMerged) == 0) {
        object@sampleLabelsMerged <- sampleLabels
    }
    if (inputFilesType == "bam" | inputFilesType == "bamPairedEnd") {
        tss <- .getTSS_from_bam(
            object@inputFiles,
            Genome,
            sampleLabels,
            inputFilesType,
            sequencingQualityThreshold,
            mappingQualityThreshold,
            softclippingAllowed
        )
    } else if (inputFilesType == "bed") {
        tss <- .getTSS_from_bed(object@inputFiles, Genome, sampleLabels)
    } else if (inputFilesType == "BigWig") {
        tss <- .getTSS_from_BigWig(object@inputFiles, Genome, sampleLabels)
    } else if (inputFilesType == "tss") {
        tss <- .getTSS_from_tss(object@inputFiles, sampleLabels)
    } else if (inputFilesType == "TSStable") {
        tss <- .getTSS_from_TSStable(object@inputFiles, sampleLabels)
    }
    tss[, pos := .asIntegerCoordinate(pos)]
    setorder(tss, "strand", "chr", "pos")
    # get library sizes
    object@librarySizes <- colSums(tss[, 4:ncol(tss), drop = FALSE], na.rm = TRUE)

    object@TSSrawMatrix <- tss
    object@TSSprocessedMatrix <- tss
    object@normalizationStatus <- "raw"
    return(object)
})
