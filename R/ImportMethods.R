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
#' "bamPairedEnd". Default is FALSE.
#' @return A modified TSSr object with updated \code{TSSrawMatrix},
#'   \code{TSSprocessedMatrix}, and \code{librarySizes} slots.
#'
#' @export
#'
#' @examples
#' # getTSS requires input files to exist.
#' # The exampleTSSr object has pre-loaded TSS data.
#' data(exampleTSSr)
#' head(slot(exampleTSSr, "TSSrawMatrix"))
setGeneric("getTSS", function(
  object,
  sequencingQualityThreshold = 10,
  mappingQualityThreshold = 20,
  softclippingAllowed = FALSE
) {
    standardGeneric("getTSS")
})
#' @rdname getTSS
#' @export
setMethod("getTSS", signature(object = "TSSr"), function(
  object,
  sequencingQualityThreshold,
  mappingQualityThreshold,
  softclippingAllowed
) {
    ## initialize values
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
    objName <- deparse(substitute(object))
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
    setorder(tss, "strand", "chr", "pos")
    # get library sizes
    object@librarySizes <- colSums(tss[, 4:ncol(tss), drop = FALSE], na.rm = TRUE)

    object@TSSrawMatrix <- tss
    object@TSSprocessedMatrix <- tss
    assign(objName, object, envir = parent.frame())
})
