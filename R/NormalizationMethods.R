################################################################################################
#' Normalize raw TSS counts
#'
#' @description Normalizes raw TSS counts in all samples by tags per million (TPM)
#'
#' @usage normalizeTSS(object)
#'
#' @param object A TSSr object.
#'
#' @export
#'
#' @examples
#' data(exampleTSSr)
#' exampleTSSr <- mergeSamples(exampleTSSr)
#' exampleTSSr <- normalizeTSS(exampleTSSr)
setGeneric(
    "normalizeTSS",
    function(object) standardGeneric("normalizeTSS"),
    signature = "object"
)
#' @rdname normalizeTSS
#' @return A modified TSSr object with updated \code{TSSprocessedMatrix}
#'   slot containing normalized TPM values. The input object is not modified;
#'   assign the returned object to retain the changes.
#' @export
setMethod("normalizeTSS", signature(object = "TSSr"), function(object) {
    message("\nNormalizing TSS matrix...")
    ## initialize values
    sampleLabelsMerged <- object@sampleLabelsMerged

    tss.dt <- object@TSSprocessedMatrix
    if (isTRUE(.isNormalized(object))) {
        stop("\tStopping... data is already normalized")
    }

    library.size <- object@librarySizes
    # if library size is empty, get library size
    # calculate size of genome
    # genomeSize <- 0
    # for (chrom in seq(Genome)) {
    #  genomeSize <- genomeSize + length(Genome[[chrom]])
    # }
    ## normalize tss data
    tss.new <- lapply(as.list(seq(sampleLabelsMerged)), function(i) {
        temp <- tss.dt[, .SD, .SDcols = sampleLabelsMerged[i]]
        sizePerMillion <- library.size[i] / 1e6
        setnames(temp, colnames(temp)[[1]], "tags")
        temp[, tags := round(tags / sizePerMillion, 6)]
        setnames(temp, colnames(temp)[[1]], sampleLabelsMerged[i])
        return(temp)
    })
    re <- NULL
    for (i in seq(sampleLabelsMerged)) {
        re <- cbind(re, tss.new[[i]])
    }
    re <- cbind(tss.dt[, c(1, 2, 3)], re)
    setorder(re, "strand", "chr", "pos")
    object@TSSprocessedMatrix <- re
    object@normalizationStatus <- "normalized"
    return(object)
})
