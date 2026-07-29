#' Filter raw TSS counts or normalized TSS
#'
#' @description Filters transcriptional or sequencing noise.
#'
#' @usage filterTSS(object, method = "poisson", normalization = TRUE,
#' pVal =0.01, tpmLow = 0.1)
#'
#' @param object A TSSr object.
#' @param method Method to be used for TSS filtering: "poisson" or "TPM". "poisson" can be used
#' only if the input TSS data in raw number of counts.
#' @param normalization Define whether normalization data to TPM.
#' Used only if method = “poisson”. Default is TRUE.
#' @param pVal Used only if method = "poisson". Default value is 0.01.
#' @param tpmLow Used only if method = "TPM". Default value is 0.1. Filtering
#'   removes a nonzero count only when this threshold is greater than that
#'   sample's smallest possible nonzero TPM value, approximately
#'   \code{1e6 / librarySize}.
#' @return A modified TSSr object with updated \code{TSSprocessedMatrix}
#'   slot after filtering. The input object is not modified; assign the
#'   returned object to retain the changes.
#'
#'
#' @export
#'
#' @examples
#' data(exampleTSSr)
#' exampleTSSr <- filterTSS(exampleTSSr, method = "TPM", tpmLow = 2)
setGeneric("filterTSS", function(
  object, method = "poisson", normalization = TRUE,
  pVal = 0.01, tpmLow = 0.1
) standardGeneric("filterTSS"), signature = "object")
#' @rdname filterTSS
#' @export
setMethod("filterTSS", signature(object = "TSSr"), function(object, method, normalization, pVal, tpmLow) {
    ## initialize values
    sampleLabelsMerged <- object@sampleLabelsMerged
    library.size <- object@librarySizes

    tss.dt <- object@TSSprocessedMatrix
    is.normalized <- .isNormalized(object)
    ## define variable as a NULL value
    tags <- NULL

    ## filter tss data
    if (method == "poisson") {
        if (isTRUE(is.normalized)) {
            stop("Raw count data required for poisson method.")
        }
        Genome <- .getGenome(object@genomeName)
        # calculate size of genome
        genomeSize <- 0
        for (chrom in seq_len(length(Genome)))
        {
            genomeSize <- genomeSize + length(Genome[[chrom]])
        }
        message("\nFiltering data with ", method, " method...")
        tss.new <- lapply(as.list(seq(sampleLabelsMerged)), function(i) {
            temp <- tss.dt[, .SD, .SDcols = sampleLabelsMerged[i]]
            setnames(temp, colnames(temp)[[1]], "tags")
            temp <- .filterWithPoisson(temp, library.size[i], genomeSize, pVal)
            if (isTRUE(normalization)) {
                sizePerMillion <- library.size[i] / 1e6
                temp[, tags := round(tags / sizePerMillion, 6)]
            }
            setnames(temp, colnames(temp)[[1]], sampleLabelsMerged[i])
            return(temp)
        })
        re <- NULL
        for (i in seq(sampleLabelsMerged)) {
            re <- cbind(re, tss.new[[i]])
        }
        re <- cbind(tss.dt[, c(1, 2, 3)], re)
        ## removes filtered rows
        re <- re[rowSums(re[, 4:ncol(re)]) > 0, ]
        setorder(re, "strand", "chr", "pos")
        object@TSSprocessedMatrix <- re
        object@normalizationStatus <- if (isTRUE(normalization)) {
            "normalized"
        } else {
            "raw"
        }
    } else if (method == "TPM") {
        message("\nFiltering data with ", method, " method...")
        if (!isTRUE(is.normalized)) {
            stop("Data must be normalized.")
        }
        tss.new <- lapply(as.list(seq(sampleLabelsMerged)), function(i) {
            temp <- tss.dt[, .SD, .SDcols = sampleLabelsMerged[i]]
            setnames(temp, colnames(temp)[[1]], "tags")
            temp[tags < tpmLow, ] <- 0
            setnames(temp, colnames(temp)[[1]], sampleLabelsMerged[i])
            return(temp)
        })
        re <- NULL
        for (i in seq(sampleLabelsMerged)) {
            re <- cbind(re, tss.new[[i]])
        }
        re <- cbind(tss.dt[, c(1, 2, 3)], re)
        ## removes filtered rows
        re <- re[rowSums(re[, 4:ncol(re)]) > 0, ]
        setorder(re, "strand", "chr", "pos")
        object@TSSprocessedMatrix <- re
        object@normalizationStatus <- "normalized"
    } else {
        message("\tNo filtering method is defined...")
    }
    return(object)
})
