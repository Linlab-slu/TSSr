#' Select genes which have core promoter shift across different experiments.
#'
#' @description Selects genes which have multiple core promoters and undergo core promoter
#' shifting across different experiments. Generates gene list with Ds (degree of shift)
#' value (Lu et al., 2019), p value and adjusted p value.
#'
#' @usage shiftPromoter(object, comparePairs, pval=0.01)
#'
#' @param object A TSSr object.
#' @param comparePairs Specified list of sample pairs for comparison.
#' @param pval Genes with adjusted p value >= pval will be returned. Default value = 0.01.
#' @return A modified TSSr object with updated \code{PromoterShift} slot. The
#'   input object is not modified; assign the returned object to retain the
#'   changes.
#'
#' @export
#'
#' @examples
#' data(exampleTSSr)
#' # Promoter shift results are pre-computed in the @PromoterShift slot
#' head(slot(exampleTSSr, "PromoterShift")[["control_VS_treat"]])
setGeneric(
    "shiftPromoter",
    function(object, comparePairs, pval = 0.01) standardGeneric("shiftPromoter"),
    signature = "object"
)
#' @rdname shiftPromoter
#' @export
setMethod("shiftPromoter", signature(object = "TSSr"), function(
  object,
  comparePairs,
  pval
) {
    ## initialize data
    message("\nCalculating core promoter shifts...")
    sampleLabelsMerged <- object@sampleLabelsMerged

    D <- lapply(as.list(seq(comparePairs)), function(i) {
        sampleOne <- comparePairs[[i]][1]
        sampleTwo <- comparePairs[[i]][2]
        cx <- object@assignedClusters[[sampleOne]]
        cy <- object@assignedClusters[[sampleTwo]]
        tss.raw <- object@TSSrawMatrix
        librarySizex <- object@librarySizes[which(sampleLabelsMerged == sampleOne)]
        librarySizey <- object@librarySizes[which(sampleLabelsMerged == sampleTwo)]
        DS <- .Ds(cx, cy, librarySizex, librarySizey, useRawCount = TRUE, pval)
        return(DS)
    })
    D.names <- vapply(as.list(seq(comparePairs)), function(i) {
        paste0(comparePairs[[i]][1], "_VS_", comparePairs[[i]][2], sep = "")
    }, character(1))
    names(D) <- D.names
    object@PromoterShift <- D
    return(object)
})
