#' Select genes which have core promoter shift across different experiments.
#'
#' @description Selects genes which have multiple core promoters and undergo core promoter
#' shifting across different experiments. Generates gene list with Ds (degree of shift)
#' value (Lu et al., 2019), p value and adjusted p value.
#'
#' @usage shiftPromoter(object,comparePairs=list(c("control","treat")), pval = 0.01)
#'
#' @param object A TSSr object.
#' @param comparePairs Specified list of sample pairs for comparison.
#' @param pVal Genes with adjusted p value >= pVal will be returned. Default value = 0.01.
#'
#' @return
#' @export
#'
#' @examples
#' shiftPromoter(myTSSr,comparePairs=list(c("control","treat")), pval = 0.01)
setGeneric("shiftPromoter",function(object,...)standardGeneric("shiftPromoter"))
#' @rdname shiftPromoter
#' @export
setMethod("shiftPromoter",signature(object = "TSSr"), function(object
                                           ,comparePairs
                                           ,pval
){
  ##initialize data
  message("\nCalculating core promoter shifts...")
  objName <- deparse(substitute(object))
  sampleLabelsMerged <- object@sampleLabelsMerged

  D <- lapply(as.list(seq(comparePairs)), function(i){
    sampleOne <- comparePairs[[i]][1]
    sampleTwo <- comparePairs[[i]][2]
    cx <- object@assignedClusters[[sampleOne]]
    cy <- object@assignedClusters[[sampleTwo]]
    tss.raw <- object@TSSrawMatrix
    librarySizex <- object@librarySizes[which(sampleLabelsMerged == sampleOne)]
    librarySizey <- object@librarySizes[which(sampleLabelsMerged == sampleTwo)]
    DS <- .Ds(cx,cy, librarySizex, librarySizey, useRawCount = TRUE, pval)
    return(DS)
  })
  D.names <- sapply(as.list(seq(comparePairs)), function(i){
    paste0(comparePairs[[i]][1],"_VS_",comparePairs[[i]][2], sep ="")
  })
  names(D) <- D.names
  object@PromoterShift <- D
  assign(objName, object, envir = parent.frame())
})
