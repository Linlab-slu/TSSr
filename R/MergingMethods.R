################################################################################################
#' Merge TSS samples
#'
#' @description Merges individual samples within TSSr object into specified groups.
#' @usage mergeSamples(object, mergeIndex)
#'
#' @param object A TSSr object
#' @param mergeIndex Integer vector specifying which samples to be merged
#' @export
#'
#' @examples
#' \donttest{
#' mergeSamples(exampleTSSr, mergeIndex = c(1,1,2,2))
#' }
setGeneric("mergeSamples",function(object, mergeIndex = NULL)standardGeneric("mergeSamples"))
#' @rdname mergeSamples
#' @export
setMethod("mergeSamples",signature(object = "TSSr"), function(object, mergeIndex){
  if(is.null(mergeIndex)){
    mergeIndex <- as.integer(object@mergeIndex)
  }else{
    object@mergeIndex <- mergeIndex
  }
  sampleLabels <- object@sampleLabels
  sampleLabelsMerged <- object@sampleLabelsMerged

  tss <- object@TSSrawMatrix

  objName <- deparse(substitute(object))
  if(length(mergeIndex) != length(sampleLabels)){
    stop("Length of mergeIndex must match number of samples.")
  }
  if(length(unique(mergeIndex)) != length(sampleLabelsMerged)){
    stop("Number of unique mergeIndex must match number of sampleLabelsMerged.")
  }

  tss.new <- lapply(as.list(seq(unique(mergeIndex))), function(i){
    tss.sub <- tss[, .SD, .SDcols = sampleLabels[which(mergeIndex == i)]]
    tss.sub[,sampleLabelsMerged[i] := rowSums(tss.sub)]
    return(tss.sub[, .SD, .SDcols =sampleLabelsMerged[i]])
  })
  re <- NULL
  for(i in seq(sampleLabelsMerged)){re <- cbind(re, tss.new[[i]])}
  re <- cbind(tss[,1:3],re)

  #object@mergeIndex <- mergeIndex
  object@TSSprocessedMatrix <- re
  object@librarySizes <- colSums(re[,4:ncol(re), drop = F], na.rm = T)
  assign(objName, object, envir = parent.frame())
})
