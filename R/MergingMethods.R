################################################################################################
setGeneric("mergeSamples",function(object,...)standardGeneric("mergeSamples"))
setMethod("mergeSamples","TSSr", function(object
                                      ,sampleLabels
                                      ,sampleLabelsMerged
                                      ,mergeIndex
){
  sampleLabels <- object@sampleLabels
  sampleLabelsMerged <- object@sampleLabelsMerged
  mergeIndex <- as.integer(object@mergeIndex)
  
  tss <- object@TSSprocessedMatrix
  
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
  #object@librarySizes <- as.integer(colSums(re[,4:ncol(re), drop = F], na.rm = T))
  assign(objName, object, envir = parent.frame())
})