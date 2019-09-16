################################################################################################
setGeneric("shiftPromoter",function(object,...)standardGeneric("shiftPromoter"))
setMethod("shiftPromoter","TSSr", function(object
                                           ,sampleOne
                                           ,sampleTwo
){
  ##initialize data
  objName <- deparse(substitute(myTSSr))
  sampleLabelsMerged <- object@sampleLabelsMerged
  cx <- object@assignedClusters[[sampleOne]]
  cy <- object@assignedClusters[[sampleTwo]]
  tss.raw <- object@TSSrawMatrix
  mergeIndex <- object@mergeIndex
  samplex <- sampleLabels[which(mergeIndex ==which(sampleLabelsMerged == sampleOne))]
  sampley <- sampleLabels[which(mergeIndex ==which(sampleLabelsMerged == sampleTwo))]
  

  
  
  cat("\n")
  object@PromoterShift <- Ds
  assign(objName, object, envir = parent.frame())
})