################################################################################################
setGeneric("shiftPromoter",function(object,...)standardGeneric("shiftPromoter"))
setMethod("shiftPromoter","TSSr", function(object
                                           ,comparePairs
){
  ##initialize data
  message("\nCalculating core promoter shifts...")
  objName <- deparse(substitute(myTSSr))
  sampleLabelsMerged <- object@sampleLabelsMerged
  
  D <- lapply(as.list(seq(comparePairs)), function(i){
    sampleOne <- comparePairs[[i]][1]
    sampleTwo <- comparePairs[[i]][2]
    cx <- object@assignedClusters[[sampleOne]]
    cy <- object@assignedClusters[[sampleTwo]]
    tss.raw <- object@TSSrawMatrix
    librarySizex <- object@librarySizes[which(sampleLabelsMerged == sampleOne)] 
    librarySizey <- object@librarySizes[which(sampleLabelsMerged == sampleTwo)] 
    DS <- .Ds(cx,cy, librarySizex, librarySizey, useRawCount = TRUE)
    return(DS)
  })
  D.names <- sapply(as.list(seq(comparePairs)), function(i){
    paste0(comparePairs[[i]][1]," VS. ",comparePairs[[i]][2], sep ="")
  })
  
  cat("\n")
  object@PromoterShift <- D
  assign(objName, object, envir = parent.frame())
})