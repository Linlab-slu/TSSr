################################################################################################
setGeneric("deGene",function(object,...)standardGeneric("deGene"))
setMethod("deGene","TSSr", function(object
                                    ,comparePairs
                                    ,pval=0.01
                                    ,useMultiCore
                                    ,numCores
){
  ##initialize data
  message("\nCalculating gene differential expression...")
  objName <- deparse(substitute(object))
  sampleLabels <- object@sampleLabels
  sampleLabelsMerged <- object@sampleLabelsMerged
  
  D <- lapply(as.list(seq(comparePairs)), function(i){
    sampleOne <- comparePairs[[i]][1]
    sampleTwo <- comparePairs[[i]][2]
    cx <- object@assignedClusters[[sampleOne]]
    cy <- object@assignedClusters[[sampleTwo]]
    tss.raw <- object@TSSrawMatrix
    mergeIndex <- object@mergeIndex
    samplex <- sampleLabels[which(mergeIndex ==which(sampleLabelsMerged == sampleOne))]
    sampley <- sampleLabels[which(mergeIndex ==which(sampleLabelsMerged == sampleTwo))]
    DE.dt <- .deseq2(cx,cy,tss.raw,samplex, sampley, sampleOne, sampleTwo,useMultiCore, numCores)
    DE.sig <- subset(DE.dt, padj < pval)
    DE.dt$gene <- row.names(DE.dt)
    DE.sig$gene <- row.names(DE.sig)
    DE.dt <- DE.dt[,c(ncol(DE.dt), 1:(ncol(DE.dt)-1))]
    DE.sig <- DE.sig[,c(ncol(DE.sig), 1:(ncol(DE.sig)-1))]
    setDT(DE.dt)
    setDT(DE.sig)
    DE <- list("DEtable" = DE.dt, "DEsig" = DE.sig)
    return(DE)
  })
  D.names <- sapply(as.list(seq(comparePairs)), function(i){
    paste0(comparePairs[[i]][1],"_VS_",comparePairs[[i]][2], sep ="")
  })
  names(D) <- D.names
  
  
  
  cat("\n")
  object@DEtables <- D
  assign(objName, object, envir = parent.frame())
})