################################################################################################
setGeneric("deGene",function(object,...)standardGeneric("deGene"))
setMethod("deGene","TSSr", function(object
                                    ,sampleOne
                                    ,sampleTwo
                                    ,pval=0.01
){
  ##initialize data
  objName <- deparse(substitute(myTSSr))
  sampleLabels <- object@sampleLabels
  sampleLabelsMerged <- object@sampleLabelsMerged
  cx <- object@assignedClusters[[sampleOne]]
  cy <- object@assignedClusters[[sampleTwo]]
  tss.raw <- object@TSSrawMatrix
  mergeIndex <- object@mergeIndex
  samplex <- sampleLabels[which(mergeIndex ==which(sampleLabelsMerged == sampleOne))]
  sampley <- sampleLabels[which(mergeIndex ==which(sampleLabelsMerged == sampleTwo))]
  DE.dt <- .deseq2(cx,cy,tss.raw,samplex, sampley, sampleOne, sampleTwo)
  DE.sig <- subset(DE.dt, padj < pval)
  DE.dt$gene <- row.names(DE.dt)
  DE.sig$gene <- row.names(DE.sig)
  DE.dt <- DE.dt[,c(ncol(DE.dt), 1:(ncol(DE.dt)-1))]
  DE.sig <- DE.sig[,c(ncol(DE.sig), 1:(ncol(DE.sig)-1))]
  setDT(DE.dt)
  setDT(DE.sig)
  DE <- list("DEtable" = DE.dt, "DEsig" = DE.sig)
  cat("\n")
  object@DEtables <- DE
  assign(objName, object, envir = parent.frame())
})