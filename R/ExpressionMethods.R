################################################################################################
#' Analysis of gene differential expression.
#'
#' @description Analyzes gene-level differential expression using DESeq2 method (Love et al., 2014).
#' @usage deGene(object,comparePairs=list(c("control","treat")), pval = 0.01,
#'  useMultiCore=FALSE, numCores = NULL)
#'
#' @param object A TSSr object.
#' @param comparePairs Specified list of sample pairs for comparison with DESeq2 method.
#' @param pval Genes with adjusted p value >= pVal will be returned. Default value = 0.01.
#' @param useMultiCore Logical indicating whether multiple cores are used (TRUE) or not (FALSE). Default is FALSE.
#' @param numCores Number of cores are used in clustering step. Used only if useMultiCore = TRUE. Default is NULL.
#'
#'
#' @return
#' @export
#'
#' @examples
#' deGene(exampleTSSr,comparePairs=list(c("control","treat")), pval = 0.01)
setGeneric("deGene",function(object, comparePairs=list(c("control","treat")), pval=0.01, useMultiCore=FALSE, numCores = NULL)standardGeneric("deGene"))
#' @rdname deGene
#' @export
setMethod("deGene",signature(object = "TSSr"), function(object, comparePairs, pval, useMultiCore, numCores){
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
  object@DEtables <- D
  assign(objName, object, envir = parent.frame())
})
