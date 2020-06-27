################################################################################################
#' Make consensus clusters across multiple samples.
#'
#' @description Makes consensus clusters from multiple samples in TSSr object and calculates
#' inter-quantile positions within consensus clusters for each sample.
#'
#' @usage consensusCluster(object, dis = 50
#' , useMultiCore=TRUE, numCores = NULL)
#'
#' @param object A TSSr object.
#' @param dis Minimum distance between two peaks to be aggregated together into the same consensus cluster.
#' @param useMultiCore Logical indicating whether multiple cores are used (TRUE) or not (FALSE). Default is FALSE.
#' @param numCores Number of cores are used in clustering step. Used only if useMultiCore = TRUE. Default is NULL.
#'
#' @export
#'
#' @examples
#' \donttest{
#' consensusCluster(exampleTSSr,useMultiCore=FALSE)
#' }
setGeneric("consensusCluster",function(object, dis = 50,useMultiCore=TRUE, numCores = NULL)standardGeneric("consensusCluster"))
#' @rdname consensusCluster
#' @export
setMethod("consensusCluster",signature(object = "TSSr"), function(object, dis, useMultiCore, numCores
){
  message("\nCreating consensus clusters...")

  ##initialize data
  tss.dt <- object@TSSprocessedMatrix

  ##define variable as a NULL value
  dominant_tss = NULL

  sampleLabelsMerged <- object@sampleLabelsMerged
  objName <- deparse(substitute(object))
  cs <- object@tagClusters
  if(length(cs) == 0){
    stop("Error: You must have tagClusters data in order to proceed.")
  }
  ##get consensus peak range
  cx <- cs[[sampleLabelsMerged[1]]]
  colnames(cx)[3:4] <- c("start.c","end.c")
  cx[,start := dominant_tss-round(dis/2)]
  cx[,end := dominant_tss + round(dis/2)]
  gr1 <- makeGRangesFromDataFrame(cx, keep.extra.columns= F)
  gr <- BiocGenerics::union(gr1,gr1)
  for(i in 2:length(sampleLabelsMerged)){
    gr <- .getConsensus(gr, cs[[sampleLabelsMerged[[i]]]], dis)
  }
  gr <- as.data.frame(gr)
  gr[,c(1,5)] <- sapply(gr[,c(1,5)], as.character)
  setDT(gr)
  setnames(gr,colnames(gr)[[1]],"chr")
  setorder(gr, "strand","chr","start")
  gr[, consensusCluster := .I]
  ##get consensus quantile
  cs.consensus <- lapply(as.list(seq(sampleLabelsMerged)), function(i){
    tss.temp <- tss.dt[,.SD, .SDcols = c("chr","pos","strand",sampleLabelsMerged[i])]
    setnames(tss.temp, colnames(tss.temp)[[4]], "tags")
    tss.temp <- tss.temp[tags >0,]
    tc <- cs[[sampleLabelsMerged[[i]]]]
    new <- .getConsensusQuantile(tc, gr, tss.temp,useMultiCore, numCores)
    return(new)
    })
  names(cs.consensus) <- sampleLabelsMerged
  object@consensusClusters <- cs.consensus
  assign(objName, object, envir = parent.frame())
})
