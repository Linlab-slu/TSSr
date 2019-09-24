################################################################################################
setGeneric("consensusCluster",function(object,...)standardGeneric("consensusCluster"))
setMethod("consensusCluster","TSSr", function(object, data = "filtered", dis = 50
){
  message("\nCreating consensus clusters...")
  
  ##initialize data
  if(data == "raw"){
    tss.dt <- object@TSSrawMatrix
  }else if(data == "filtered"){
    tss.dt <- object@TSSfilteredMatrix
  }
  sampleLabelsMerged <- object@sampleLabelsMerged
  objName <- deparse(substitute(myTSSr))
  cs <- object@tagClusters
  ##get consensus peak range
  cx <- cs[[sampleLabelsMerged[1]]]
  colnames(cx)[3:4] <- c("start.c","end.c")
  cx[,start := dominant_tss-round(dis/2)]
  cx[,end := dominant_tss + round(dis/2)]
  gr1 <- makeGRangesFromDataFrame(cx, keep.extra.columns= F)
  gr <- union(gr1,gr1)
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
    new <- .getConsensusQuantile(tc, gr, tss.temp)
    })
  names(cs.consensus) <- sampleLabelsMerged
  cat("\n")
  object@consensusClusters <- cs.consensus
  assign(objName, object, envir = parent.frame())
})