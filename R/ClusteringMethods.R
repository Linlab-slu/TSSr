################################################################################################
##
################################################################################################
setGeneric("clusterTSS",function(object,...)standardGeneric("clusterTSS"))
setMethod("clusterTSS","TSSr", function(object,data = "filtered", method = "peakclu"
                                        ,peakDistance=100, localThreshold = 0.02
                                        ,extensionDistance=30,clusterThreshold = 1
                                        , nonOverLapping=TRUE
                                        ,useMultiCore=FALSE, numCores=NULL
){
  message("\nClustering TSS data with ", method, " method...")
  ##initialize values
  Genome <- .getGenome(object@genomeName)
  sampleLabelsMerged <- object@sampleLabelsMerged
  objName <- deparse(substitute(object))
  
  ##pick which data will be used for clustering
  if(data == "raw"){
    tss.dt <- object@TSSrawMatrix
  }else if(data == "filtered"){
    tss.dt <- object@TSSfilteredMatrix
  }
  # pass sub datatables to peak-caller and clustering functions
  if (useMultiCore) {
    library(parallel)
    if (is.null(numCores)) {
      numCores <- detectCores()
    }
    print(paste("process is running on", numCores, "cores..."))
    ##separate tss table by sampleLables
    cs <- lapply(as.list(seq(sampleLabelsMerged)), function(i){
      temp <- tss.dt[,.SD, .SDcols = c("chr","pos","strand",sampleLabelsMerged[i])]
      setnames(temp, colnames(temp)[[4]], "tags")
      temp <- temp[tags >0,]
      temp[, "subset" := paste0(chr,"_",strand)]
      setkey(temp,subset)
      clusters <- mclapply(as.list(temp[,unique(subset)]), function(x) {
        tss <- temp[.(x)]
        setkey(tss, NULL)
        setorder(tss, pos)
        if(method == "peakclu"){
          cluster.data <- .clusterByPeak(tss, peakDistance, localThreshold, extensionDistance, nonOverLapping)
        }
      }, mc.cores = numCores)
      tss.clusters <- rbindlist(clusters, use.names=TRUE, fill=TRUE)
      tss.clusters <- tss.clusters[tags >clusterThreshold,]
      setorder(tss.clusters, "strand","chr","start")
      tss.clusters[, cluster := .I]
      return(tss.clusters)
    })
  }else{
    cs <- lapply(as.list(seq(sampleLabelsMerged)), function(i){
      temp <- tss.dt[,.SD, .SDcols = c("chr","pos","strand",sampleLabelsMerged[i])]
      setnames(temp, colnames(temp)[[4]], "tags")
      temp <- temp[tags >0,]
      temp[, "subset" := paste0(chr,"_",strand)]
      setkey(temp,subset)
      clusters <- lapply(as.list(temp[,unique(subset)]), function(x) {
        tss <- temp[.(x)]
        setkey(tss, NULL)
        setorder(tss, pos)
        if(method == "peakclu"){
          cluster.data <- .clusterByPeak(tss, peakDistance, localThreshold, extensionDistance, nonOverLapping)
        }
      })
      tss.clusters <- rbindlist(clusters, use.names=TRUE, fill=TRUE)
      tss.clusters <- tss.clusters[tags >clusterThreshold,]
      setorder(tss.clusters, "strand","chr","start")
      tss.clusters[, cluster := .I]
      return(tss.clusters)
    })
  }
  names(cs) <- sampleLabelsMerged
  object@tagClusters <- cs
  assign(objName, object, envir = parent.frame())
})
