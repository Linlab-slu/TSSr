################################################################################################
#' Cluster TSSs into tag clusters
#'
#' @description Clusters TSSs within small genomic regions into tag clusters (TCs) using "peakclu" method.
#'  "peakclu" method is an implementation of peak-based clustering. The minimum distance of two neighboring
#'   peaks can be specified.
#'
#' @usage clusterTSS(object, method = "peakclu", peakDistance=100,extensionDistance=30
#' , localThreshold = 0.02,clusterThreshold = 1, useMultiCore=FALSE, numCores=NULL)
#'
#'
#' @param object  A TSSr object
#' @param method  Clustering method to be used for clustering: "peakclu". Default is "peakclu".
#' @param peakDistance  Minimum distance of two neighboring peaks. Default value = 100.
#' @param extensionDistance Maximal distance between peak and its neighboring TSS or two
#' neighboring TSSs to be grouped in the same cluster. Default value = 30.
#' @param localThreshold  Ignore downstream TSSs with signal < localThreshold*peak within
#' clusters, which is used to filter TSS signals brought from possible recapping events,
#' or sequencing noise. Default value = 0.02.
#' @param clusterThreshold  Ignore clusters if signal < clusterThreshold. Default value = 1.
#' @param useMultiCore Logical indicating whether multiple cores are used (TRUE) or not (FALSE). Default is FALSE.
#' @param numCores Number of cores are used in clustering step. Used only if useMultiCore = TRUE. Default is NULL.
#'
#' @return
#' @export
#'
#' @examples
#' clusterTSS(myTSSr, method = "peakclu",clusterThreshold = 1, useMultiCore=TRUE, numCores = NULL)

setGeneric("clusterTSS",function(object, method = "peakclu"
                                 ,peakDistance=100,extensionDistance=30
                                 ,localThreshold = 0.02,clusterThreshold = 1
                                 ,useMultiCore=FALSE, numCores=NULL)standardGeneric("clusterTSS"))
#' @rdname clusterTSS
#' @export
setMethod("clusterTSS",signature(object = "TSSr"), function(object, method, peakDistance, extensionDistance
                                                            , localThreshold,clusterThreshold,useMultiCore, numCores){
  message("\nClustering TSS data with ", method, " method...")
  ##initialize values
  Genome <- .getGenome(object@genomeName)
  sampleLabelsMerged <- object@sampleLabelsMerged
  objName <- deparse(substitute(object))

  # initialize data
  tss.dt <- object@TSSprocessedMatrix

  # pass sub datatables to peak-caller and clustering functions
  if (useMultiCore) {
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
          cluster.data <- .clusterByPeak(tss, peakDistance, localThreshold, extensionDistance)
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
          cluster.data <- .clusterByPeak(tss, peakDistance, localThreshold, extensionDistance)
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
