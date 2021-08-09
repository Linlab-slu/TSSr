#' Analysis of core promoter shape
#'
#' @description Calculates core promoter shape based on the distributions of TSSs within core
#'  promoters using Shape Index (SI) algorithm (Hoskins et al. 2011) or Promoter Shape Score (PSS)
#'   algorithm (Lu et al. 2019).
#'
#' @usage shapeCluster(object, clusters = "consensusClusters", method = "PSS",
#'  useMultiCore=FALSE, numCores = NULL)
#'
#' @param object A TSSr object.
#' @param clusters Clusters to be used for calculating shape score: "tagClusters" or "consensusClusters".
#'  Default is "consensusClusters".
#' @param method Method to be used for calculating core promoter shape score: "SI" or "PSS". Default is "PSS".
#' @param useMultiCore Logical indicating whether multiple cores are used (TRUE) or not (FALSE). Default is FALSE.
#' @param numCores Number of cores are used in clustering step. Used only if useMultiCore = TRUE. Default is NULL.

#'
#' @export
#'
#' @examples
#' \donttest{
#'  data(exampleTSSr)
#' 	#shapeCluster(exampleTSSr,clusters = "consensusClusters" , method = "PSS")
#' 	#shapeCluster(exampleTSSr,clusters = "tagClusters" , method = "SI")
#' }
setGeneric("shapeCluster",function(object, clusters = "consensusClusters"
                                   , method = "PSS", useMultiCore = FALSE, numCores = NULL)standardGeneric("shapeCluster"))
#' @rdname shapeCluster
#' @export
setMethod("shapeCluster",signature(object = "TSSr"), function(object, clusters, method, useMultiCore, numCores){
  message("\nCalculating ", clusters," shape with ",method," method...")
  ##define variable as a NULL value
  pos = interquantile_width = chr= NULL

  ##initialize data
  tss.dt <- object@TSSprocessedMatrix

  if(clusters == "tagClusters"){
    cs.dt <- object@tagClusters
  }else if(clusters == "consensusClusters"){
    cs.dt <- object@consensusClusters
  }

  sampleLabelsMerged <- object@sampleLabelsMerged
  objName <- deparse(substitute(object))

  if (useMultiCore) {
    if (is.null(numCores)) {
      numCores <- detectCores()
    }
    print(paste("process is running on", numCores, "cores..."))
    ##
    cs.shape <- lapply(as.list(seq(sampleLabelsMerged)), function(i){
      cs <- cs.dt[[sampleLabelsMerged[i]]]
      tss <- tss.dt[,.SD, .SDcols = c("chr","pos","strand",sampleLabelsMerged[i])]
      setnames(tss, colnames(tss)[[4]], "tags")
      tss <- tss[tags >0,]
      ce <- mclapply(seq_len(cs[,.N]), function(x) {
        data <- cs[x,]
        tss.sub <- tss[chr == data$chr & strand == data$strand & pos >= data$q_0.1 & pos <= data$q_0.9,]
        temp <- sum(tss.sub[,tags])
        if(method == "PSS"){
          data$shape.score <- -sum(sapply(tss.sub[,tags],function(y){y/temp*log(y/temp,2)}))*log(data[,interquantile_width],2)
        }else if(method == "SI"){
          data$shape.score <- 2+sum(sapply(tss.sub[,tags],function(y){y/temp*log(y/temp,2)}))
        }else{
          message("\nNo shape method is provided...")
        }
        return(data)
      }, mc.cores = numCores)
      ce <- rbindlist(ce)
    })

  }else{
    cs.shape <- lapply(as.list(seq(sampleLabelsMerged)), function(i){
      cs <- cs.dt[[sampleLabelsMerged[i]]]
      tss <- tss.dt[,.SD, .SDcols = c("chr","pos","strand",sampleLabelsMerged[i])]
      setnames(tss, colnames(tss)[[4]], "tags")
      tss <- tss[tags >0,]
      ce <- lapply(as.list(seq_len(cs[,.N])), function(x) {
        data <- cs[x,]
        tss.sub <- tss[chr == data$chr & strand == data$strand & pos >= data$q_0.1 & pos <= data$q_0.9,]
        temp <- sum(tss.sub[,tags])
        if(method == "PSS"){
          data$shape.score <- -sum(sapply(tss.sub[,tags],function(y){y/temp*log(y/temp,2)}))*log(data[,interquantile_width],2)
        }else if(method == "SI"){
          data$shape.score <- 2+sum(sapply(tss.sub[,tags],function(y){y/temp*log(y/temp,2)}))
        }else{
          message("\nNo shape method is provided...")
        }
        return(data)
      })
      ce <- rbindlist(ce)
    })
  }
  names(cs.shape) <- sampleLabelsMerged
  object@clusterShape <- cs.shape
  assign(objName, object, envir = parent.frame())
})
