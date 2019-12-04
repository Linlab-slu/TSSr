################################################################################################
setGeneric("shapeCluster",function(object,...)standardGeneric("shapeCluster"))
setMethod("shapeCluster","TSSr", function(object, clusters = "consensusClusters", method = "PSS", useMultiCore= FALSE, numCores = NULL
){
  message("\nCalculating ", clusters," shape with ",method," method...")
  
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
    library(parallel)
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