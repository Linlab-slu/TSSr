################################################################################################
setGeneric("annotateCluster",function(object,...)standardGeneric("annotateCluster"))
setMethod("annotateCluster","TSSr", function(object
                                             ,clusters = "consensusClusters"
                                             ,filterCluster = TRUE
                                             ,filterClusterThreshold = 0.02
                                             ,annotationType = "genes"
                                             ,upstream=1000
                                             ,upstreamOverlap = 500
                                             ,downstream = 0
){
  message("\nAnnotating...")
  Genome <- .getGenome(object@genomeName)
  sampleLabelsMerged <- object@sampleLabelsMerged
  objName <- deparse(substitute(object))
  refGFF <- object@refSource
  organismName <- object@organismName
  
  ##prepare annotation file
  txdb <- suppressWarnings(makeTxDbFromGFF(refGFF, organismName, format = "auto"))
  if(annotationType == "genes"){
    ref <- setDT(as.data.frame(genes(txdb)))
  }else if(annotationType == "transcripts"){
    ref <- setDT(as.data.frame(transcripts(txdb)))
  }
  
  ##prepare clusters
  if(clusters == "tagClusters"){
    cs.dt <- object@tagClusters
  }else if(clusters == "consensusClusters"){
    cs.dt <- object@consensusClusters
  }
  
  ##
  asn <- lapply(as.list(seq(sampleLabelsMerged)), function(i){
    cs.temp <- cs.dt[[sampleLabelsMerged[i]]]
    cs.asn <- .assign2gene(cs.temp,ref,upstream, upstreamOverlap, downstream,filterCluster)
    return(cs.asn)
  })
  names(asn) <- sampleLabelsMerged
  ##subset assigned
  asned <- lapply(as.list(seq(sampleLabelsMerged)), function(i){
    cs <- asn[[sampleLabelsMerged[i]]]
    cs <- cs[!is.na(cs$gene),]
    setDT(cs)
    return(cs)
  })
  ##subset unassigned
  unasn <- lapply(as.list(seq(sampleLabelsMerged)), function(i){
    cs <- asn[[sampleLabelsMerged[i]]]
    cs <- cs[is.na(cs$gene),]
    setDT(cs)
    return(cs)
  })
  ##filter clusters
  if(filterCluster == "TRUE"){
    asn.filtered <- lapply(as.list(seq(sampleLabelsMerged)), function(i){
      cs <- asn[[sampleLabelsMerged[i]]]
      m <- cs[is.na(gene) & is.na(inCoding),]
      n <- cs[!is.na(gene) | !is.na(inCoding),]
      n[,inCoding := ifelse(!is.na(inCoding) & !is.na(gene), NA, inCoding)]
      n[, r := ifelse(is.na(gene), inCoding, gene)]
      setkey(n, r)
      new <- lapply(as.list(n[,unique(r)]), function(i){
        temp <- n[.(i)]
        max.tags <- temp[,max(tags)]
        if(temp[1,strand] == "+"){
          temp[,f := ifelse(dominant_tss > temp[which.max(tags),dominant_tss] & tags < max.tags * filterClusterThreshold, 0,1)]
        }else{
          temp[,f := ifelse(dominant_tss < temp[which.max(tags),dominant_tss] & tags < max.tags * filterClusterThreshold, 0,1)]
        }
        return(temp)
      })
      new <- rbindlist(new)
      new <- new[f==1,]
      cs <- rbind(m[,1:12],new[,1:12])
      return(cs)
    })
  }
  cat("\n")
  names(asned) <- sampleLabelsMerged
  names(unasn) <- sampleLabelsMerged
  object@assignedClusters <- asned
  object@unassignedClusters <- unasn
  if(filterCluster == TRUE){
    names(asn.filtered) <- sampleLabelsMerged
    object@filteredClusters <- asn.filtered}
  assign(objName, object, envir = parent.frame())
})