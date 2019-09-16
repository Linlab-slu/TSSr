################################################################################################
setGeneric("annotateCluster",function(object,...)standardGeneric("annotateCluster"))
setMethod("annotateCluster","TSSr", function(object
                                             ,clusters = "tagClusters"
                                             ,reference
                                             ,organim
                                             ,annotationType = "genes"
                                             ,upstream=1000
                                             ,upstreamOverlap = 500
                                             ,downstream = 0
){
  message("\nAnnotating...")
  Genome <- .getGenome(object@genomeName)
  sampleLabelsMerged <- object@sampleLabelsMerged
  objName <- deparse(substitute(myTSSr))
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
    cs <- cs.dt[[sampleLabelsMerged[i]]]
    cs.asn <- .assign2gene(cs,ref,upstream, upstreamOverlap, downstream)
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
  cat("\n")
  names(asned) <- sampleLabelsMerged
  names(unasn) <- sampleLabelsMerged
  object@assignedClusters <- asned
  object@unassignedClusters <- unasn
  assign(objName, object, envir = parent.frame())
})