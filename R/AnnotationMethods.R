###############################################################################
#' Annotate clusters with GFF annotation file.
#'
#' @description Annotates clusters with gene or transcript names from GFF annotation file.
#'
#' @usage    annotateCluster(object,clusters = "consensusClusters",filterCluster = TRUE
#' , filterClusterThreshold = 0.02, annotationType = "genes",upstream=1000
#' , upstreamOverlap = 500,downstream = 0)
#'
#' @param object  A TSSr object
#' @param clusters Clusters to be annotated: "consensusClusters" or "tagClusters".
#' Default is "consensusClusters".
#' @param filterCluster Logical indicating whether clusters downstream of a highly
#' expressed cluster are filtered. Setting filterCluster as "TRUE" would reduce weak
#' clusters brought from recapping, transcriptional or sequencing noise. Default is TRUE.
#' @param filterClusterThreshold  Ignore downstream clusters if signal < filterClusterThreshold*the
#' strongest clusters within the same gene promoter region. Default value = 0.02.
#' @param annotationType  Specify annotation feature to be associated with: "genes"
#' or "transcripts". Default is "genes".
#' @param upstream  Upstream distance to the start position of annotation feature.
#' Default value = 1000.
#' @param upstreamOverlap Upstream distance to the start position of annotation
#' feature if overlapped with the upstream neighboring feature. Default value = 500.
#' @param downstream  Downstream distance to the start position of annotation feature.
#' Default value = 0. Note: if annotationType == "transctipt" or the gene annotations
#' start from transcription start sites (TSSs), the recommended value = 500.
#' @return Large List of elements - one element for each sample
#'
#' @export
#'
#' @examples
#' data(exampleTSSr)
#' annotateCluster(exampleTSSr,clusters = "consensusClusters", filterCluster = TRUE
#' , filterClusterThreshold = 0.02, annotationType = "genes", upstream=1000
#' , upstreamOverlap = 500, downstream = 0)
#'
#'
setGeneric("annotateCluster",function(object, clusters = "consensusClusters"
                                      ,filterCluster = TRUE
                                      ,filterClusterThreshold = 0.02
                                      ,annotationType = "genes"
                                      ,upstream=1000
                                      ,upstreamOverlap = 500
                                      ,downstream = 0)standardGeneric("annotateCluster"))
#' @rdname annotateCluster
#' @export
setMethod("annotateCluster",signature(object = "TSSr"), function(object, clusters, filterCluster
                                                                 , filterClusterThreshold,annotationType
                                                                 , upstream, upstreamOverlap, downstream
){
  message("\nAnnotating...")
  Genome <- .getGenome(object@genomeName)
  sampleLabelsMerged <- object@sampleLabelsMerged
  objName <- deparse(substitute(object))
  refGFF <- object@refSource
  refTable <- object@refTable
  organismName <- object@organismName

  ##check whether there is refTable provided
  if(length(refTable) != 0){
    ref <- refTable
  }else{
    ##check whether there is refSource provided
    if(length(refGFF) == 0 ){
      stop("Please provide correct refSource file!")
    }
    ##define variable as a NULL value
    inCoding = r = f = dominant_tss = NULL
    ##prepare annotation file
    txdb <- suppressWarnings(makeTxDbFromGFF(refGFF, organismName, format = "auto"))
    if(annotationType == "genes"){
      ref <- setDT(as.data.frame(genes(txdb)))
    }else if(annotationType == "transcripts"){
      ref <- setDT(as.data.frame(transcripts(txdb)))
    }
    object@refTable <- ref
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
        temp <- n[list(i)]
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
  names(asned) <- sampleLabelsMerged
  names(unasn) <- sampleLabelsMerged
  object@assignedClusters <- asned
  object@unassignedClusters <- unasn
  if(filterCluster == TRUE){
    names(asn.filtered) <- sampleLabelsMerged
    object@filteredClusters <- asn.filtered}
  assign(objName, object, envir = parent.frame())
})
