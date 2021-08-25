#' Identification of enhancers
#'
#' @description Calculates enhancer candidates, which are caracterized by 
#' bidirectional clusters as described in Andersson et al. 2014.
#'
#' @usage callEnhancer(object, flanking = 400)
#'
#' @param object A TSSr object.
#' @param flanking The flanking region range where bidirectional clusters 
#' composing a enhancer candidate. Default is 400.
#' @return Large List of elements - one element for each sample
#'
#' @export
#'
#' @examples
#' \donttest{
#'  data(exampleTSSr)
#' 	#callEnhancer(exampleTSSr,flanking = 400)
#' }
setGeneric("callEnhancer",function(object, flanking = 400)standardGeneric("callEnhancer"))
#' @rdname callEnhancer
#' @export
setMethod("callEnhancer",signature(object = "TSSr"), 
          function(object, flanking){
  message("\nCalculating enhancers...")
  ##define variable as a NULL value
  inCoding = dominant_tss = strand.m = strand.p = cluster = D = chr= NULL
  dominant_tss.m = dominant_tss.p = tags.p = tags.m = NULL
  ##initialize data
  cs.dt <- object@unassignedClusters
  
  if(length(cs.dt) == 0){
    stop("Warning! Clusters must be annotated before calling enhancers.")
  }
  
  sampleLabelsMerged <- object@sampleLabelsMerged
  objName <- deparse(substitute(object))
  
  cs.en <- lapply(as.list(seq(sampleLabelsMerged)), function(i){
    cs <- cs.dt[[sampleLabelsMerged[i]]]
    cs <- cs[is.na(inCoding)]
    if(nrow(cs) >0){
      cs[,gene:= NULL]
      cs[,inCoding:= NULL]
      setkey(cs,chr)
      ce <- lapply(as.list(as.character(cs[,unique(chr)])), function(x) {
        data <- cs[x]
        # data[, strand_no := ifelse(strand == "+",1,0)]
        setorder(data, dominant_tss)
        ##
        data[, strand.m := ifelse(strand == "-" 
                                  & data.table::shift(strand,1,type = "lead") == "+" 
                                  & data.table::shift(dominant_tss, 1, type = "lead") - dominant_tss <= flanking
                                  ,1,0)]
        data[, strand.p := ifelse(strand == "+" 
                                  & data.table::shift(strand,1,type = "lag") == "-" 
                                  & dominant_tss - data.table::shift(dominant_tss, 1, type = "lag") <= flanking
                                  ,1,0)]
        data <- data[strand.m == 1 | strand.p == 1]
        if(nrow(data) >0){
          en <- lapply(as.list(seq(nrow(data)/2)),function(p){
            temp <- data[c(2*p-1,2*p)]
            data.table(enhancer = p,
                       cluster.m = temp[1,cluster],
                       cluster.p = temp[2,cluster],
                       chr = temp[1,chr],
                       dominant_tss.m = temp[1,dominant_tss],
                       dominant_tss.p = temp[2,dominant_tss],
                       tags.m = temp[1,tags],
                       tags.p = temp[2,tags])
          })
          en <- rbindlist(en)
          en <- en[, D:= (tags.p-tags.m)/(tags.p+tags.m)]
          en[D > -0.8 & D < 0.8]
        }
      })
      ce <- rbindlist(ce)
    }
  })

  names(cs.en) <- sampleLabelsMerged
  object@enhancers <- cs.en
  assign(objName, object, envir = parent.frame())
})
##------------------------------------------------------------------------------
