#' Filter raw TSS counts or normalized TSS
#'
#' @description Filters transcriptional or sequencing noise.
#'
#' @usage filterTSS(object, method = "poisson", normalization = TRUE,
#' pVal =0.01, tpmLow = 0.1)
#'
#' @param object A TSSr object.
#' @param method Method to be used for TSS filtering: "poisson" or "TPM". "poisson" can be used
#' only if the input TSS data in raw number of counts.
#' @param normalization Define whether normalization data to TPM. Used only if method = “poisson”. Default is TRUE.
#' @param pVal Used only if method = "poisson". Default value is 0.01.
#' @param tpmLow Used only if method = "TPM". Default value is 0.1.
#'
#'
#' @export
#'
#' @examples
#' \donttest{
#' filterTSS(exampleTSSr, method = "TPM", tpmLow=0.1)
#' filterTSS(exampleTSSr, method = "poisson", pVal = 0.01)
#' }
setGeneric("filterTSS",function(object, method = "poisson", normalization = TRUE
                                , pVal =0.01, tpmLow = 0.1)standardGeneric("filterTSS"))
#' @rdname filterTSS
#' @export
setMethod("filterTSS",signature(object = "TSSr"), function(object, method, normalization, pVal, tpmLow){
  ##initialize values
  Genome <- .getGenome(object@genomeName)
  sampleLabelsMerged <- object@sampleLabelsMerged
  objName <- deparse(substitute(object))
  library.size <- object@librarySizes
  # calculate size of genome
  genomeSize <- 0
  for (chrom in 1:length(Genome)) {
    genomeSize <- genomeSize + length(Genome[[chrom]])
  }
  tss.dt <- object@TSSprocessedMatrix
  ##define variable as a NULL value
  tags = NULL

  ##filter tss data
  if(method  == "poisson"){
    message("\nFiltering data with ", method," method...")
    if(any(tss.dt[,4] > 0 & tss.dt[,4] < 1)){
      stop("Warning! Raw count data required for poisson method.")
    }
    tss.new <- lapply(as.list(seq(sampleLabelsMerged)), function(i){
      temp <- tss.dt[,.SD, .SDcols = sampleLabelsMerged[i]]
      setnames(temp, colnames(temp)[[1]], "tags")
      temp <- .filterWithPoisson(temp, library.size[i], genomeSize, pVal)
      if(normalization == "TRUE"){
        sizePerMillion <- library.size[i] / 1e6
        temp[, tags := round(tags / sizePerMillion, 6)]
      }
      setnames(temp, colnames(temp)[[1]], sampleLabelsMerged[i])
      return(temp)
    })
    re <- NULL
    for(i in seq(sampleLabelsMerged)){re <- cbind(re, tss.new[[i]])}
    re <- cbind(tss.dt[,1:3],re)
    ##removes filtered rows
    re <- re[rowSums(re[,4:ncol(re)]) >0,]
    setorder(re, "strand","chr","pos")
    object@TSSprocessedMatrix <- re
    assign(objName, object, envir = parent.frame())
  }else if(method  == "TPM"){
    message("\nFiltering data with ", method," method...")
    if(any(tss.dt[,4] > 0 & tss.dt[,4] < 1) == FALSE){
      stop("Warning! Data must be normalized.")
    }
    tss.new <- lapply(as.list(seq(sampleLabelsMerged)), function(i){
      temp <- tss.dt[,.SD, .SDcols = sampleLabelsMerged[i]]
      setnames(temp, colnames(temp)[[1]], "tags")
      temp[tags < tpmLow,]=0
      setnames(temp, colnames(temp)[[1]], sampleLabelsMerged[i])
      return(temp)
    })
    re <- NULL
    for(i in seq(sampleLabelsMerged)){re <- cbind(re, tss.new[[i]])}
    re <- cbind(tss.dt[,1:3],re)
    ##removes filtered rows
    re <- re[rowSums(re[,4:ncol(re)]) >0,]
    setorder(re, "strand","chr","pos")
    object@TSSprocessedMatrix <- re
    assign(objName, object, envir = parent.frame())
  }else{
    message("\tNo filtering method is defined...")
  }
})

################################################################################################
.filterWithPoisson <- function(data, coverageDepth, genomeSize, pVal) {
  # calculate lambda value (average)
  lambda <- coverageDepth / (genomeSize * 2)
  # get cutoff value
  cutoff <- qpois(pVal, lambda, lower.tail=FALSE, log.p=FALSE)
  # fiter tss table
  data[data<cutoff,] =0
  return(data)
}
