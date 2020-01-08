################################################################################################
#' Normalize raw TSS counts
#'
#' @description Normalizes raw TSS counts in all samples by tags per million (TPM)
#'
#' @usage normalizeTSS(object, Genome)
#'
#' @param object A TSSr object.
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#' normalizeTSS(myTSSr)
setGeneric("normalizeTSS",function(object,...)standardGeneric("normalizeTSS"))
#' @rdname normalizeTSS
#' @export
setMethod("normalizeTSS",signature(object = "TSSr"), function(object
                                       ,Genome
                                       ,Normalization = "TRUE"
){
  message("\nNormalizing TSS matrix...")
  ##initialize values
  Genome <- .getGenome(object@genomeName)
  sampleLabelsMerged <- object@sampleLabelsMerged
  objName <- deparse(substitute(object))

  tss.dt <- object@TSSprocessedMatrix
  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  if (all(is.wholenumber(object@librarySizes)) == FALSE) {
    stop('\tStopping... data is already normalized')
  }

  library.size <- object@librarySizes
  # if library size is empty, get library size
  # calculate size of genome
  genomeSize <- 0
  for (chrom in seq(Genome)) {
    genomeSize <- genomeSize + length(Genome[[chrom]])
  }
  ##normalize tss data
  if(Normalization == "TRUE"){
    tss.new <- lapply(as.list(seq(sampleLabelsMerged)), function(i){
      temp <- tss.dt[,.SD, .SDcols = sampleLabelsMerged[i]]
      sizePerMillion <- library.size[i] / 1e6
      setnames(temp, colnames(temp)[[1]], "tags")
      temp[, tags := round(tags / sizePerMillion, 6)]
      setnames(temp, colnames(temp)[[1]], sampleLabelsMerged[i])
      return(temp)
    })
  }else{
    message("\tNo normalization is defined...")
  }
  re <- NULL
  for(i in seq(sampleLabelsMerged)){re <- cbind(re, tss.new[[i]])}
  re <- cbind(tss.dt[,1:3],re)
  setorder(re, "strand","chr","pos")
  object@TSSprocessedMatrix <- re
  assign(objName, object, envir = parent.frame())
})
