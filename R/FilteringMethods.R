setGeneric("filterTSS",function(object,...)standardGeneric("filterTSS"))
setMethod("filterTSS","TSSr", function(object
                                       ,Genome
                                       ,method = "poisson"
                                       ,Normalization = TRUE
                                       ,pVal=0.01
                                       ,tpmLow = 0.1
){
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
  ##filter tss data
  if(method  == "poisson"){
    message("\nFiltering data with ", method," method...")
    tss.raw <- object@TSSmergedMatrix
    tss.new <- lapply(as.list(seq(sampleLabelsMerged)), function(i){
      temp <- tss.raw[,.SD, .SDcols = sampleLabelsMerged[i]]
      setnames(temp, colnames(temp)[[1]], "tags")
      temp <- .filterWithPoisson(temp, library.size[i], genomeSize, pVal)
      if(Normalization == "TRUE"){
        sizePerMillion <- library.size[i] / 1e6
        temp[, tags := round(tags / sizePerMillion, 6)]
      }
      setnames(temp, colnames(temp)[[1]], sampleLabelsMerged[i])
      return(temp)
    })
    re <- NULL
    for(i in seq(sampleLabelsMerged)){re <- cbind(re, tss.new[[i]])}
    re <- cbind(tss.raw[,1:3],re)
    ##removes filtered rows
    re <- re[rowSums(re[,4:ncol(re)]) >0,]
    setorder(re, "strand","chr","pos")
    object@TSSfilteredMatrix <- re
    assign(objName, object, envir = parent.frame())
  }else if(method  == "TPM"){
    message("\nFiltering data with ", method," method...")
    tss.tpm <- object@TSSnormalizedMatrix
    tss.new <- lapply(as.list(seq(sampleLabelsMerged)), function(i){
      temp <- tss.tpm[,.SD, .SDcols = sampleLabelsMerged[i]]
      setnames(temp, colnames(temp)[[1]], "tags")
      temp[sampleLabelsMerged[i] < tpmLow,]=0
      setnames(temp, colnames(temp)[[1]], sampleLabelsMerged[i])
      return(temp)
    })
  re <- NULL
  for(i in seq(sampleLabelsMerged)){re <- cbind(re, tss.new[[i]])}
  re <- cbind(tss.tpm[,1:3],re)
  ##removes filtered rows
  re <- re[rowSums(re[,4:ncol(re)]) >0,]
  setorder(re, "strand","chr","pos")
  object@TSSfilteredMatrix <- re
  assign(objName, object, envir = parent.frame())
  }else{
    message("\tNo filtering method is defined...")
  }
})