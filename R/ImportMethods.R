################################################################################################
##myTSSr class hold all data and metadata about one TSS experiment
setClass(Class = "TSSr",
         representation(genomeName = "character" 
                        ,inputFiles = "character"
                        ,inputFilesType = "character"
                        ,sampleLabels = "character"
                        ,sampleLabelsMerged = "character"
                        ,librarySizes = "integer"
                        ,TSSrawMatrix = "data.frame"
                        ,mergeIndex = "integer"
                        ,TSSmergedMatrix = "data.frame"
                        ,TSSnormalizedMatrix = "data.frame"
                        ,TSSfilteredMatrix = "data.frame"
                        ,tagClusters = "list"
                        ,consensusClusters = "list"
                        ,clusterShape = "list"
                        ,refSource = "character"
                        ,organismName = "character"
                        ,assignedClusters = "list"
                        ,unassignedClusters = "list"
                        ,DEtables = "list"
                        ,PromoterShift = "list"
                        ),
         prototype(genomeName = character()
                   ,inputFiles = character()
                   ,inputFilesType = character()
                   ,sampleLabels = character()
                   ,sampleLabelsMerged = character()
                   ,librarySizes = integer()
                   ,TSSrawMatrix = data.frame()
                   ,mergeIndex = integer()
                   ,TSSmergedMatrix = data.frame()
                   ,TSSnormalizedMatrix = data.frame()
                   ,TSSfilteredMatrix = data.frame()
                   ,tagClusters = list()
                   ,consensusClusters = list()
                   ,clusterShape = list()
                   ,refSource = character()
                   ,organismName = character()
                   ,assignedClusters = list()
                   ,unassignedClusters = list()
                   ,DEtables = list()
                   ,PromoterShift = list()
                   ),
         validity=function(object){
           supportedTypes <- c("bam", "bamPairedEnd", "bed", "ctss", "CTSStable", "BigWig")
           if(!(object@inputFilesType %in% supportedTypes))
             return(paste(sQuote("inputFilesType"), "must be one of supported input file types:",
                          paste(sQuote(supportedTypes), collapse = ", "), "."))
         })

################################################################################################
setGeneric("getTSS",function(object,...)standardGeneric("getTSS"))
setMethod("getTSS","TSSr", function(object
                                    ,sequencingQualityThreshold = 10
                                    ,mappingQualityThreshold = 20
                                    ){
  ##initialize values
  Genome <- .getGenome(object@genomeName)
  sampleLabels <- object@sampleLabels
  objName <- deparse(substitute(myTSSr))
  if(object@inputFilesType == "bam" | object@inputFilesType == "bamPairedEnd"){
    tss <- .getTSS_from_bam(object@inputFiles
                     ,Genome
                     ,sampleLabels
                     ,sequencingQualityThreshold
                     ,mappingQualityThreshold)
  }else if(object@inputFilesType == "bed"){
    tss <- .getTSS_from_bed(object@inputFiles, Genome, sampleLabels)
  }else if(object@inputFilesType == "BigWig"){
    tss <- .getTSS_from_BigWig(object@inputFiles,Genome, sampleLabels)
  }else if(object@inputFilesType == "ctss"){
    tss <- .getTSS_from_tss(object@inputFiles, sampleLabels)
  }else if(object@inputFilesType == "ctssTable"){
    tss <- .getTSS_from_TSStable(object@inputFiles)
  }
  setorder(tss, "strand","chr","pos")
  cat("\n")
  object@TSSrawMatrix <- tss
  assign(objName, object, envir = parent.frame())
})


