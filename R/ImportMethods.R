################################################################################################
#' Precisely identify TSSs from bam files, paired end bam files, bed files, BigWig files, tss files, or tss tables.
#'
#' @description getTSS function is used to precisely identify TSSs from multiple input file formats.
#' The files include users' home-made alignment files (bam format) or downloaded files from public databases.
#'  See inputFilesType for details on the supported input file formats.
#'
#' @usage 	getTSS(object, sequencingQualityThreshold = 10, mappingQualityThreshold = 20
#' , removeNewG = TRUE, correctG = TRUE)
#' @usage .getTSS_from_bam(bam.files, Genome, sampleLabels,inputFilesType
#' , sequencingQualityThreshold, mappingQualityThreshold)
#' @usage .getTSS_from_bed(bed.files, Genome, sampleLabels)
#' @usage .getTSS_from_BigWig(BigWig.files, Genome, sampleLabels)
#' @usage .getTSS_from_tss(tss.files, sampleLabels)
#' @usage .getTSS_from_TSStable(TSStable.file, sampleLabels)
#'
#' @param object A TSSr object.
#' @param sequencingQualityThreshold Used only if inputFilesType == "bam" or "bamPairedEnd", otherwise ignored.
#' @param mappingQualityThreshold Used only if inputFilesType == "bam" or "bamPairedEnd", otherwise ignored.
#'
#' @return
#' @export
#'
#' @importFrom Rsamtools ScanBamParam
#' @importFrom Rsamtools scanBamFlag
#' @importFrom Rsamtools scanBam
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomicRanges GRangesList
#' @importFrom GenomicRanges findOverlaps
#' @importFrom GenomicRanges elementMetadata
#' @importFrom GenomicAlignments seqnames
#' @importFrom GenomeInfoDb seqlengths
#'
#' @importFrom GenomicRanges resize
#' @importFrom BiocGenerics union
#'
#' @importFrom IRanges IPos
#' @importFrom IRanges IRanges
#' @importFrom IRanges RleList
#' @importFrom IRanges Views
#' @importFrom IRanges extractList
#' @importFrom IRanges reduce
#' @importFrom IRanges reverse
#' @importFrom IRanges viewApply
#'
#'
#' @importFrom stringr str_extract_all
#' @examples
#' getTSS(myTSSr)
setGeneric("getTSS",function(object
                             ,sequencingQualityThreshold = 10
                             ,mappingQualityThreshold = 20)standardGeneric("getTSS"))
#' @rdname getTSS
#' @export
setMethod("getTSS",signature(object = "TSSr"), function(object
                                    ,sequencingQualityThreshold = 10
                                    ,mappingQualityThreshold = 20
                                    ){
  ##initialize values
  Genome <- .getGenome(object@genomeName)
  sampleLabels <- object@sampleLabels
  inputFilesType <- object@inputFilesType
  if (length(object@sampleLabelsMerged) == 0) {
    object@sampleLabelsMerged <- sampleLabels
  }
  objName <- deparse(substitute(object))
  if(inputFilesType == "bam" | inputFilesType == "bamPairedEnd"){
    tss <- .getTSS_from_bam(object@inputFiles
                     ,Genome
                     ,sampleLabels
                     ,inputFilesType
                     ,sequencingQualityThreshold
                     ,mappingQualityThreshold
                     )
  }else if(inputFilesType == "bed"){
    tss <- .getTSS_from_bed(object@inputFiles, Genome, sampleLabels)
  }else if(inputFilesType == "BigWig"){
    tss <- .getTSS_from_BigWig(object@inputFiles,Genome, sampleLabels)
  }else if(inputFilesType == "tss"){
    tss <- .getTSS_from_tss(object@inputFiles, sampleLabels)
  }else if(inputFilesType == "TSStable"){
    tss <- .getTSS_from_TSStable(object@inputFiles, sampleLabels)
  }
  setorder(tss, "strand","chr","pos")
  # get library sizes
  object@librarySizes <- colSums(tss[,4:ncol(tss), drop = F], na.rm = T)

  object@TSSrawMatrix <- tss
  object@TSSprocessedMatrix <- tss
  assign(objName, object, envir = parent.frame())
})


