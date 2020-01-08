################################################################################################
##
#' @name TSSr-class
#' @docType class
#' @noRd
#' @export
#'
#' @import stringr
#' @import data.table
#' @import rtracklayer
#' @import ggplot2
#' @importFrom GenomicAlignments qwidth
#' @importFrom GenomicAlignments readGAlignments
#' @import GenomicFeatures
#' @import Gviz
#' @import DESeq2
#' @import BSgenome.Scerevisiae.UCSC.sacCer3
#' @import calibrate
#' @import ggfortify
#'
#'
#' @slot genomeName character.
#' @slot inputFiles character.
#' @slot inputFilesType character.
#' @slot sampleLabels character.
#' @slot sampleLabelsMerged character.
#' @slot librarySizes numeric.
#' @slot TSSrawMatrix data.frame.
#' @slot mergeIndex numeric.
#' @slot TSSprocessedMatrix data.frame.
#' @slot tagClusters list.
#' @slot consensusClusters list.
#' @slot clusterShape list.
#' @slot refSource character.
#' @slot organismName character.
#' @slot assignedClusters list.
#' @slot unassignedClusters list.
#' @slot filteredClusters list.
#' @slot DEtables list.
#' @slot PromoterShift list.
#'
#'
TSSr <- setClass(Class = "TSSr",
         representation(genomeName = "character"
                        ,inputFiles = "character"
                        ,inputFilesType = "character"
                        ,sampleLabels = "character"
                        ,sampleLabelsMerged = "character"
                        ,librarySizes = "numeric"
                        ,TSSrawMatrix = "data.frame"
                        ,mergeIndex = "numeric"
                        ,TSSprocessedMatrix = "data.frame"
                        ,tagClusters = "list"
                        ,consensusClusters = "list"
                        ,clusterShape = "list"
                        ,refSource = "character"
                        ,organismName = "character"
                        ,assignedClusters = "list"
                        ,unassignedClusters = "list"
                        ,filteredClusters = "list"
                        ,DEtables = "list"
                        ,PromoterShift = "list"
                        ),
         prototype(genomeName = character()
                   ,inputFiles = character()
                   ,inputFilesType = character()
                   ,sampleLabels = character()
                   ,sampleLabelsMerged = character()
                   ,librarySizes = numeric()
                   ,TSSrawMatrix = data.frame()
                   ,mergeIndex = numeric()
                   ,TSSprocessedMatrix = data.frame()
                   ,tagClusters = list()
                   ,consensusClusters = list()
                   ,clusterShape = list()
                   ,refSource = character()
                   ,organismName = character()
                   ,assignedClusters = list()
                   ,unassignedClusters = list()
                   ,filteredClusters = list()
                   ,DEtables = list()
                   ,PromoterShift = list()
                   ),
         validity=function(object){
           supportedTypes <- c("bam", "bamPairedEnd", "bed", "tss", "TSStable", "BigWig")
           if(!(object@inputFilesType %in% supportedTypes))
             return(paste(sQuote("inputFilesType"), "must be one of supported input file types:",
                          paste(sQuote(supportedTypes), collapse = ", "), "."))
         })

