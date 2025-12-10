#' TSSr: A package for transcription start site sequencing data analyses.
#'
#' TSSr is designed to analyze transcription start sites (TSSs) and core promoters
#' with most types of 5’end sequencing data, such as cap analysis of gene expression (CAGE)
#' (Takahashi, Lassmann et al. 2012), no-amplification non-tagging CAGE libraries for
#' Illumina next-generation sequencers (nAnT-iCAGE) (Murata, Nishiyori-Sueki et al. 2014),
#' a Super-Low Input Carrier-CAGE (SLIC-CAGE) (Cvetesic, Leitch et al. 2018),
#' NanoCAGE (Cumbie, Ivanchenko et al. 2015), TSS-seq (Malabat, Feuerbach et al. 2015),
#' transcript isoform sequencing (TIF-seq) (Pelechano, Wei et al. 2013),
#' transcript-leaders sequencing (TL-seq) (Arribere and Gilbert 2013),
#' precision nuclear run-on sequencing (PRO-Cap) (Mahat, Kwak et al. 2016),
#' and GRO-Cap/5’GRO-seq (Kruesi, Core et al. 2013).
#'
#' TSSr package provides a comprehensive workflow on TSS data starts from identification
#' of accurate TSS locations, clustering TSSs within small genomic regions corresponding
#' to core promoters, and transcriptional activity quantifications, as well as specialized
#' downstream analyses including core promoter shape, cluster annotation, gene differential
#' expression, core promoter shift. TSSr can take multiple formats of files as input, such
#' as Binary Sequence Alignment Map (BAM) files (single-ended or paired-ended), Browser
#' Extension Data (bed) files, BigWig files, ctss files or tss tables. TSSr also generates
#' various types of TSS or core promoter track files which can be visualized in the UCSC
#' Genome Browser or Integrative Genomics Viewer (IGV). TSSr also exports downstream analyses
#' result tables and plots. Multiple cores are supported on Linux or Mac platforms.
#' @importFrom methods setClass setGeneric setMethod setRefClass
#' @import stringr
#' @import rtracklayer
#' @import ggplot2
#' @import GenomicFeatures
#' @import Gviz
#' @import DESeq2
#' @import calibrate
#' @import ggfortify
#' @import parallel
#' @importFrom IRanges IRanges
#' @importFrom IRanges findOverlaps
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom GenomicRanges elementMetadata
#' @importFrom GenomicRanges resize
#' @importFrom GenomicRanges mcols
#' @importFrom Rsamtools scanBam
#' @importFrom Rsamtools ScanBamParam
#' @importFrom Rsamtools bamFlag
#' @importFrom Rsamtools scanBamFlag
#' @importFrom GenomicAlignments qwidth
#' @importFrom GenomicAlignments readGAlignments
#' @importFrom GenomicAlignments seqnames
#' @importFrom dplyr filter
#' @import GenomeInfoDb
#' @importFrom BiocGenerics union
#' @import data.table
#' @import rmarkdown
#' @import parallel
#' @importFrom grDevices dev.off pdf rainbow
#' @importFrom graphics hist pairs par plot points strwidth text
#' @importFrom methods as new
#' @importFrom stats chisq.test cor p.adjust prcomp qpois
#' @importFrom utils installed.packages read.table write.table
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
#' @slot refTable data.frame.
#' @slot organismName character.
#' @slot assignedClusters list.
#' @slot unassignedClusters list.
#' @slot filteredClusters list.
#' @slot enhancers list.
#' @slot DEtables list.
#' @slot TAGtables list
#' @slot PromoterShift list.
#' @docType class
#' @name TSSr-class
#' @rdname TSSr-class
#' @exportClass TSSr
#'
setClass(Class = "TSSr",
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
                        ,refTable = "data.frame"
                        ,organismName = "character"
                        ,assignedClusters = "list"
                        ,unassignedClusters = "list"
                        ,filteredClusters = "list"
                        ,enhancers = "list"
                        ,DEtables = "list"
                        ,TAGtables = "list"
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
                   ,refTable = data.frame()
                   ,organismName = character()
                   ,assignedClusters = list()
                   ,unassignedClusters = list()
                   ,filteredClusters = list()
                   ,enhancers = list()
                   ,DEtables = list()
                   ,TAGtables = list()
                   ,PromoterShift = list()
         ),
         validity=function(object){
           supportedTypes <- c("bam", "bamPairedEnd", "bed", "tss", "TSStable", "BigWig")
           if(!(object@inputFilesType %in% supportedTypes))
             return(paste(sQuote("inputFilesType"), "must be one of supported input file types:",
                          paste(sQuote(supportedTypes), collapse = ", "), "."))
         })
