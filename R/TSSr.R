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
#' as Binary Sequence Alignment Mao (BAM) files (single-ended or paired-ended), Browser
#' Extension Data (bed) files, BigWig files, ctss files or tss tables. TSSr also generates
#' various types of TSS or core promoter track files which can be visualized in the UCSC
#' Genome Browser or Integrative Genomics Viewer (IGV). TSSr also exports downstream analyses
#' result tables and plots. Multiple cores are supported on Linux or Mac platforms.
#'
#' @docType package
#' @name TSSr-package
NULL
#> NULL
#'
#' TSSr object
#'
#' @aliases TSSr
#' @exportClass TSSr
#'
