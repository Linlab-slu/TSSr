#' TSSr example data
#'
#' A pre-processed TSSr object containing a subset of the CAGE dataset from
#' Lu and Lin, Genome Research 2019 Jul;29(7):1198-1210.
#'
#' This example object contains TSS data that has already been processed
#' through the TSSr workflow. The original BAM files are not included, so its
#' \code{inputFiles} and \code{refSource} slots are empty and \code{getTSS()}
#' cannot be called on this object. The \code{inputFilesType = "bam"} value is
#' retained only as provenance, and the bundled annotation is stored in
#' \code{refTable}.
#' However, all downstream analysis functions (such as \code{mergeSamples},
#' \code{normalizeTSS}, \code{filterTSS}, \code{clusterTSS}, etc.) can be
#' used with this object.
#'
#' @docType data
#'
#' @usage data(exampleTSSr)
#'
#' @format An S4 object of class \code{\linkS4class{TSSr}} containing:
#' \describe{
#'   \item{genomeName}{BSgenome reference name}
#'   \item{inputFiles}{An empty character vector because source BAM files are
#'     not bundled}
#'   \item{inputFilesType}{The original input type, retained as provenance}
#'   \item{sampleLabels}{Sample identifiers}
#'   \item{TSSrawMatrix}{Raw TSS count matrix}
#'   \item{TSSprocessedMatrix}{Processed TSS count matrix}
#'   \item{refSource}{An empty character vector; no machine-specific GFF path
#'     is stored}
#'   \item{refTable}{The bundled chromosome-I reference annotation}
#' }
#'
#' @source \url{https://genome.cshlp.org/content/29/7/1198.short}
#'
#' @examples
#' data(exampleTSSr)
#' exampleTSSr
#'
#' @keywords datasets
#' @name exampleTSSr
"exampleTSSr"
