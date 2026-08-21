#' Select genes which have core promoter shift across different experiments.
#'
#' @description Selects genes which have multiple core promoters and undergo core promoter
#' shifting across different experiments. Generates gene list with Ds (degree of shift)
#' value (Lu et al., 2019), p value and adjusted p value.
#' Chi-squared tests use raw counts reconstructed from TPM values and the
#' corresponding library sizes. When sparse 2-by-2 tables have small expected
#' counts, the function emits one summary warning for all affected gene-level
#' tests. Their approximate p values may be unreliable; the warning should not
#' be treated as safe to ignore for inferential use.
#'
#' @usage shiftPromoter(object, comparePairs, pval=0.01)
#'
#' @param object A TSSr object.
#' @param comparePairs Specified list of sample pairs for comparison.
#' @param pval Genes with adjusted p value less than or equal to \code{pval}
#'   will be returned. Default value = 0.01.
#' @return A modified TSSr object with updated \code{PromoterShift} slot. The
#'   input object is not modified; assign the returned object to retain the
#'   changes.
#'
#' @export
#'
#' @examples
#' data(exampleTSSr)
#' exampleTSSr <- shiftPromoter(
#'     exampleTSSr, comparePairs = list(c("control", "treat")), pval = 0.01
#' )
#' head(PromoterShift(exampleTSSr, comparison = "control_VS_treat"))
setGeneric(
    "shiftPromoter",
    function(object, comparePairs, pval = 0.01) standardGeneric("shiftPromoter"),
    signature = "object"
)
#' @rdname shiftPromoter
#' @export
setMethod("shiftPromoter", signature(object = "TSSr"), function(
  object,
  comparePairs,
  pval
) {
    .requireWorkflowArtifact(
        object,
        "assignedClusters",
        "run annotateCluster() first"
    )
    .requireWorkflowArtifact(
        object,
        "librarySizes",
        "run getTSS() and mergeSamples() first"
    )
    ## initialize data
    message("\nCalculating core promoter shifts...")
    sampleLabelsMerged <- object@sampleLabelsMerged
    warningState <- new.env(parent = emptyenv())
    warningState$count <- 0L

    D <- lapply(as.list(seq(comparePairs)), function(i) {
        sampleOne <- comparePairs[[i]][1]
        sampleTwo <- comparePairs[[i]][2]
        cx <- object@assignedClusters[[sampleOne]]
        cy <- object@assignedClusters[[sampleTwo]]
        tss.raw <- object@TSSrawMatrix
        librarySizex <- object@librarySizes[which(sampleLabelsMerged == sampleOne)]
        librarySizey <- object@librarySizes[which(sampleLabelsMerged == sampleTwo)]
        DS <- withCallingHandlers(
            .Ds(
                cx,
                cy,
                librarySizex,
                librarySizey,
                useRawCount = TRUE,
                pval
            ),
            warning = function(warning_condition) {
                if (identical(
                    conditionMessage(warning_condition),
                    "Chi-squared approximation may be incorrect"
                )) {
                    warningState$count <- warningState$count + 1L
                    invokeRestart("muffleWarning")
                }
            }
        )
        return(DS)
    })
    if (warningState$count > 0L) {
        warning(
            "Chi-squared approximation may be incorrect for ",
            warningState$count,
            " gene-level test(s) with small expected counts; the affected ",
            "p-values may be unreliable.",
            call. = FALSE
        )
    }
    D.names <- vapply(as.list(seq(comparePairs)), function(i) {
        paste0(comparePairs[[i]][1], "_VS_", comparePairs[[i]][2], sep = "")
    }, character(1))
    names(D) <- D.names
    object@PromoterShift <- D
    return(object)
})
