#' Analysis of core promoter shape
#'
#' @description Calculates core promoter shape based on the distributions of TSSs within core
#'  promoters using Shape Index (SI) algorithm (Hoskins et al. 2011) or Promoter Shape Score (PSS)
#'   algorithm (Lu et al. 2019).
#'
#' @usage shapeCluster(object, clusters = "consensusClusters", method = "PSS",
#'  useMultiCore=FALSE, numCores = NULL)
#'
#' @param object A TSSr object.
#' @param clusters Clusters to be used for calculating shape score: "tagClusters" or "consensusClusters".
#'  Default is "consensusClusters".
#' @param method Method to be used for calculating core promoter shape score: "SI" or "PSS". Default is "PSS".
#' @param useMultiCore Logical indicating whether multiple cores are used (TRUE) or not (FALSE). Default is FALSE.
#' @param numCores Number of cores are used in clustering step. Used only if useMultiCore = TRUE. Default is NULL.
#' @return A modified TSSr object with updated \code{clusterShape} slot. The
#'   input object is not modified; assign the returned object to retain the
#'   changes.
#'
#' @export
#'
#' @examples
#' example_input <- system.file(
#'     "extdata", "example-tss-table.tsv", package = "TSSr",
#'     mustWork = TRUE
#' )
#' example_object <- TSSr(
#'     genomeName = "BSgenome.Scerevisiae.UCSC.sacCer3",
#'     inputFiles = example_input,
#'     inputFilesType = "TSStable",
#'     sampleLabels = c("SL01", "SL02", "SL03", "SL04"),
#'     sampleLabelsMerged = c("control", "treat"),
#'     mergeIndex = c(1, 1, 2, 2)
#' )
#' example_object <- getTSS(example_object)
#' example_object <- mergeSamples(example_object)
#' example_object <- normalizeTSS(example_object)
#' example_object <- clusterTSS(example_object, useMultiCore = FALSE)
#' example_object <- consensusCluster(
#'     example_object, useMultiCore = FALSE
#' )
#' example_object <- shapeCluster(
#'     example_object, clusters = "consensusClusters", method = "PSS",
#'     useMultiCore = FALSE
#' )
#' head(clusterShape(example_object, sample = "control"))
setGeneric("shapeCluster", function(
  object, clusters = "consensusClusters",
  method = "PSS", useMultiCore = FALSE, numCores = NULL
) standardGeneric("shapeCluster"), signature = "object")
#' @rdname shapeCluster
#' @export
setMethod("shapeCluster", signature(object = "TSSr"), function(object, clusters, method, useMultiCore, numCores) {
    message("\nCalculating ", clusters, " shape with ", method, " method...")
    ## define variable as a NULL value
    pos <- interquantile_width <- chr <- NULL

    ## initialize data
    tss.dt <- object@TSSprocessedMatrix

    if (clusters == "tagClusters") {
        cs.dt <- object@tagClusters
    } else if (clusters == "consensusClusters") {
        cs.dt <- object@consensusClusters
    }

    sampleLabelsMerged <- object@sampleLabelsMerged

    calculate.sample.shape <- function(i) {
        cs <- data.table::copy(cs.dt[[sampleLabelsMerged[[i]]]])
        if (nrow(cs) == 0L) {
            return(data.table())
        }
        tss <- tss.dt[, .SD, .SDcols = c(
            "chr", "pos", "strand", sampleLabelsMerged[[i]]
        )]
        setnames(tss, colnames(tss)[[4L]], "tags")
        tss <- tss[tags > 0, ]

        if (!method %in% c("PSS", "SI")) {
            message("\nNo shape method is provided...")
            return(cs)
        }

        hits <- .mapPointsToIntervals(
            tss,
            cs,
            intervalStart = "q_0.1",
            intervalEnd = "q_0.9"
        )
        default.score <- if (method == "PSS") 0 else 2
        score <- rep(default.score, nrow(cs))
        if (nrow(hits) > 0L) {
            matched <- hits[, .(
                interval_id,
                tags = tss[["tags"]][point_id]
            )]
            contributions <- matched[, {
                total <- sum(tags)
                proportions <- tags / total
                .(contribution = sum(
                    proportions * log(proportions, 2)
                ))
            }, by = interval_id]
            if (method == "PSS") {
                score[contributions[["interval_id"]]] <-
                    -contributions[["contribution"]] *
                    log(
                        cs[["interquantile_width"]][
                            contributions[["interval_id"]]
                        ],
                        2
                    )
            } else {
                score[contributions[["interval_id"]]] <-
                    2 + contributions[["contribution"]]
            }
        }
        cs[, shape.score := score]
        setkey(cs, NULL)
        cs
    }

    sample.indices <- seq_along(sampleLabelsMerged)
    if (useMultiCore) {
        if (is.null(numCores)) {
            numCores <- detectCores()
        }
        message("process is running on ", numCores, " cores...")
        cs.shape <- mclapply(
            sample.indices,
            calculate.sample.shape,
            mc.cores = numCores
        )
    } else {
        cs.shape <- lapply(sample.indices, calculate.sample.shape)
    }
    names(cs.shape) <- sampleLabelsMerged
    object@clusterShape <- cs.shape
    return(object)
})
