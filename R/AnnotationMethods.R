###############################################################################
#' Annotate clusters with GFF annotation file.
#'
#' @description Annotates clusters with gene or transcript names from GFF annotation file.
#'
#' @usage    annotateCluster(object,clusters = "consensusClusters",filterCluster = TRUE
#' , filterClusterThreshold = 0.02, annotationType = "genes",upstream=1000
#' , upstreamOverlap = 500,downstream = 0)
#'
#' @param object A TSSr object containing the selected cluster tables and
#'   either a populated \code{refTable} slot or an existing GFF file in
#'   \code{refSource}.
#' @param clusters Clusters to be annotated: "consensusClusters" or "tagClusters".
#' Default is "consensusClusters".
#' @param filterCluster Logical indicating whether clusters downstream of a highly
#' expressed cluster are filtered. Setting filterCluster as "TRUE" would reduce weak
#' clusters brought from recapping, transcriptional or sequencing noise. Default is TRUE.
#' @param filterClusterThreshold  Ignore downstream clusters if signal < filterClusterThreshold*the
#' strongest clusters within the same gene promoter region. Default value = 0.02.
#' @param annotationType  Specify annotation feature to be associated with: "genes"
#' or "transcripts". Default is "genes".
#' @param upstream  Upstream distance to the start position of annotation feature.
#' Default value = 1000.
#' @param upstreamOverlap Upstream distance to the start position of annotation
#' feature if overlapped with the upstream neighboring feature. Default value = 500.
#' @param downstream  Downstream distance to the start position of annotation feature.
#' Default value = 0. Note: if annotationType == "transctipt" or the gene annotations
#' start from transcription start sites (TSSs), the recommended value = 500.
#' @return A modified TSSr object with updated \code{assignedClusters},
#'   \code{unassignedClusters}, and optionally \code{filteredClusters} slots.
#'   The input object is not modified; assign the returned object to retain
#'   the changes.
#'
#' @export
#'
#' @examples
#' data(exampleTSSr)
#' exampleTSSr <- annotateCluster(
#'     exampleTSSr, clusters = "consensusClusters", filterCluster = TRUE
#' )
#' head(assignedClusters(exampleTSSr, sample = "control"))
#'
setGeneric("annotateCluster", function(
  object, clusters = "consensusClusters",
  filterCluster = TRUE,
  filterClusterThreshold = 0.02,
  annotationType = "genes",
  upstream = 1000,
  upstreamOverlap = 500,
  downstream = 0
) standardGeneric("annotateCluster"), signature = "object")
#' @rdname annotateCluster
#' @export
setMethod("annotateCluster", signature(object = "TSSr"), function(
  object, clusters, filterCluster,
  filterClusterThreshold, annotationType,
  upstream, upstreamOverlap, downstream
) {
    if (length(clusters) != 1L ||
        !clusters %in% c("tagClusters", "consensusClusters")) {
        stop(
            "'clusters' must be either 'tagClusters' or 'consensusClusters'.",
            call. = FALSE
        )
    }
    if (length(annotationType) != 1L ||
        !annotationType %in% c("genes", "transcripts")) {
        stop(
            "'annotationType' must be either 'genes' or 'transcripts'.",
            call. = FALSE
        )
    }
    if (length(filterCluster) != 1L || is.na(filterCluster) ||
        !is.logical(filterCluster)) {
        stop("'filterCluster' must be TRUE or FALSE.", call. = FALSE)
    }
    next_steps <- if (identical(clusters, "tagClusters")) {
        "run clusterTSS() first"
    } else {
        "run clusterTSS() and consensusCluster() first"
    }
    cs.dt <- .requireWorkflowArtifact(object, clusters, next_steps)

    message("\nAnnotating...")
    sampleLabelsMerged <- object@sampleLabelsMerged
    if (length(sampleLabelsMerged) == 0L) {
        stop(
            "'sampleLabelsMerged' is empty; define sample groups before annotation.",
            call. = FALSE
        )
    }
    missing_samples <- setdiff(sampleLabelsMerged, names(cs.dt))
    if (length(missing_samples) > 0L) {
        stop(
            "'", clusters, "' has no table for merged sample(s): ",
            paste(sQuote(missing_samples), collapse = ", "), ".",
            call. = FALSE
        )
    }
    refGFF <- object@refSource
    refTable <- object@refTable

    ## check whether there is refTable provided
    if (length(refTable) != 0) {
        ref <- refTable
    } else {
        ## check whether there is refSource provided
        if (length(refGFF) == 0) {
            stop("Please provide correct refSource file!")
        }
        .validateFilePaths(
            refGFF,
            "refSource",
            requireSingle = TRUE
        )
        ## define variable as a NULL value
        inCoding <- r <- f <- dominant_tss <- gene_id <- tx_name <- tx_id <- NULL
        ## prepare annotation file
        txdb <- withCallingHandlers(
            makeTxDbFromGFF(refGFF, format = "auto"),
            warning = function(warning_condition) {
                if (identical(
                    conditionMessage(warning_condition),
                    paste(
                        "genome version information is not available for this",
                        "TxDb object"
                    )
                )) {
                    invokeRestart("muffleWarning")
                }
            }
        )
        if (annotationType == "genes") {
            ref <- setDT(as.data.frame(genes(txdb)))
        } else if (annotationType == "transcripts") {
            ref <- setDT(as.data.frame(transcripts(txdb)))
            ref[, gene_id := ifelse(
                is.na(tx_name) | !nzchar(as.character(tx_name)),
                paste0("transcript_", tx_id),
                as.character(tx_name)
            )]
        }
        object@refTable <- ref
    }
    ref <- data.table::as.data.table(data.table::copy(ref))
    required_reference_columns <- c(
        "seqnames", "start", "end", "width", "strand", "gene_id"
    )
    missing_reference_columns <- setdiff(required_reference_columns, names(ref))
    if (length(missing_reference_columns) > 0L) {
        stop(
            "Reference annotation is missing required column(s): ",
            paste(sQuote(missing_reference_columns), collapse = ", "), ".",
            call. = FALSE
        )
    }

    required_cluster_columns <- c(
        "cluster", "chr", "start", "end", "strand", "dominant_tss"
    )
    if (isTRUE(filterCluster)) {
        required_cluster_columns <- c(
            required_cluster_columns,
            "tags", "tags.dominant_tss", "q_0.1", "q_0.9",
            "interquantile_width"
        )
    }
    for (sample_label in sampleLabelsMerged) {
        missing_cluster_columns <- setdiff(
            required_cluster_columns,
            names(cs.dt[[sample_label]])
        )
        if (length(missing_cluster_columns) > 0L) {
            stop(
                "Cluster table for ", sQuote(sample_label),
                " is missing required column(s): ",
                paste(sQuote(missing_cluster_columns), collapse = ", "), ".",
                call. = FALSE
            )
        }
    }
    cluster_chromosomes <- unique(unlist(lapply(
        cs.dt[sampleLabelsMerged],
        function(cluster_table) as.character(cluster_table$chr)
    )))
    if (length(intersect(cluster_chromosomes, as.character(ref$seqnames))) == 0L) {
        stop(
            "No chromosome names are shared between clusters and the reference ",
            "annotation; check chromosome naming conventions.",
            call. = FALSE
        )
    }

    ##
    asn <- lapply(as.list(seq(sampleLabelsMerged)), function(i) {
        cs.temp <- data.table::as.data.table(
            data.table::copy(cs.dt[[sampleLabelsMerged[i]]])
        )
        cs.asn <- .assign2gene(cs.temp, ref, upstream, upstreamOverlap, downstream, filterCluster)
        return(cs.asn)
    })
    names(asn) <- sampleLabelsMerged
    ## subset assigned
    asned <- lapply(as.list(seq(sampleLabelsMerged)), function(i) {
        cs <- asn[[sampleLabelsMerged[i]]]
        cs <- cs[!is.na(cs$gene), ]
        setDT(cs)
        return(cs)
    })
    ## subset unassigned
    unasn <- lapply(as.list(seq(sampleLabelsMerged)), function(i) {
        cs <- asn[[sampleLabelsMerged[i]]]
        cs <- cs[is.na(cs$gene), ]
        setDT(cs)
        return(cs)
    })
    ## filter clusters
    if (isTRUE(filterCluster)) {
        filteredColumns <- c(
            "cluster", "chr", "start", "end", "strand", "dominant_tss",
            "tags", "tags.dominant_tss", "q_0.1", "q_0.9",
            "interquantile_width", "gene"
        )
        asn.filtered <- lapply(as.list(seq(sampleLabelsMerged)), function(i) {
            cs <- asn[[sampleLabelsMerged[i]]]
            m <- cs[is.na(gene) & is.na(inCoding), ]
            n <- cs[!is.na(gene) | !is.na(inCoding), ]
            n[, inCoding := ifelse(!is.na(inCoding) & !is.na(gene), NA, inCoding)]
            n[, r := ifelse(is.na(gene), inCoding, gene)]
            setkey(n, r)
            if (nrow(n) == 0L) {
                new <- n
            } else {
                new <- lapply(as.list(n[, unique(r)]), function(i) {
                    temp <- n[list(i)]
                    max.tags <- temp[, max(tags)]
                    if (temp[1, strand] == "+") {
                        temp[, f := ifelse(dominant_tss > temp[which.max(tags), dominant_tss] & tags < max.tags * filterClusterThreshold, 0, 1)]
                    } else {
                        temp[, f := ifelse(dominant_tss < temp[which.max(tags), dominant_tss] & tags < max.tags * filterClusterThreshold, 0, 1)]
                    }
                    return(temp)
                })
                new <- rbindlist(new)
                new <- new[f == 1, ]
            }
            cs <- rbindlist(list(
                m[, .SD, .SDcols = filteredColumns],
                new[, .SD, .SDcols = filteredColumns]
            ), use.names = TRUE)
            return(cs)
        })
    }
    names(asned) <- sampleLabelsMerged
    names(unasn) <- sampleLabelsMerged
    object@assignedClusters <- asned
    object@unassignedClusters <- unasn
    if (filterCluster == TRUE) {
        names(asn.filtered) <- sampleLabelsMerged
        object@filteredClusters <- asn.filtered
    }
    return(object)
})
