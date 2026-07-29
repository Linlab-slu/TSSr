.formatTSSrCount <- function(value) {
    format(value, big.mark = ",", scientific = FALSE, trim = TRUE)
}

.formatTSSrScalar <- function(value) {
    if (length(value) == 0L || is.na(value[[1L]]) || value[[1L]] == "") {
        return("<unset>")
    }
    as.character(value[[1L]])
}

.formatTSSrItems <- function(items, separator = ", ", maximum = 4L) {
    if (length(items) > maximum) {
        items <- c(
            items[seq_len(maximum)],
            sprintf("... +%d", length(items) - maximum)
        )
    }
    paste(items, collapse = separator)
}

.formatTSSrLabels <- function(labels) {
    if (length(labels) == 0L) {
        return("<none>")
    }
    .formatTSSrItems(labels)
}

.formatTSSrTableSummary <- function(tables) {
    if (length(tables) == 0L) {
        return("<not run>")
    }

    labels <- names(tables)
    if (is.null(labels) || any(labels == "")) {
        labels <- paste0("result", seq_along(tables))
    }
    row_counts <- vapply(tables, function(x) {
        if (is.data.frame(x) || is.matrix(x)) nrow(x) else NA_integer_
    }, integer(1))
    formatted_counts <- ifelse(
        is.na(row_counts),
        "unknown",
        vapply(row_counts, .formatTSSrCount, character(1))
    )

    .formatTSSrItems(paste(labels, formatted_counts), separator = "; ")
}

.formatTSSrComparisons <- function(tables) {
    if (length(tables) == 0L) {
        return("<not run>")
    }

    labels <- names(tables)
    if (is.null(labels) || any(labels == "")) {
        labels <- paste0("comparison", seq_along(tables))
    }
    sprintf("%d (%s)", length(tables), .formatTSSrItems(labels))
}

#' Display a compact summary of a TSSr object
#'
#' Shows the genome, sample and TSS counts, and concise row counts for completed
#' analyses without expanding the full contents of the object's slots.
#'
#' @param object A TSSr object.
#'
#' @return The input TSSr object, invisibly.
#'
#' @export
#'
#' @examples
#' data(exampleTSSr)
#' show(exampleTSSr)
setMethod(
    "show",
    signature(object = "TSSr"),
    function(object) {
        cat("TSSr object\n")
        cat(sprintf("  Genome: %s\n", .formatTSSrScalar(object@genomeName)))
        cat(sprintf(
            "  Samples: %d (%s)\n",
            length(object@sampleLabels),
            .formatTSSrLabels(object@sampleLabels)
        ))
        cat(sprintf(
            "  Merged samples: %d (%s)\n",
            length(object@sampleLabelsMerged),
            .formatTSSrLabels(object@sampleLabelsMerged)
        ))
        cat(sprintf(
            "  TSSs: raw %s; processed %s\n",
            .formatTSSrCount(nrow(object@TSSrawMatrix)),
            .formatTSSrCount(nrow(object@TSSprocessedMatrix))
        ))
        cat("  Analyses:\n")
        result_slots <- c(
            "tagClusters" = "Tag clusters",
            "consensusClusters" = "Consensus clusters",
            "clusterShape" = "Cluster shapes",
            "assignedClusters" = "Assigned clusters",
            "unassignedClusters" = "Unassigned clusters",
            "filteredClusters" = "Filtered clusters",
            "enhancers" = "Enhancers"
        )
        for (slot_name in names(result_slots)) {
            cat(sprintf(
                "    %s: %s\n",
                result_slots[[slot_name]],
                .formatTSSrTableSummary(methods::slot(object, slot_name))
            ))
        }
        cat(sprintf(
            "    DE comparisons: %s\n",
            .formatTSSrComparisons(object@DEtables)
        ))
        cat(sprintf(
            "    TAG tables: %s\n",
            .formatTSSrTableSummary(object@TAGtables)
        ))
        cat(sprintf(
            "    Promoter shifts: %s\n",
            .formatTSSrTableSummary(object@PromoterShift)
        ))
        invisible(object)
    }
)
