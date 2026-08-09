###############################################################################
##
.getConsensus <- function(gr1, cy, dis) {
    ## define variable as a NULL value
    dominant_tss <- NULL

    colnames(cy)[3:4] <- c("start.c", "end.c")
    cy[, start := dominant_tss - round(dis / 2)]
    cy[, end := dominant_tss + round(dis / 2)]
    gr2 <- makeGRangesFromDataFrame(cy, keep.extra.columns = FALSE)
    hits <- findOverlaps(gr1, gr2)
    new <- c(BiocGenerics::union(gr1[unname(as.data.frame(hits)[, 1])], gr2[unname(as.data.frame(hits)[, 2])]), gr1[-unname(as.data.frame(hits)[, 1])], gr2[-unname(as.data.frame(hits)[, 2])])
    return(new)
}

###############################################################################
##
.getConsensusQuantile <- function(tc, gr, tss.temp, useMultiCore, numCores) {
    interval_id <- point_id <- interquantile_width <- NULL
    `q_0.1` <- `q_0.9` <- NULL

    tc.hits <- .mapPointsToIntervals(
        tc,
        gr,
        pointColumn = "dominant_tss"
    )
    if (nrow(tc.hits) == 0L) {
        return(data.table())
    }

    spans <- tc.hits[, list(
        start = min(tc[["start"]][point_id]),
        end = max(tc[["end"]][point_id])
    ), by = interval_id]
    span.intervals <- data.table(
        chr = gr[["chr"]][spans[["interval_id"]]],
        strand = gr[["strand"]][spans[["interval_id"]]],
        start = spans[["start"]],
        end = spans[["end"]],
        consensus_id = spans[["interval_id"]]
    )
    tss.hits <- .mapPointsToIntervals(tss.temp, span.intervals)
    if (nrow(tss.hits) == 0L) {
        return(data.table())
    }
    matched <- tss.hits[, list(
        interval_id = span.intervals[["consensus_id"]][interval_id],
        point_id,
        pos = tss.temp[["pos"]][point_id],
        tags = tss.temp[["tags"]][point_id]
    )]

    summarize.one <- function(data) {
        reverse.order <- order(-data[["pos"]])
        q1.index <- .firstCumulativeFractionIndex(data[["tags"]], 0.1)
        q9.index <- .firstCumulativeFractionIndex(
            data[["tags"]][reverse.order], 0.1
        )
        list(
            start = min(data[["pos"]]),
            end = max(data[["pos"]]),
            dominant_tss = data[["pos"]][which.max(data[["tags"]])],
            tags = sum(data[["tags"]]),
            tags.dominant_tss = max(data[["tags"]]),
            q_0.1 = data[["pos"]][q1.index],
            q_0.9 = data[["pos"]][reverse.order[q9.index]]
        )
    }

    if (useMultiCore) {
        if (is.null(numCores)) {
            numCores <- detectCores()
        }
        message("process is running on ", numCores, " cores...")
        matched.groups <- split(
            matched,
            by = "interval_id",
            keep.by = TRUE
        )
        summaries <- mclapply(
            matched.groups,
            summarize.one,
            mc.cores = numCores
        )
        interval.ids <- vapply(
            matched.groups,
            function(data) data[["interval_id"]][[1L]],
            integer(1)
        )
        summary.table <- rbindlist(Map(
            function(interval.id, summary) {
                c(list(interval_id = interval.id), summary)
            },
            interval.ids,
            summaries
        ))
    } else {
        summary.table <- matched[, summarize.one(.SD), by = interval_id]
    }
    interval.ids <- summary.table[["interval_id"]]
    tc_clusters <- data.table(
        cluster = gr[["consensusCluster"]][interval.ids],
        chr = gr[["chr"]][interval.ids],
        start = as.integer(summary.table[["start"]]),
        end = as.integer(summary.table[["end"]]),
        strand = gr[["strand"]][interval.ids],
        dominant_tss = as.integer(summary.table[["dominant_tss"]]),
        tags = summary.table[["tags"]],
        tags.dominant_tss = summary.table[["tags.dominant_tss"]],
        q_0.1 = as.integer(summary.table[["q_0.1"]]),
        q_0.9 = as.integer(summary.table[["q_0.9"]])
    )
    tc_clusters[, interquantile_width := q_0.9 - q_0.1 + 1L]
    setorder(tc_clusters, "strand", "chr", "start")
    tc_clusters
}
