###############################################################################

.firstCumulativeFractionIndex <- function(values, fraction) {
    cumulative <- cumsum(values)
    target <- fraction * sum(values)
    tolerance <- 8 * .Machine$double.eps * max(1, length(values)) *
        max(1, abs(target), abs(cumulative))

    which(cumulative >= target - tolerance)[[1L]]
}

.rangeMaxFirstIndex <- function(values) {
    n <- length(values)
    if (n == 0L) {
        return(list(query = function(left, right) integer(length(left))))
    }

    choose_best <- function(a, b) {
        replace_a <- values[b] > values[a] | (values[b] == values[a] & b < a)
        a[replace_a] <- b[replace_a]
        a
    }

    max_log <- floor(log2(n))
    table <- vector("list", max_log + 1L)
    table[[1L]] <- seq_len(n)
    if (max_log > 0L) {
        for (level in seq_len(max_log)) {
            span <- 2L^(level - 1L)
            prev <- table[[level]]
            len <- n - 2L^level + 1L
            left <- prev[seq_len(len)]
            right <- prev[seq_len(len) + span]
            table[[level + 1L]] <- choose_best(left, right)
        }
    }

    query <- function(left, right) {
        len <- right - left + 1L
        level <- floor(log2(len))
        left_idx <- integer(length(left))
        right_idx <- integer(length(left))
        for (current_level in unique(level)) {
            current <- which(level == current_level)
            span <- 2L^current_level
            current_table <- table[[current_level + 1L]]
            left_idx[current] <- current_table[left[current]]
            right_idx[current] <- current_table[right[current] - span + 1L]
        }
        choose_best(left_idx, right_idx)
    }

    list(query = query)
}

.selectPeakMaxIndices <- function(positions, tags, peakDistance) {
    n <- length(positions)
    if (n == 0L) {
        return(integer())
    }

    priority <- order(-tags, positions)
    left <- findInterval(positions - peakDistance, positions) + 1L
    right <- findInterval(
        positions + peakDistance,
        positions,
        left.open = TRUE
    )
    suppressed <- logical(n)
    selected <- integer(n)
    selected_count <- 0L

    for (index in priority) {
        if (suppressed[[index]]) {
            next
        }
        selected_count <- selected_count + 1L
        selected[[selected_count]] <- index
        suppressed[seq.int(left[[index]], right[[index]])] <- TRUE
    }

    sort(selected[seq_len(selected_count)])
}

.localFilterMask <- function(
    positions, tags, peakID, peakDistance, localThreshold, strand
) {
    n <- length(positions)
    if (n == 0L || !any(peakID > 0L)) {
        return(logical(n))
    }

    peak_tags <- rep.int(-Inf, n)
    peak_tags[peakID > 0L] <- tags[peakID > 0L]
    if (strand == "+") {
        left <- findInterval(
            positions - peakDistance,
            positions,
            left.open = TRUE
        ) + 1L
        right <- seq_len(n)
    } else {
        left <- seq_len(n)
        right <- findInterval(positions + peakDistance, positions)
    }
    strongest_peak <- .rangeMaxFirstIndex(peak_tags)$query(left, right)
    thresholds <- peak_tags[strongest_peak] * localThreshold
    remove <- is.finite(peak_tags[strongest_peak]) & tags < thresholds
    remove[is.na(remove)] <- FALSE
    remove
}

.summarizeClusterIntervals <- function(tss.dt, clusters) {
    positions <- tss.dt[["pos"]]
    tags <- tss.dt[["tags"]]
    starts <- clusters[["V1"]]
    ends <- clusters[["V2"]]
    left <- findInterval(starts, positions, left.open = TRUE) + 1L
    right <- findInterval(ends, positions)
    dominant <- .rangeMaxFirstIndex(tags)$query(left, right)

    summaries <- lapply(seq_along(left), function(index) {
        rows <- seq.int(left[[index]], right[[index]])
        cluster_tags <- tags[rows]
        q1_relative <- .firstCumulativeFractionIndex(cluster_tags, 0.1)
        q9_reverse <- .firstCumulativeFractionIndex(
            rev(cluster_tags), 0.1
        )
        q1 <- positions[rows[[q1_relative]]]
        q9 <- positions[rows[[length(rows) - q9_reverse + 1L]]]

        list(
            index,
            tss.dt[["chr"]][[rows[[1L]]]],
            starts[[index]],
            ends[[index]],
            tss.dt[["strand"]][[rows[[1L]]]],
            positions[dominant[[index]]],
            sum(cluster_tags),
            tags[dominant[[index]]],
            q1,
            q9,
            q9 - q1 + 1
        )
    })
    rbindlist(summaries)
}

.clusterByPeak <- function(
    tss.dt, peakDistance, localThreshold, extensionDistance
) {
    # The public clusterTSS() entry point requires peakDistance > 0.
    tss.dt <- copy(tss.dt)
    setkeyv(tss.dt, "pos")
    pos_vec <- tss.dt[["pos"]]
    tags_vec <- tss.dt[["tags"]]
    left <- findInterval(pos_vec - peakDistance, pos_vec) + 1L
    right <- findInterval(pos_vec + peakDistance, pos_vec)
    equal_upper <- right > 0L & pos_vec[right] == pos_vec + peakDistance
    right[equal_upper] <- right[equal_upper] - 1L
    window_peak <- .rangeMaxFirstIndex(tags_vec)$query(left, right)
    peakID <- ifelse(pos_vec == pos_vec[window_peak], seq_along(pos_vec), 0L)

    .clusterFromPeakIDs(
        tss.dt, peakID, peakDistance, localThreshold, extensionDistance
    )
}

.clusterByPeakMax <- function(
    tss.dt, peakDistance, localThreshold, extensionDistance
) {
    # The public clusterTSS() entry point requires peakDistance > 0.
    tss.dt <- copy(tss.dt)
    setkeyv(tss.dt, "pos")
    selected <- .selectPeakMaxIndices(
        tss.dt[["pos"]], tss.dt[["tags"]], peakDistance
    )
    peakID <- integer(nrow(tss.dt))
    peakID[selected] <- selected

    .clusterFromPeakIDs(
        tss.dt, peakID, peakDistance, localThreshold, extensionDistance
    )
}

.clusterFromPeakIDs <- function(
    tss.dt, peakID, peakDistance, localThreshold, extensionDistance
) {
    copied.dt <- copy(tss.dt)
    tss.dt <- copy(tss.dt)
    setkey(tss.dt, pos)
    ## define variable as a NULL value
    pos <- peak <- ID <- forward <- reverse <- V1 <- V2 <- chr <- NULL

    # manipulate data.table to collapse clustered rows
    tss.dt[, peak := peakID]
    tss.dt[, ID := .I]
    ###############################################################################
    ## local filtering
    ###############################################################################
    local_filter <- .localFilterMask(
        tss.dt[["pos"]], tss.dt[["tags"]], peakID,
        peakDistance, localThreshold, unique(tss.dt$strand)
    )
    if (any(local_filter)) {
        tss.dt <- tss.dt[!local_filter, ]
    }
    ###############################################################################
    ###############################################################################
    tss.dt[, forward := ifelse(data.table::shift(pos, 1, type = "lead") < pos + extensionDistance, 1, 0)] #
    tss.dt[, reverse := ifelse(data.table::shift(pos, 1, type = "lag") > pos - extensionDistance, 1, 0)]
    tss.dt <- tss.dt[, list(peak = max(peak), start = min(pos), end = max(pos), tags = sum(tags)), by = list(rleid(peak, forward, reverse))] ## ZL?

    # get start and end boundaries for clusters
    # TODO revisit this code for better optimization
    clusters <- lapply(as.list(tss.dt[peak > 0, rleid]), function(x) {
        start <- tss.dt[x, start]
        end <- tss.dt[x, end]

        if (x - 1 > 0 && tss.dt[x - 1, !peak > 0] && tss.dt[x - 1, end] > start - extensionDistance) {
            start <- tss.dt[x - 1, start]
            if (x - 2 > 0 && tss.dt[x - 2, !peak > 0] && tss.dt[x - 2, end] > start - extensionDistance) {
                start <- tss.dt[x - 2, start]
            }
        }
        if (x + 1 < tss.dt[, .N] && tss.dt[x + 1, !peak > 0] && tss.dt[x + 1, start] < end + extensionDistance) {
            end <- tss.dt[x + 1, end]
            if (x + 2 < tss.dt[, .N] && tss.dt[x + 2, !peak > 0] && tss.dt[x + 2, start] < end + extensionDistance) {
                end <- tss.dt[x + 2, end]
            }
        }
        list(start, end)
    })

    clusters <- rbindlist(clusters)

    # deal with overlapping clusters here
    # TODO this section needs some more optimization/work

    rowVec <- which(clusters$V2 >= data.table::shift(clusters$V1, 1, type = "lead"))
    if (length(rowVec) > 0) {
        ###############################################################################
        ###############################################################################
        for (i in seq_len(length(rowVec))) {
            clusters$V1[rowVec[i] + 1] <- clusters$V1[rowVec[i]]
        }
        clusters <- clusters[-rowVec, ]
    } ##


    #  clusters <- unique(clusters)

    # get full clustering data
    # core promoter boundaries are calculated here (i.e. cumsum distribution)
    tss_clusters <- .summarizeClusterIntervals(copied.dt, clusters)
    setnames(tss_clusters, c(
        "cluster",
        "chr", "start", "end", "strand",
        "dominant_tss", "tags", "tags.dominant_tss",
        "q_0.1", "q_0.9", "interquantile_width"
    ))
    return(tss_clusters)
}
