###############################################################################

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

.clusterByPeak <- function(tss.dt, peakDistance, localThreshold, extensionDistance) {
    # create copy for reference later
    copied.dt <- copy(tss.dt)
    setkey(tss.dt, pos)
    ## define variable as a NULL value
    pos <- peak <- ID <- forward <- reverse <- V1 <- V2 <- chr <- NULL
    # get peakID
    pos_vec <- tss.dt[["pos"]]
    tags_vec <- tss.dt[["tags"]]
    left <- findInterval(pos_vec - peakDistance, pos_vec) + 1L
    right <- findInterval(pos_vec + peakDistance, pos_vec)
    equal_upper <- right > 0L & pos_vec[right] == pos_vec + peakDistance
    right[equal_upper] <- right[equal_upper] - 1L
    window_peak <- .rangeMaxFirstIndex(tags_vec)$query(left, right)
    peakID <- ifelse(pos_vec == pos_vec[window_peak], seq_along(pos_vec), 0L)

    # manipulate data.table to collapse clustered rows
    tss.dt[, peak := peakID]
    tss.dt[, ID := .I]
    ###############################################################################
    ## local filtering
    ###############################################################################
    if (unique(tss.dt$strand) == "+") {
        localF <- lapply(peakID[peakID > 0], function(i) {
            temp <- tss.dt[pos >= tss.dt$pos[i] & pos <= tss.dt$pos[i] + peakDistance, ]
            temp$ID[which(temp$tag < tss.dt$tags[i] * localThreshold)]
        })
    } else {
        localF <- lapply(peakID[peakID > 0], function(i) {
            temp <- tss.dt[pos >= tss.dt$pos[i] - peakDistance & pos <= tss.dt$pos[i], ]
            temp$ID[which(temp$tag < tss.dt$tags[i] * localThreshold)]
        })
    }
    if (length(unlist(localF)) > 0) {
        tss.dt <- tss.dt[-unlist(localF), ]
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
    tss_clusters <- lapply(as.list(seq_len(clusters[, .N])), function(i) {
        start <- clusters[i, V1]
        end <- clusters[i, V2]
        # copied.dt[, ID := .I]##NEW April18
        cluster.data <- copied.dt[pos >= start & pos <= end, ]
        tags.sum <- cluster.data[, sum(tags)] ## NEW Sep25, tags -> tags.sum
        q1 <- cluster.data[which(cumsum(tags) > 0.1 * tags.sum), min(pos)]
        q9 <- cluster.data[order(-pos)][which(cumsum(tags) > 0.1 * tags.sum), max(pos)]
        list(
            i,
            cluster.data[, chr[[1]]],
            start,
            end,
            cluster.data[, strand[[1]]]
            # ,cluster.data[which(ID %in% peakID), pos]  # NEW - use id column to find the intersection for the peakID vector (should hopefully only be 1!)
            , cluster.data[which.max(tags), pos],
            tags.sum,
            cluster.data[, max(tags)],
            q1,
            q9,
            q9 - q1 + 1
        )
    })

    # set names
    tss_clusters <- rbindlist(tss_clusters)
    setnames(tss_clusters, c(
        "cluster",
        "chr", "start", "end", "strand",
        "dominant_tss", "tags", "tags.dominant_tss",
        "q_0.1", "q_0.9", "interquantile_width"
    ))
    return(tss_clusters)
}
