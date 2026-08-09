###############################################################################
## .deseq2 function calcuates gene differential expression based on Deseq2 package
## .deseq2 function takes two assigned clusters and tss.raw table
## users need to provide which sample they want to compare and
## run script with the following example command:
## .deseq2(clustersx.asn,clustersy.asn, tss.raw,
##                              samplex <- c("ScerBY4741.1","ScerBY4741.2"),
##                              sampley <- c("ScerArrest.1","ScerArrest.2"),
##                              sampleOne <- "ScerBY4741",sampleTwo <- "ScerArrest")
############################################################################
## tss.raw is the raw tss merged tables, before any sums
############################################################################
.deseq2 <- function(cx, cy, tss.raw, samplex, sampley, sampleOne, sampleTwo,
                    useMultiCore, numCores, TAGtables) {
    .requireSuggestedPackage("DESeq2", "deGene()")
    ## get raw count tables
    if (sampleOne %in% names(TAGtables)) {
        xCounts <- TAGtables[[which(names(TAGtables) == sampleOne)]]
    } else {
        xCounts <- .tagCount_updated(cx, tss.raw, samplex, useMultiCore, numCores)
        ## save the tagCount results
        TAGtables[[sampleOne]] <- xCounts
    }

    if (sampleTwo %in% names(TAGtables)) {
        yCounts <- TAGtables[[which(names(TAGtables) == sampleTwo)]]
    } else {
        yCounts <- .tagCount_updated(cy, tss.raw, sampley, useMultiCore, numCores)
        TAGtables[[sampleTwo]] <- yCounts
    }
    xCounts <- xCounts[, -c(2:11)]
    yCounts <- yCounts[, -c(2:11)]
    ## tag counts by gene for sampleOne and sampleTwo
    setkey(xCounts, gene)
    setkey(yCounts, gene)
    one <- xCounts[, lapply(.SD, sum), by = gene, .SDcols = samplex]
    two <- yCounts[, lapply(.SD, sum), by = gene, .SDcols = sampley]
    ## merge the two raw count tables together by genes
    one[, (samplex) := lapply(.SD, as.integer), .SDcols = samplex]
    two[, (sampley) := lapply(.SD, as.integer), .SDcols = sampley]
    Dtable <- merge(
        as.data.frame(one), as.data.frame(two),
        by = "gene", all = TRUE
    )
    Dtable[is.na(Dtable)] <- 0
    ##
    gene.names <- Dtable[["gene"]]
    Dtable <- Dtable[, c(samplex, sampley), drop = FALSE]
    rownames(Dtable) <- gene.names
    Dtable <- data.matrix(Dtable)
    condition <- factor(c(rep(sampleOne, times = length(samplex)), rep(sampleTwo, times = length(sampley))))
    dds <- DESeq2::DESeqDataSetFromMatrix(
        countData = Dtable, data.frame(condition), ~condition
    )
    dds$condition <- factor(dds$condition, levels = c(sampleOne, sampleTwo))
    dds <- DESeq2::DESeq(dds)
    res <- DESeq2::results(dds)
    res <- res[order(res$padj), ]
    return(list(
        DEtable = as.data.frame(res),
        TAGtables = TAGtables
    ))
}

############################################################################
## .tagCount updated
.tagCount_updated <- function(cs, tss.raw, samples, useMultiCore, numCores) {
    ## define variable as a NULL value
    tag_sum <- NULL

    cols <- c("chr", "pos", "strand", samples)
    tss1 <- tss.raw[, .SD, .SDcols = cols]
    # exclude rows with no count
    tss1[, tag_sum := rowSums(tss1[, .SD, .SDcols = samples])]
    tss <- tss1[tag_sum != 0]
    tss[, tag_sum := NULL]

    tags <- matrix(
        0,
        nrow = nrow(cs),
        ncol = length(samples),
        dimnames = list(NULL, samples)
    )
    hits <- .mapPointsToIntervals(tss, cs)
    if (nrow(hits) > 0L) {
        matched <- data.table(
            interval_id = hits[["interval_id"]]
        )
        for (sample in samples) {
            matched[[sample]] <-
                tss[[sample]][hits[["point_id"]]]
        }
        counts <- matched[, lapply(.SD, sum),
            by = interval_id,
            .SDcols = samples
        ]
        tags[counts[["interval_id"]], ] <- as.matrix(
            counts[, .SD, .SDcols = samples]
        )
    }
    cbind(data.table::copy(cs), as.data.frame(tags))
}

## .tagCount is slow
.tagCount <- function(cs, tss.raw, samples, useMultiCore, numCores) {
    cols <- c("chr", "pos", "strand", samples)
    tss <- tss.raw[, .SD, .SDcols = cols]
    ## define variable as a NULL value
    chr <- strand <- start <- end <- NULL

    if (useMultiCore) {
        if (is.null(numCores)) {
            numCores <- detectCores()
        }
        message("process is running on ", numCores, " cores...")
        tags <- mclapply(seq_len(cs[, .N]), function(r) {
            data <- tss[tss$chr == cs[r, chr] & tss$strand == cs[r, strand] & tss$pos >= cs[r, start] & tss$pos <= cs[r, end], ]
            temp <- vapply(as.list(samples), function(s) {
                sum(data[, .SD, .SDcols = s])
            }, numeric(1))
            return(temp)
        }, mc.cores = numCores)
    } else {
        tags <- lapply(seq_len(cs[, .N]), function(r) {
            data <- tss[tss$chr == cs[r, chr] & tss$strand == cs[r, strand] & tss$pos >= cs[r, start] & tss$pos <= cs[r, end], ]
            temp <- vapply(as.list(samples), function(s) {
                sum(data[, .SD, .SDcols = s])
            }, numeric(1))
            return(temp)
        })
    }
    tags <- data.frame(matrix(unlist(tags), nrow = length(tags), byrow = TRUE), stringsAsFactors = FALSE)
    colnames(tags) <- samples
    cs <- cbind(cs, tags)
    return(cs)
}
