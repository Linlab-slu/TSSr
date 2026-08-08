################################################################################
.asIntegerCoordinate <- function(x, column = "pos", minimum = 1L) {
    if (!is.numeric(x)) {
        stop("Coordinate column '", column, "' must be numeric.")
    }
    if (any(!is.finite(x))) {
        stop("Coordinate column '", column, "' contains non-finite values.")
    }
    if (any(x != trunc(x))) {
        stop("Coordinate column '", column, "' must contain whole numbers.")
    }
    if (any(x < minimum | x > .Machine$integer.max)) {
        stop(
            "Coordinate column '", column, "' must be between ",
            minimum, " and ", .Machine$integer.max, "."
        )
    }
    as.integer(x)
}

################################################################################
.getGenome <- function(genomeName) {
    if (is.null(genomeName)) {
        stop("Can not run this function with a NULL genome.")
    }
    if (genomeName %in% rownames(installed.packages()) == FALSE) {
        stop("Requested genome is not installed! Please install required BSgenome package before running TSSr.")
    }
    requireNamespace(genomeName)
    getExportedValue(genomeName, genomeName)
}

################################################################################
## .getTSS_from_bam function calls TSS from bam files
.cigarReferenceWidth <- function(cigar) {
    as.integer(cigarillo::cigar_extent_along_ref(cigar, N.regions.removed = FALSE))
}

.reverseComplementText <- function(seq) {
    seq <- toupper(as.character(seq))
    complemented <- chartr(
        "ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN", seq
    )
    as.character(Biostrings::reverse(Biostrings::BStringSet(complemented)))
}

.uncodedGTrimWidths <- function(read.seq, reference.seq) {
    read.seq <- toupper(as.character(read.seq))
    reference.seq <- toupper(as.character(reference.seq))
    if (length(read.seq) != length(reference.seq)) {
        stop("read.seq and reference.seq must have equal lengths.")
    }

    trim.width <- integer(length(read.seq))
    read.width <- nchar(read.seq)
    reference.width <- nchar(reference.seq)
    active <- which(read.width > 0L & reference.width > 0L)
    position <- 1L

    while (length(active) > 0L) {
        active <- active[
            read.width[active] >= position &
                reference.width[active] >= position
        ]
        if (length(active) == 0L) {
            break
        }

        active <- active[
            substr(read.seq[active], position, position) == "G" &
                substr(reference.seq[active], position, position) != "G"
        ]
        if (length(active) == 0L) {
            break
        }

        trim.width[active] <- trim.width[active] + 1L
        position <- position + 1L
    }
    trim.width
}

.uncodedGTrimWidth <- function(read.seq, reference.seq) {
    .uncodedGTrimWidths(read.seq, reference.seq)[[1L]]
}

.trimUncodedGOneStrand <- function(readsGR, Genome, minusStrand = FALSE) {
    if (length(readsGR) == 0L) {
        return(readsGR)
    }

    read.seq <- as.character(GenomicRanges::elementMetadata(readsGR)$seq)
    terminal.width <- pmin(
        nchar(read.seq),
        pmax(width(readsGR), 0L)
    )
    valid <- which(terminal.width > 0L)
    trim.width <- integer(length(readsGR))
    if (length(valid) == 0L) {
        return(readsGR)
    }

    terminal.ranges <- GenomicRanges::resize(
        readsGR[valid],
        width = terminal.width[valid],
        fix = "start"
    )
    reference.seq <- getSeq(Genome, terminal.ranges, as.character = TRUE)

    if (minusStrand) {
        seq.end <- nchar(read.seq[valid])
        query.seq <- .reverseComplementText(substr(
            read.seq[valid],
            seq.end - terminal.width[valid] + 1L,
            seq.end
        ))
    } else {
        query.seq <- substr(read.seq[valid], 1L, terminal.width[valid])
    }

    trim.width[valid] <- .uncodedGTrimWidths(query.seq, reference.seq)

    trimmed <- which(trim.width > 0L)
    if (length(trimmed) > 0L) {
        if (minusStrand) {
            end(readsGR)[trimmed] <- end(readsGR)[trimmed] - trim.width[trimmed]
        } else {
            start(readsGR)[trimmed] <- start(readsGR)[trimmed] + trim.width[trimmed]
        }
    }
    readsGR
}

.removeNewG <- function(readsGR.p, readsGR.m, Genome) {
    readsGR.p <- .trimUncodedGOneStrand(readsGR.p, Genome, minusStrand = FALSE)
    readsGR.m <- .trimUncodedGOneStrand(readsGR.m, Genome, minusStrand = TRUE)

    .makeTSSFromGRanges(readsGR.p, readsGR.m)
}

.makeTSSFromGRanges <- function(readsGR.p, readsGR.m) {
    chr <- pos <- tag_count <- NULL

    TSS.p <- data.table(
        chr = as.character(seqnames(readsGR.p)),
        pos = as.integer(start(readsGR.p)), strand = "+"
    )
    TSS.m <- data.table(
        chr = as.character(seqnames(readsGR.m)),
        pos = as.integer(end(readsGR.m)), strand = "-"
    )
    TSS <- rbind(TSS.p, TSS.m)
    TSS <- TSS[, c("chr", "pos", "strand")]
    TSS$tag_count <- 1
    setDT(TSS)
    TSS <- TSS[, list(
        tag_count = as.integer(sum(tag_count))
    ), by = list(chr, pos, strand)]

    TSS
}

.aggregateTSSCountTables <- function(tables) {
    chr <- pos <- strand <- tag_count <- NULL

    if (length(tables) == 0L) {
        return(.emptyTSSCounts())
    }

    counts <- rbindlist(tables, use.names = TRUE)
    counts[, list(
        tag_count = as.integer(sum(tag_count))
    ), by = list(chr, pos, strand)]
}

.combineSampleTSSTables <- function(tables, sampleLabels) {
    key.columns <- c("chr", "pos", "strand")
    if (length(tables) == 0L || length(tables) != length(sampleLabels)) {
        stop(
            "The number of sample TSS tables must match the number of ",
            "sample labels."
        )
    }

    sample.info <- lapply(seq_along(tables), function(i) {
        table <- data.table::copy(data.table::as.data.table(tables[[i]]))
        expected.names <- c(key.columns, sampleLabels[[i]])
        if (!identical(names(table), expected.names)) {
            stop(
                "The TSS table for sample '", sampleLabels[[i]],
                "' must contain columns: ",
                paste(expected.names, collapse = ", "), "."
            )
        }
        data.table::set(
            table, j = "pos",
            value = .asIntegerCoordinate(table[["pos"]])
        )
        if (anyDuplicated(table, by = key.columns)) {
            stop(
                "Duplicate chr/pos/strand coordinates in the TSS table ",
                "for sample '", sampleLabels[[i]], "'."
            )
        }
        if (!is.numeric(table[[sampleLabels[[i]]]])) {
            stop(
                "TSS counts for sample '", sampleLabels[[i]],
                "' must be numeric."
            )
        }
        missing.count <- anyNA(table[[sampleLabels[[i]]]])
        data.table::setnafill(
            table, type = "const", fill = 0,
            cols = sampleLabels[[i]]
        )
        list(table = table, missing.count = missing.count)
    })
    sample.tables <- lapply(sample.info, function(info) info$table)
    missing.counts <- vapply(
        sample.info, function(info) info$missing.count, logical(1)
    )

    result <- unique(data.table::rbindlist(lapply(
        sample.tables,
        function(table) table[, key.columns, with = FALSE]
    )))
    data.table::setorderv(result, c("strand", "chr", "pos"))

    for (i in seq_along(sample.tables)) {
        sample.label <- sampleLabels[[i]]
        values <- sample.tables[[i]]
        count <- values[[sample.label]]
        needs.missing.fill <- missing.counts[[i]] ||
            nrow(values) < nrow(result)
        if (is.integer(count) && !needs.missing.fill) {
            result[, (sample.label) := rep(0L, .N)]
        } else if (is.numeric(count)) {
            result[, (sample.label) := rep(0, .N)]
        } else {
            stop("TSS counts for sample '", sample.label,
                 "' must be numeric.")
        }
        data.table::setnames(values, sample.label, "value")
        if (nrow(values) > 0L) {
            result[
                values,
                on = c("chr", "pos", "strand"),
                (sample.label) := get("i.value")
            ]
        }
    }
    data.table::setcolorder(result, c(key.columns, sampleLabels))
    result
}

.appendTSSChunkTable <- function(
  tables, chunk, compactionLimit = 8L
) {
    tables[[length(tables) + 1L]] <- chunk
    if (length(tables) >= compactionLimit) {
        return(list(.aggregateTSSCountTables(tables)))
    }
    tables
}

.getTSS_from_bam <- function(
  bam.files, Genome, sampleLabels, inputFilesType,
  sequencingQualityThreshold,
  mappingQualityThreshold,
  softclippingAllowed,
  bamYieldSize
) {
    ## define variable as a NULL value
    chr <- pos <- tag_count <- strand <- NULL

    what <- c("rname", "strand", "pos", "seq", "qual", "mapq", "flag", "cigar")
    param <- ScanBamParam(
        what = what,
        flag = scanBamFlag(
            isUnmappedQuery = FALSE,
            isNotPassingQualityControls = FALSE
        ),
        mapqFilter = mappingQualityThreshold
    )
    if (inputFilesType == "bamPairedEnd") {
        Rsamtools::bamFlag(param) <- scanBamFlag(
            isUnmappedQuery = FALSE,
            isProperPair    = TRUE,
            isFirstMateRead = TRUE
        )
    }
    first <- TRUE
    for (i in seq_len(length(bam.files))) {
        message("\nReading in file: ", bam.files[i], "...")
        message("\t-> Filtering out low quality reads...")
        if (!softclippingAllowed) {
            message("\t-> Removing the bases of the reads if mismatched 'Gs'...")
        }
        bamFile <- Rsamtools::BamFile(
            bam.files[i], yieldSize = bamYieldSize
        )
        open(bamFile)
        chunk.tables <- tryCatch({
            chunks <- list()
            repeat {
                bam <- scanBam(bamFile, param = param)
                if (length(bam[[1L]]$pos) == 0L) {
                    break
                }
                chunk <- .bamChunkToTSS(
                    bam,
                    Genome,
                    sequencingQualityThreshold,
                    mappingQualityThreshold,
                    softclippingAllowed
                )
                if (nrow(chunk) > 0L) {
                    chunks <- .appendTSSChunkTable(chunks, chunk)
                }
            }
            chunks
        }, finally = {
            close(bamFile)
        })

        TSS <- .aggregateTSSCountTables(chunk.tables)

        setnames(TSS, c("chr", "pos", "strand", sampleLabels[i]))
        setkey(TSS, chr, pos, strand)
        if (first == TRUE) {
            TSS.all.samples <- TSS
        } else {
            TSS.all.samples <- merge(TSS.all.samples, TSS, all = TRUE)
        }
        first <- FALSE
        gc()
    }
    TSS.all.samples[, 4:ncol(TSS.all.samples)][is.na(TSS.all.samples[, 4:ncol(TSS.all.samples)])] <- 0
    TSS.all.samples[, pos := .asIntegerCoordinate(pos)]
    return(TSS.all.samples)
}

.emptyTSSCounts <- function() {
    data.table(
        chr = character(),
        pos = integer(),
        strand = character(),
        tag_count = integer()
    )
}

.bamChunkToTSS <- function(
  bam, Genome, sequencingQualityThreshold, mappingQualityThreshold,
  softclippingAllowed
) {
    qa.avg <- as.integer(vapply(
        as(bam[[1L]]$qual, "IntegerList"), mean, numeric(1)
    ))
    mapped.length <- .cigarReferenceWidth(bam[[1L]]$cigar)
    readsGR <- GRanges(
        seqnames = as.vector(bam[[1L]]$rname),
        IRanges(start = bam[[1L]]$pos, width = mapped.length),
        strand = bam[[1L]]$strand,
        qual = qa.avg,
        mapq = bam[[1L]]$mapq,
        seq = bam[[1L]]$seq,
        read.length = width(bam[[1L]]$seq),
        flag = bam[[1L]]$flag
    )
    readsGR <- readsGR[as.character(seqnames(readsGR)) %in% seqnames(Genome)]
    readsGR <- readsGR[
        !(end(readsGR) >
            seqlengths(Genome)[as.character(seqnames(readsGR))])
    ]
    GenomicRanges::elementMetadata(readsGR)$mapq[
        is.na(GenomicRanges::elementMetadata(readsGR)$mapq)
    ] <- Inf
    readsGR.p <- readsGR[
        as.character(strand(readsGR)) == "+" &
            GenomicRanges::elementMetadata(readsGR)$qual >=
                sequencingQualityThreshold &
            GenomicRanges::elementMetadata(readsGR)$mapq >=
                mappingQualityThreshold
    ]
    readsGR.m <- readsGR[
        as.character(strand(readsGR)) == "-" &
            GenomicRanges::elementMetadata(readsGR)$qual >=
                sequencingQualityThreshold &
            GenomicRanges::elementMetadata(readsGR)$mapq >=
                mappingQualityThreshold
    ]
    if (softclippingAllowed) {
        .makeTSSFromGRanges(readsGR.p, readsGR.m)
    } else {
        .removeNewG(readsGR.p, readsGR.m, Genome)
    }
}

################################################################################
## .getTSS_from_bed function calls TSS from bed files
.getTSS_from_bed <- function(bed.files, Genome, sampleLabels) {
    first <- TRUE
    ## define variable as a NULL value
    chr <- pos <- tag_count <- NULL

    for (i in seq_len(length(bed.files))) {
        message("\nReading in file: ", bed.files[i], "...")
        readsGR <- import(bed.files[i], format = "BED")
        readsGR <- readsGR[as.character(seqnames(readsGR)) %in% seqnames(Genome)]
        readsGR <- readsGR[!(end(readsGR) > seqlengths(Genome)[as.character(seqnames(readsGR))])]
        readsGR.p <- readsGR[strand(readsGR) == "+"]
        readsGR.m <- readsGR[strand(readsGR) == "-"]
        message("\t-> Making TSS table...")
        TSS.plus <- data.table(chr = as.character(seqnames(readsGR.p)), pos = as.integer(start(readsGR.p)), strand = rep("+", times = length(readsGR.p)))
        TSS.minus <- data.table(chr = as.character(seqnames(readsGR.m)), pos = as.integer(end(readsGR.m)), strand = rep("-", times = length(readsGR.m)))
        TSS <- rbind(TSS.plus, TSS.minus)
        TSS$tag_count <- 1
        TSS <- TSS[, as.integer(sum(tag_count)), by = list(chr, pos, strand)]
        setnames(TSS, c("chr", "pos", "strand", sampleLabels[i]))
        setkey(TSS, chr, pos, strand)
        if (first == TRUE) {
            TSS.all.samples <- TSS
        } else {
            TSS.all.samples <- merge(TSS.all.samples, TSS, all = TRUE)
        }
        first <- FALSE
        gc()
    }
    TSS.all.samples[, 4:ncol(TSS.all.samples)][is.na(TSS.all.samples[, 4:ncol(TSS.all.samples)])] <- 0
    TSS.all.samples[, pos := .asIntegerCoordinate(pos)]
    return(TSS.all.samples)
}
################################################################################
## .getTSS_from_BigWig function calls TSS from BigWig files

.getTSS_from_BigWig <- function(BigWig.files, Genome, sampleLabels) {
    # library.sizes <- vector()
    ## define variable as a NULL value
    chr <- pos <- NULL
    first <- TRUE
    for (i in seq_len(length(BigWig.files))) {
        message("\nReading in file: ", BigWig.files[i], "...")
        readsGR <- import(BigWig.files[i], format = "BigWig")
        readsGR <- readsGR[as.character(seqnames(readsGR)) %in% seqnames(Genome)]
        readsGR <- readsGR[!(end(readsGR) > seqlengths(Genome)[as.character(seqnames(readsGR))])]
        readsGR.p <- readsGR[score(readsGR) > 0]
        readsGR.m <- readsGR[score(readsGR) < 0]
        message("\t-> Making TSS table...")
        TSS.plus <- data.table(chr = as.character(seqnames(readsGR.p)), pos = as.integer(start(readsGR.p)), strand = rep("+", times = length(readsGR.p)), score = as.numeric(abs(readsGR.p$score)))
        TSS.minus <- data.table(chr = as.character(seqnames(readsGR.m)), pos = as.integer(end(readsGR.m)), strand = rep("-", times = length(readsGR.m)), score = as.numeric(abs(readsGR.m$score)))
        TSS <- rbind(TSS.plus, TSS.minus)

        setDT(TSS)

        setnames(TSS, c("chr", "pos", "strand", sampleLabels[i]))
        setkey(TSS, chr, pos, strand)

        # library.sizes <- c(library.sizes, as.integer(sum(data.table(TSS)[,4])))
        if (first == TRUE) {
            TSS.all.samples <- TSS
        } else {
            TSS.all.samples <- merge(TSS.all.samples, TSS, all = TRUE)
        }
        first <- FALSE
        gc()
    }
    TSS.all.samples[, 4:ncol(TSS.all.samples)][is.na(TSS.all.samples[, 4:ncol(TSS.all.samples)])] <- 0
    TSS.all.samples[, pos := .asIntegerCoordinate(pos)]
    return(TSS.all.samples)
}


################################################################################
## .getTSS_from_tss function calls TSS from tss files

.getTSS_from_tss <- function(tss.files, sampleLabels) {
    pos <- NULL
    first <- TRUE

    for (i in seq_len(length(tss.files))) {
        message("\nReading in file: ", tss.files[i], "...")
        TSS <- read.table(
            file = tss.files[i], header = TRUE, sep = "\t",
            colClasses = c("character", "numeric", "character", "numeric"),
            col.names = c("chr", "pos", "strand", sampleLabels[i])
        )

        setDT(TSS)
        TSS[, pos := .asIntegerCoordinate(pos)]

        setkeyv(TSS, cols = c("chr", "pos", "strand"))
        if (first == TRUE) {
            TSS.all.samples <- TSS
        } else {
            TSS.all.samples <- merge(TSS.all.samples, TSS, all = TRUE)
        }
        first <- FALSE
        gc()
    }
    TSS.all.samples <- data.table(TSS.all.samples)
    TSS.all.samples[, 4:ncol(TSS.all.samples)][is.na(TSS.all.samples[, 4:ncol(TSS.all.samples)])] <- 0
    return(TSS.all.samples)
}

################################################################################################
## .getTSS_from_TSStable function calls TSS from one TSStable file

.getTSS_from_TSStable <- function(TSStable.file, sampleLabels) {
    pos <- NULL
    if (length(TSStable.file) > 1) {
        stop("Only one file should be provided when inputFilesType = \"TSStable\"!")
    }
    if (file.exists(TSStable.file) == FALSE) {
        stop("Could not locate input file ", TSStable.file)
    }

    TSS.all.samples <- read.table(
        file = TSStable.file, header = TRUE, stringsAsFactors = FALSE,
        colClasses = c(
            "character", "numeric", "character",
            rep("numeric", length(sampleLabels))
        ),
        col.names = c("chr", "pos", "strand", sampleLabels)
    )
    if (ncol(TSS.all.samples) != (length(sampleLabels) + 3)) {
        stop("Number of provided sample labels must match the number of samples in the TSS table!")
    }
    setDT(TSS.all.samples)
    TSS.all.samples[, pos := .asIntegerCoordinate(pos)]
    TSS.all.samples[, 4:ncol(TSS.all.samples)][is.na(TSS.all.samples[, 4:ncol(TSS.all.samples)])] <- 0
    return(TSS.all.samples)
}
