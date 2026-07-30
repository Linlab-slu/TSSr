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
    vapply(strsplit(seq, ""), function(bases) {
        paste(rev(chartr("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN", bases)), collapse = "")
    }, character(1))
}

.uncodedGTrimWidth <- function(read.seq, reference.seq) {
    read.seq <- toupper(as.character(read.seq))
    reference.seq <- toupper(as.character(reference.seq))
    limit <- min(nchar(read.seq), nchar(reference.seq))
    trim.width <- 0L
    if (limit == 0L) {
        return(trim.width)
    }

    for (i in seq_len(limit)) {
        read.base <- substr(read.seq, i, i)
        reference.base <- substr(reference.seq, i, i)
        if (read.base == "G" && reference.base != "G") {
            trim.width <- trim.width + 1L
        } else {
            break
        }
    }
    trim.width
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

    trim.width[valid] <- vapply(seq_along(query.seq), function(i) {
        .uncodedGTrimWidth(query.seq[i], reference.seq[i])
    }, integer(1))

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
    message("\t-> Removing the bases of the reads if mismatched 'Gs'...")
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
    TSS <- TSS[, as.integer(sum(tag_count)), by = list(chr, pos, strand)]

    TSS
}

.getTSS_from_bam <- function(
  bam.files, Genome, sampleLabels, inputFilesType,
  sequencingQualityThreshold,
  mappingQualityThreshold,
  softclippingAllowed
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
    chunksize <- 1e6
    first <- TRUE
    for (i in seq_len(length(bam.files))) {
        message("\nReading in file: ", bam.files[i], "...")
        bam <- scanBam(bam.files[i], param = param)
        message("\t-> Filtering out low quality reads...")
        qual <- bam[[1]]$qual
        start <- 1
        # chunksize <- 1e6
        qa.avg <- vector(mode = "integer")
        repeat {
            if (start + chunksize <= length(qual)) {
                end <- start + chunksize
            } else {
                end <- length(qual)
            }
            qa.avg <- c(qa.avg, as.integer(vapply(as(qual[start:end], "IntegerList"), mean, numeric(1))))
            if (end == length(qual)) {
                break
            } else {
                start <- end + 1
            }
        }
        cigar <- bam[[1]]$cigar
        start <- 1
        # chunksize <- 1e6
        mapped.length <- vector(mode = "integer")
        repeat {
            if (start + chunksize <= length(cigar)) {
                end <- start + chunksize
            } else {
                end <- length(cigar)
            }
            mapped.length <- c(mapped.length, .cigarReferenceWidth(cigar[start:end]))
            if (end == length(cigar)) {
                break
            } else {
                start <- end + 1
            }
        }
        readsGR <- GRanges(
            seqnames = as.vector(bam[[1]]$rname), IRanges(start = bam[[1]]$pos, width = mapped.length),
            strand = bam[[1]]$strand, qual = qa.avg, mapq = bam[[1]]$mapq, seq = bam[[1]]$seq, read.length = width(bam[[1]]$seq),
            flag = bam[[1]]$flag
        )
        readsGR <- readsGR[as.character(seqnames(readsGR)) %in% seqnames(Genome)]
        readsGR <- readsGR[!(end(readsGR) > seqlengths(Genome)[as.character(seqnames(readsGR))])]
        GenomicRanges::elementMetadata(readsGR)$mapq[is.na(GenomicRanges::elementMetadata(readsGR)$mapq)] <- Inf
        readsGR.p <- readsGR[(as.character(strand(readsGR)) == "+" & GenomicRanges::elementMetadata(readsGR)$qual >=
            sequencingQualityThreshold) & GenomicRanges::elementMetadata(readsGR)$mapq >= mappingQualityThreshold]
        readsGR.m <- readsGR[(as.character(strand(readsGR)) == "-" & GenomicRanges::elementMetadata(readsGR)$qual >=
            sequencingQualityThreshold) & GenomicRanges::elementMetadata(readsGR)$mapq >= mappingQualityThreshold]
        if (softclippingAllowed) {
            TSS <- .makeTSSFromGRanges(readsGR.p, readsGR.m)
        } else {
            TSS <- .removeNewG(readsGR.p, readsGR.m, Genome)
        }

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
