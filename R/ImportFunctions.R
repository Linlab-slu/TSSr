################################################################################################
.getGenome <- function(genomeName) {
  if (is.null(genomeName)){
    stop("Can not run this function with a NULL genome.")}
  if(genomeName %in% rownames(installed.packages()) == FALSE){
    stop("Requested genome is not installed! Please install required BSgenome package before running TSSr.")}
  requireNamespace(genomeName)
  getExportedValue(genomeName, genomeName)
}

################################################################################################
##.getTSS_from_bam function calls TSS from bam files
.getTSS_from_bam <- function(bam.files, Genome, sampleLabels, inputFilesType
                             ,sequencingQualityThreshold
                             ,mappingQualityThreshold){

  what <- c("rname", "strand", "pos", "seq", "qual", "mapq","flag","cigar")
  param <- ScanBamParam( what = what
                         , flag = scanBamFlag(isUnmappedQuery = FALSE,
                                              isNotPassingQualityControls = FALSE)
                         , mapqFilter = mappingQualityThreshold)
  if (inputFilesType == "bamPairedEnd"){
    bamFlag(param) <- scanBamFlag( isUnmappedQuery = FALSE
                                   , isProperPair    = TRUE
                                   , isFirstMateRead = TRUE)}
  first <- TRUE
  for(i in 1:length(bam.files)) {
    message("\nReading in file: ", bam.files[i], "...")
    bam <- scanBam(bam.files[i], param = param)
    message("\t-> Filtering out low quality reads...")
    qual <- bam[[1]]$qual
    start <- 1
    chunksize <- 1e6
    qa.avg <- vector(mode = "integer")
    repeat {
      if (start + chunksize <= length(qual)) {
        end <- start + chunksize
      } else {
        end <- length(qual)
      }
      qa.avg <- c(qa.avg, as.integer(sapply(as(qual[start:end], "IntegerList"),mean)))
      if (end == length(qual)) {
        break
      } else {
        start <- end + 1
      }
    }
    cigar <- bam[[1]]$cigar
    start <- 1
    chunksize <- 1e6
    mapped.length <- vector(mode = "integer")
    repeat {
      if (start + chunksize <= length(cigar)) {
        end <- start + chunksize
      } else {
        end <- length(cigar)
      }
      mapped.length <- c(mapped.length, as.integer(sum(as(str_extract_all(bam[[1]]$cigar[start:end], "([0-9]+)"),"IntegerList"))))
      if (end == length(cigar)) {
        break
      } else {
        start <- end + 1
      }
    }
    readsGR <- GRanges(seqnames = as.vector(bam[[1]]$rname), IRanges(start = bam[[1]]$pos, width = mapped.length),
                       strand = bam[[1]]$strand, qual = qa.avg, mapq = bam[[1]]$mapq, seq = bam[[1]]$seq, read.length = width(bam[[1]]$seq),
                       flag = bam[[1]]$flag)
    readsGR <- readsGR[unique(as.character(readsGR@seqnames)) %in% seqnames(Genome)]
    readsGR <- readsGR[!(end(readsGR) > seqlengths(Genome)[as.character(seqnames(readsGR))])]
    GenomicRanges::elementMetadata(readsGR)$mapq[is.na(GenomicRanges::elementMetadata(readsGR)$mapq)] <- Inf
    readsGR.p <- readsGR[(as.character(strand(readsGR)) == "+" & GenomicRanges::elementMetadata(readsGR)$qual >= sequencingQualityThreshold) & GenomicRanges::elementMetadata(readsGR)$mapq >= mappingQualityThreshold]
    readsGR.m <- readsGR[(as.character(strand(readsGR)) == "-" & GenomicRanges::elementMetadata(readsGR)$qual >= sequencingQualityThreshold) & GenomicRanges::elementMetadata(readsGR)$mapq >= mappingQualityThreshold]
    # remove G mismatch
    TSS <- .removeNewG(readsGR.p, readsGR.m, Genome)

    setnames(TSS, c("chr", "pos", "strand", sampleLabels[i]))
    setkey(TSS, chr, pos, strand)
    if(first == TRUE) {
      TSS.all.samples <- TSS
    }else{
      TSS.all.samples <- merge(TSS.all.samples, TSS, all = TRUE)
    }
    first <- FALSE
  }
  TSS.all.samples[,4:ncol(TSS.all.samples)][is.na(TSS.all.samples[,4:ncol(TSS.all.samples)])] =0
  return(TSS.all.samples)
}

################################################################################################
.removeNewG <- function(readsGR.p, readsGR.m, Genome) {

  message("\t-> Removing the first base of the reads if mismatched 'G'...")
  G.reads.p <- which(substr(GenomicRanges::elementMetadata(readsGR.p)$seq, start = 1, stop = 1) == "G")
  G.reads.m <- which(substr(GenomicRanges::elementMetadata(readsGR.m)$seq, start = GenomicRanges::elementMetadata(readsGR.m)$read.length, stop = GenomicRanges::elementMetadata(readsGR.m)$read.length) == "C")
  if(length(G.reads.p)>0){
    G.mismatch.p <- G.reads.p[getSeq(Genome, GenomicRanges::resize(readsGR.p[G.reads.p], width = 1, fix = "start"), as.character = TRUE) != "G"]
    GenomicRanges::elementMetadata(readsGR.p)$removedG <- FALSE
    GenomicRanges::elementMetadata(readsGR.p)$removedG[G.mismatch.p] <- TRUE
    start(readsGR.p)[G.mismatch.p] <- start(readsGR.p)[G.mismatch.p] + as.integer(1)
    TSS.p <- data.table(chr = as.character(seqnames(readsGR.p)), pos = start(readsGR.p), strand = "+", removedG = GenomicRanges::elementMetadata(readsGR.p)$removedG, stringsAsFactors = FALSE)
  }else{
    G.mismatch.p <- NULL
    TSS.p <- data.table()
  }
  if(length(G.reads.m)>0){
    G.mismatch.m <- G.reads.m[getSeq(Genome, GenomicRanges::resize(readsGR.m[G.reads.m], width = 1, fix = "start"), as.character = TRUE) != "G"]
    GenomicRanges::elementMetadata(readsGR.m)$removedG <- FALSE
    GenomicRanges::elementMetadata(readsGR.m)$removedG[G.mismatch.m] <- TRUE
    end(readsGR.m)[G.mismatch.m] <- end(readsGR.m)[G.mismatch.m] - as.integer(1)
    TSS.m <- data.table(chr = as.character(seqnames(readsGR.m)), pos = end(readsGR.m), strand = "-", removedG = GenomicRanges::elementMetadata(readsGR.m)$removedG, stringsAsFactors = FALSE)
  }else{
    G.mismatch.m <- NULL
    TSS.m <- data.table()
  }
  TSS <- rbind(TSS.p, TSS.m)
  TSS <- TSS[,c("chr", "pos", "strand")]
  TSS$tag_count <- 1
  setDT(TSS)
  TSS <- TSS[, as.integer(sum(tag_count)), by = list(chr, pos, strand)]

  return(TSS)
}

################################################################################################
##.getTSS_from_bed function calls TSS from bed files
.getTSS_from_bed <- function(bed.files, Genome, sampleLabels){
  first <- TRUE
  for(i in 1:length(inputFiles)) {
    message("\nReading in file: ", inputFiles[i], "...")
    readsGR <- import(inputFiles[i], format = "BED")
    readsGR <- readsGR[seqnames(readsGR) %in% seqnames(Genome)]
    readsGR <- readsGR[!(end(readsGR) > seqlengths(Genome)[as.character(seqnames(readsGR))])]
    readsGR.p <- readsGR[strand(readsGR) == "+"]
    readsGR.m <- readsGR[strand(readsGR) == "-"]
    message("\t-> Making TSS table...")
    TSS.plus <- data.table(chr = as.character(seqnames(readsGR.p)), pos = as.integer(start(readsGR.p)), strand = rep("+", times = length(readsGR.p)), stringsAsFactors = F)
    TSS.minus <- data.table(chr = as.character(seqnames(readsGR.m)), pos = as.integer(end(readsGR.m)), strand = rep("-", times = length(readsGR.m)), stringsAsFactors = F)
    TSS <- rbind(TSS.plus, TSS.minus)
    TSS$tag_count <- 1
    TSS <- data.table(TSS)
    TSS <- TSS[, as.integer(sum(tag_count)), by = list(chr, pos, strand)]
    setnames(TSS, c("chr", "pos", "strand", sampleLabels[i]))
    setkey(TSS, chr, pos, strand)
    if(first == TRUE) {
      TSS.all.samples <- TSS
    }else{
      TSS.all.samples <- merge(TSS.all.samples, TSS, all = TRUE)
    }
    first <- FALSE
  }
  TSS.all.samples[,4:ncol(TSS.all.samples)][is.na(TSS.all.samples[,4:ncol(TSS.all.samples)])] =0
  return(TSS.all.samples)
}
################################################################################################
##.getTSS_from_BigWig function calls TSS from BigWig files

.getTSS_from_BigWig <- function(BigWig.files, Genome, sampleLabels){
  #library.sizes <- vector()
  first <- TRUE
  for(i in 1:length(inputFiles)) {
    message("\nReading in file: ", inputFiles[i], "...")
    readsGR <- import(inputFiles[i], format = "BigWig")
    readsGR <- readsGR[seqnames(readsGR) %in% seqnames(Genome)]
    readsGR <- readsGR[!(end(readsGR) > seqlengths(Genome)[as.character(seqnames(readsGR))])]
    readsGR.p <- readsGR[score(readsGR) > 0]
    readsGR.m <- readsGR[score(readsGR) < 0]
    message("\t-> Making TSS table...")
    TSS.plus <- data.table(chr = as.character(seqnames(readsGR.p)), pos = as.integer(start(readsGR.p)), strand = rep("+", times = length(readsGR.p)), score = as.numeric(abs(readsGR.p$score)), stringsAsFactors = F)
    TSS.minus <- data.table(chr = as.character(seqnames(readsGR.m)), pos = as.integer(end(readsGR.m)), strand = rep("-", times = length(readsGR.m)), score = as.numeric(abs(readsGR.m$score)), stringsAsFactors = F)
    TSS <- rbind(TSS.plus, TSS.minus)

    setDT(TSS)

    setnames(TSS, c("chr", "pos", "strand", sampleLabels[i]))
    setkey(TSS, chr, pos, strand)

    #library.sizes <- c(library.sizes, as.integer(sum(data.table(TSS)[,4])))
    if(first == TRUE) {
      TSS.all.samples <- TSS
    }else{
      TSS.all.samples <- merge(TSS.all.samples, TSS, all = TRUE)
    }
    first <- FALSE
  }
  TSS.all.samples[,4:ncol(TSS.all.samples)][is.na(TSS.all.samples[,4:ncol(TSS.all.samples)])] =0
  return(TSS.all.samples)
}


################################################################################################
##.getTSS_from_tss function calls TSS from tss files

.getTSS_from_tss <- function(tss.files, sampleLabels){
  first <- TRUE
  for(i in 1:length(inputFiles)) {
    message("\nReading in file: ", inputFiles[i], "...")
    TSS <- read.table(file = inputFiles[i], header = F, sep = "\t"
                      ,colClasses = c("character", "integer", "character", "integer")
                      ,col.names = c("chr", "pos", "strand", sampleLabels[i]))

    setDT(TSS)

    setkeyv(TSS, cols = c("chr", "pos", "strand"))
    if(first == TRUE) {
      TSS.all.samples <- TSS
    }else{
      TSS.all.samples <- merge(TSS.all.samples, TSS, all = TRUE)
    }
    first <- FALSE
  }
  TSS.all.samples <- data.table(TSS.all.samples)
  TSS.all.samples[,4:ncol(TSS.all.samples)][is.na(TSS.all.samples[,4:ncol(TSS.all.samples)])] =0
  return(TSS.all.samples)
}

################################################################################################
##.getTSS_from_TSStable function calls TSS from one TSStable file

.getTSS_from_TSStable <- function(TSStable.file, sampleLabels){
  if(length(inputFiles) > 1){
    stop("Only one file should be provided when inputFilesType = \"TSStable\"!")
  }
  if(file.exists(inputFiles) == FALSE){
    stop("Could not locate input file ", inputFiles)
  }

  TSS.all.samples <- read.table(file = inputFiles, header = T, stringsAsFactors = FALSE
                                ,colClasses = c("character", "integer", "character", rep("integer", length(sampleLabels)))
                                ,col.names = c("chr", "pos", "strand", sampleLabels))
  if(ncol(TSS.all.samples) != (length(sampleLabels) + 3)){
    stop("Number of provided sample labels must match the number of samples in the TSS table!")
  }
  setDT(TSS.all.samples)
  TSS.all.samples[,4:ncol(TSS.all.samples)][is.na(TSS.all.samples[,4:ncol(TSS.all.samples)])] =0
  return(TSS.all.samples)
}

