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
##.getTSS_from_bam function calls TSS from bam file
##.getTSS_from_bam function takes two input files, bam and Genome
##Genome is BSgenome
##run script with the following example command:
##.getTSS_from_bam(bam, Genome, sampleLabels,sequencingQualityThreshold = 10,
##                            mappingQualityThreshold = 20,
##                            removeNewG = TRUE,
##                            correctG = TRUE)

.getTSS_from_bam <- function(bam.files, Genome, sampleLabels, inputFilesType
                             ,sequencingQualityThreshold = 10
                             ,mappingQualityThreshold = 20
                             ,removeNewG = TRUE
                             ,correctG = TRUE){
  what <- c("rname", "strand", "pos", "seq", "qual", "mapq","flag","cigar")
  param <- ScanBamParam( what = what
                         , flag = scanBamFlag(isUnmappedQuery = FALSE)
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
      qa.avg <- c(qa.avg, as.integer(mean(as(qual[start:end], "IntegerList"))))
      if (end == length(qual)) {
        break
      } else {
        start <- end + 1
      }
    }
    ##repeat cigar, correct intron length proble
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
    readsGR <- GRanges(seqnames = as.vector(bam[[1]]$rname), IRanges(start = bam[[1]]$pos, width = mapped.length), strand = bam[[1]]$strand, qual = qa.avg, mapq = bam[[1]]$mapq, seq = bam[[1]]$seq, read.length = width(bam[[1]]$seq), flag = bam[[1]]$flag)
    readsGR <- readsGR[seqnames(readsGR) %in% seqnames(Genome)]
    readsGR <- readsGR[!(end(readsGR) > seqlengths(Genome)[as.character(seqnames(readsGR))])]
    elementMetadata(readsGR)$mapq[is.na(elementMetadata(readsGR)$mapq)] <- Inf
    readsGR.p <- readsGR[(as.character(strand(readsGR)) == "+" & elementMetadata(readsGR)$qual >= sequencingQualityThreshold) & elementMetadata(readsGR)$mapq >= mappingQualityThreshold]
    readsGR.m <- readsGR[(as.character(strand(readsGR)) == "-" & elementMetadata(readsGR)$qual >= sequencingQualityThreshold) & elementMetadata(readsGR)$mapq >= mappingQualityThreshold]
    if(removeNewG == TRUE){
      TSS <- .removeNewG(readsGR.p, readsGR.m, Genome, correctG = correctG)
    }else{
      TSS.p <- data.frame(chr = as.character(seqnames(readsGR.p)), pos = as.integer(start(readsGR.p)), strand = rep("+", times = length(readsGR.p)), stringsAsFactors = F)
      TSS.m <- data.frame(chr = as.character(seqnames(readsGR.m)), pos = as.integer(end(readsGR.m)), strand = rep("-", times = length(readsGR.m)), stringsAsFactors = F)
      TSS <- rbind(TSS.p, TSS.m)
      TSS$tag_count <- 1
      TSS <- data.table(TSS)
      TSS <- TSS[, as.integer(sum(tag_count)), by = list(chr, pos, strand)]
    }
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
.removeNewG <- function(readsGR.p, readsGR.m, Genome, correctG = TRUE) {
  message("\t-> Removing the first base of the reads if mismatched 'G'...")
  G.reads.p <- which(substr(elementMetadata(readsGR.p)$seq, start = 1, stop = 1) == "G")
  G.reads.m <- which(substr(elementMetadata(readsGR.m)$seq, start = elementMetadata(readsGR.m)$read.length, stop = elementMetadata(readsGR.m)$read.length) == "C")
  if(length(G.reads.p)>0){
    G.mismatch.p <- G.reads.p[getSeq(Genome, resize(readsGR.p[G.reads.p], width = 1, fix = "start"), as.character = TRUE) != "G"]
    elementMetadata(readsGR.p)$removedG <- FALSE
    elementMetadata(readsGR.p)$removedG[G.mismatch.p] <- TRUE
    start(readsGR.p)[G.mismatch.p] <- start(readsGR.p)[G.mismatch.p] + as.integer(1)
    TSS.p <- data.frame(chr = as.character(seqnames(readsGR.p)), pos = start(readsGR.p), strand = "+", removedG = elementMetadata(readsGR.p)$removedG, stringsAsFactors = FALSE)
  }else{
    G.mismatch.p <- NULL
    TSS.p <- data.frame()
  }
  if(length(G.reads.m)>0){
    G.mismatch.m <- G.reads.m[getSeq(Genome, resize(readsGR.m[G.reads.m], width = 1, fix = "start"), as.character = TRUE) != "G"]
    elementMetadata(readsGR.m)$removedG <- FALSE
    elementMetadata(readsGR.m)$removedG[G.mismatch.m] <- TRUE
    end(readsGR.m)[G.mismatch.m] <- end(readsGR.m)[G.mismatch.m] - as.integer(1)
    TSS.m <- data.frame(chr = as.character(seqnames(readsGR.m)), pos = end(readsGR.m), strand = "-", removedG = elementMetadata(readsGR.m)$removedG, stringsAsFactors = FALSE)
  }else{
    G.mismatch.m <- NULL
    TSS.m <- data.frame()
  }
  if(correctG){
    message("\t-> Correcting the systematic bias...")
    G.chance <- (length(G.mismatch.p) + length(G.mismatch.m)) / ((length(readsGR.p) - length(G.reads.p)) + (length(readsGR.m) - length(G.reads.m)) + length(G.mismatch.p) + length(G.mismatch.m))
    if(nrow(TSS.p)>0){
      TSS.G.p <- TSS.p[G.reads.p,]
      TSS.G.p.corrected <- lapply(as.list(unique(TSS.G.p$chr)), function(x) {tss.corrected <- .newGcorrect(tss = subset(TSS.G.p, chr == x), G.chance = G.chance, correction.orientation = 1); tss.corrected$chr = x; return(tss.corrected)})
      TSS.G.p.corrected <- do.call(rbind, TSS.G.p.corrected)
      TSS.G.p.corrected$strand <- "+"
      TSS.G.p.corrected <- TSS.G.p.corrected[,c("chr", "pos", "strand", "nr_tags")]
      TSS.no.G.p <- data.table(TSS.p[-G.reads.p,])
      TSS.no.G.p <- TSS.no.G.p[, length(removedG), by = list(chr, pos, strand)]
      setnames(TSS.no.G.p, c("chr", "pos", "strand", "nr_tags"))
      TSS.p.final <- rbind(TSS.G.p.corrected, as.data.frame(TSS.no.G.p))
      TSS.p.final <- data.table(TSS.p.final)
      TSS.p.final <- TSS.p.final[, sum(nr_tags), by = list(chr, pos, strand)]
    }else{
      TSS.p.final <- data.table()
    }
    if(nrow(TSS.m)>0){
      TSS.G.m <- TSS.m[G.reads.m,]
      TSS.G.m.corrected <- lapply(as.list(unique(TSS.G.m$chr)), function(x) {tss.corrected <- .newGcorrect(tss = subset(TSS.G.m, chr == x), G.chance = G.chance, correction.orientation = -1); tss.corrected$chr = x; return(tss.corrected)})
      TSS.G.m.corrected <- do.call(rbind, TSS.G.m.corrected)
      TSS.G.m.corrected$strand <- "-"
      TSS.G.m.corrected <- TSS.G.m.corrected[,c("chr", "pos", "strand", "nr_tags")]
      TSS.no.G.m <- data.table(TSS.m[-G.reads.m,])
      TSS.no.G.m <- TSS.no.G.m[, length(removedG), by = list(chr, pos, strand)]
      setnames(TSS.no.G.m, c("chr", "pos", "strand", "nr_tags"))
      TSS.m.final <- rbind(TSS.G.m.corrected, as.data.frame(TSS.no.G.m))
      TSS.m.final <- data.table(TSS.m.final)
      TSS.m.final <- TSS.m.final[, sum(nr_tags), by = list(chr, pos, strand)]
    }else{
      TSS.m.final <- data.table()
    }
    TSS <- data.table(rbind(as.data.frame(TSS.p.final), as.data.frame(TSS.m.final)))
  }else{
    TSS <- rbind(TSS.p, TSS.m)
    TSS <- TSS[,c("chr", "pos", "strand")]
    TSS$tag_count <- 1
    TSS <- data.table(TSS)
    TSS <- TSS[, as.integer(sum(tag_count)), by = list(chr, pos, strand)]		
  }
  return(TSS)
}

################################################################################################
## Correcting systematic G nucleotide addition bias to CAGE tags (as described in Carninci et al., Nature Genetics 2006, Supplementary Information, section 3-e)

.newGcorrect <- function(tss, G.chance, correction.orientation) {
  tss.dt <- data.table(tss)
  tss.count <- tss.dt[, list(length(removedG), sum(removedG)), by = pos]
  setkey(tss.count, pos)
  tss.to.correct <- data.frame(subset(tss.count, V1 != V2))
  if(nrow(tss.to.correct) > 0){
    if(correction.orientation > 0){
      tss.gap <- c(Inf, diff(tss.to.correct$pos))
    }else if(correction.orientation < 0){
      tss.gap <- c(diff(tss.to.correct$pos), Inf)
    }
    G.start <- which(tss.gap != 1)
    G.follow <- which(tss.gap == 1)
    tss.to.append <- data.frame()
    while(length(G.start) > 0) {
      F <- as.integer(pmax(round(tss.to.correct$V1[G.start] - tss.to.correct$V2[G.start]/G.chance), 0))
      tss.to.correct$V1[G.start] <- tss.to.correct$V1[G.start] - F
      idx <- G.start + correction.orientation
      if(correction.orientation > 0){
        idx[idx == (nrow(tss.to.correct) + 1)] <- 1			
      }else if(correction.orientation < 0){
        idx[idx == 0] <- nrow(tss.to.correct)
      }
      G.start.followed <- tss.to.correct$pos[idx] %in% tss.to.correct$pos[G.follow]
      tss.to.correct$V1[G.start[G.start.followed] + correction.orientation] <- tss.to.correct$V1[G.start[G.start.followed] + correction.orientation] + F[G.start.followed]
      tss.to.correct$V2[G.start[G.start.followed] + correction.orientation] <- F[G.start.followed]
      tss.to.append <- rbind(tss.to.append, data.frame(pos = tss.to.correct$pos[G.start[!G.start.followed]] + correction.orientation, V1 = F[!G.start.followed], V2 = F[!G.start.followed]))
      G.start <- (G.start + correction.orientation)[G.start.followed]
    }
    tss.final <- rbind(tss.to.correct, tss.to.append)
    tss.final <- rbind(tss.final, as.data.frame(tss.count[V1 == V2]))
    tss.final <- tss.final[order(tss.final$pos),]
    tss.final <- data.frame(pos = tss.final$pos, nr_tags = tss.final$V1)
  }else{
    tss.final <- data.frame(pos = tss.count$pos, nr_tags = tss.count$V1)
  }
  return(subset(tss.final, nr_tags > 0))
}
################################################################################################
##.getTSS_from_bed function calls TSS from bed files
##.getTSS_from_tss function takes bed files and Genome
##run script with the following example command:
##.getTSS_from_bed(inputFiles, Genome)

.getTSS_from_bed <- function(inputFiles, Genome, sampleLabels){
  first <- TRUE
  for(i in 1:length(inputFiles)) {
    message("\nReading in file: ", inputFiles[i], "...")
    readsGR <- import(inputFiles[i], format = "BED")
    readsGR <- readsGR[seqnames(readsGR) %in% seqnames(Genome)]
    readsGR <- readsGR[!(end(readsGR) > seqlengths(Genome)[as.character(seqnames(readsGR))])]
    readsGR.p <- readsGR[strand(readsGR) == "+"]
    readsGR.m <- readsGR[strand(readsGR) == "-"]
    message("\t-> Making TSS table...")
    TSS.plus <- data.frame(chr = as.character(seqnames(readsGR.p)), pos = as.integer(start(readsGR.p)), strand = rep("+", times = length(readsGR.p)), stringsAsFactors = F)
    TSS.minus <- data.frame(chr = as.character(seqnames(readsGR.m)), pos = as.integer(end(readsGR.m)), strand = rep("-", times = length(readsGR.m)), stringsAsFactors = F)
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
  return(TSS.all.samples)
}
################################################################################################
##.getTSS_from_BigWig function calls TSS from BigWig files
##.getTSS_from_tss function takes BigWig files and Genome
##run script with the following example command:
##.getTSS_from_BigWig(inputFiles, Genome)

.getTSS_from_BigWig <- function(inputFiles, Genome, sampleLabels){
  library.sizes <- vector()
  first <- TRUE
  for(i in 1:length(inputFiles)) {
    message("\nReading in file: ", inputFiles[i], "...")
    readsGR <- import(inputFiles[i], format = "BigWig")
    readsGR <- readsGR[seqnames(readsGR) %in% seqnames(Genome)]
    readsGR <- readsGR[!(end(readsGR) > seqlengths(Genome)[as.character(seqnames(readsGR))])]
    readsGR.p <- readsGR[strand(readsGR) == "+"]
    readsGR.m <- readsGR[strand(readsGR) == "-"]
    message("\t-> Making TSS table...")
    TSS.plus <- data.frame(chr = as.character(seqnames(readsGR.p)), pos = as.integer(start(readsGR.p)), strand = rep("+", times = length(readsGR.p)), stringsAsFactors = F)
    TSS.minus <- data.frame(chr = as.character(seqnames(readsGR.m)), pos = as.integer(end(readsGR.m)), strand = rep("-", times = length(readsGR.m)), stringsAsFactors = F)
    TSS <- rbind(TSS.plus, TSS.minus)
    TSS$tag_count <- 1
    TSS <- data.table(TSS)
    TSS <- TSS[, as.integer(sum(tag_count)), by = list(chr, pos, strand)]
    setnames(TSS, c("chr", "pos", "strand", sampleLabels[i]))
    setkey(TSS, chr, pos, strand)
    
    library.sizes <- c(library.sizes, as.integer(sum(data.frame(TSS)[,4])))
    if(first == TRUE) {
      TSS.all.samples <- TSS
    }else{
      TSS.all.samples <- merge(TSS.all.samples, TSS, all = TRUE)
    }
    first <- FALSE
  }   
  return(TSS.all.samples)
}


################################################################################################
##.getTSS_from_tss function calls TSS from bam file
##.getTSS_from_tss function takes tss files
##run script with the following example command:
##.getTSS_from_tss(inputFiles)

.getTSS_from_tss <- function(inputFiles, sampleLabels){
  first <- TRUE
  for(i in 1:length(inputFiles)) {
    message("\nReading in file: ", inputFiles[i], "...")
    TSS <- read.table(file = inputFiles[i], header = F, sep = "\t"
                      ,colClasses = c("character", "integer", "character", "integer")
                      ,col.names = c("chr", "pos", "strand", sampleLabels[i]))
    TSS <- data.table(TSS)
    setkeyv(TSS, cols = c("chr", "pos", "strand"))
    if(first == TRUE) {
      TSS.all.samples <- TSS
    }else{
      TSS.all.samples <- merge(TSS.all.samples, TSS, all = TRUE)
    }			
    first <- FALSE
  }
  TSS.all.samples <- data.frame(TSS.all.samples)
  return(TSS.all.samples)
}

################################################################################################
##.getTSS_from_TSStable function calls TSS from one TSStable file
##.getTSS_from_tss function takes one TSS table file
##run script with the following example command:
##.getTSS_from_TSStable(inputFiles)

.getTSS_from_TSStable <- function(inputFiles){
  if(length(inputFiles) > 1){
    stop("Only one file should be provided when inputFilesType = \"TSStable\"!")
  }
  if(file.exists(inputFiles) == FALSE){
    stop("Could not locate input file ", inputFiles)
  }			
  TSS.all.samples <- read.table(file = inputFiles, header = F, stringsAsFactors = FALSE)
  if(ncol(TSS.all.samples) != (length(sampleLabels) + 3)){
    stop("Number of provided sample labels must match the number of samples in the TSS table!")
  }
  return(TSS.all.samples)
}

