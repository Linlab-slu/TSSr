###################################################################################################################
##.getTSS_from_BigWig function calls TSS from BigWig files
##.getTSS_from_tss function takes BigWig files and genome
##run script with the following example command:
##.getTSS_from_BigWig(BigWig.files, genome)


.getTSS_from_BigWig <- function(BigWig.files, genome){
  library.sizes <- vector()
  first <- TRUE
  for(i in 1:length(BigWig.files)) {
    message("\nReading in file: ", BigWig.files[i], "...")
    readsGR <- import(BigWig.files[i], format = "BigWig")
    readsGR <- readsGR[seqnames(readsGR) %in% seqnames(genome)]
    readsGR <- readsGR[!(end(readsGR) > seqlengths(genome)[as.character(seqnames(readsGR))])]
    readsGR.p <- readsGR[strand(readsGR) == "+"]
    readsGR.m <- readsGR[strand(readsGR) == "-"]
    message("\t-> Making TSS table...")
    TSS.plus <- data.frame(chr = as.character(seqnames(readsGR.p)), pos = as.integer(start(readsGR.p)), strand = rep("+", times = length(readsGR.p)), stringsAsFactors = F)
    TSS.minus <- data.frame(chr = as.character(seqnames(readsGR.m)), pos = as.integer(end(readsGR.m)), strand = rep("-", times = length(readsGR.m)), stringsAsFactors = F)
    TSS <- rbind(TSS.plus, TSS.minus)
    TSS$tag_count <- 1
    TSS <- data.table(TSS)
    TSS <- TSS[, as.integer(sum(tag_count)), by = list(chr, pos, strand)]
    setnames(TSS, c("chr", "pos", "strand", sample.labels[i]))
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