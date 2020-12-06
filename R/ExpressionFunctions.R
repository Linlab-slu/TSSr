###################################################################################################################
##.deseq2 function calcuates gene differential expression based on Deseq2 package
##.deseq2 function takes two assigned clusters and tss.raw table
##users need to provide which sample they want to compare and
##run script with the following example command:
##.deseq2(clustersx.asn,clustersy.asn, tss.raw,
##                              samplex <- c("ScerBY4741.1","ScerBY4741.2"),
##                              sampley <- c("ScerArrest.1","ScerArrest.2"),
##                              sampleOne <- "ScerBY4741",sampleTwo <- "ScerArrest")
############################################################################
##tss.raw is the raw tss merged tables, before any sums
.deseq2 <- function(cx,cy, tss.raw, samplex,sampley, sampleOne,sampleTwo,useMultiCore, numCores){
  ##get raw count tables
  xCounts <-.tagCount(cx, tss.raw,samplex,useMultiCore, numCores)
  yCounts <-.tagCount(cy, tss.raw,sampley,useMultiCore, numCores)
  xCounts <- xCounts[,-c(2:11)]
  yCounts <- yCounts[,-c(2:11)]

  ##tag counts by gene for sampleOne
  setkey(xCounts, gene)
  if(useMultiCore){
    if(is.null(numCores)){
      numCores <- detectCores()
    }
    one <- mclapply(as.list(unique(xCounts$gene)), function(my.gene) {
      data <- xCounts[list(my.gene)]
      return(c(my.gene,colSums(data[,-c(1,2,3)])))
    }, mc.cores = numCores)
  }else{
    one <- lapply(as.list(unique(xCounts$gene)), function(my.gene) {
      data <- xCounts[list(my.gene)]
      return(c(my.gene,colSums(data[,-c(1,2,3)])))
    })
  }
  one <- data.frame(matrix(unlist(one), nrow=length(one), byrow=T),stringsAsFactors=FALSE)

  ##tag counts by gene for sampleTwo
  setkey(yCounts, gene)
  if(useMultiCore){
    if(is.null(numCores)){
      numCores <- detectCores()
    }
    two <- mclapply(as.list(unique(yCounts$gene)), function(my.gene) {
      data <- yCounts[list(my.gene)]
      return(c(my.gene,colSums(data[,-c(1,2,3)])))
    }, mc.cores = numCores)

  }else{
    two <- lapply(as.list(unique(yCounts$gene)), function(my.gene) {
      data <- yCounts[list(my.gene)]
      return(c(my.gene,colSums(data[,-c(1,2,3)])))
    })
  }
  two <- data.frame(matrix(unlist(two), nrow=length(two), byrow=T),stringsAsFactors=FALSE)
  ##merge the two raw count tables together by genes
  one[,2:ncol(one)] <- sapply(one[,2:ncol(one)], as.integer)
  two[,2:ncol(two)] <- sapply(two[,2:ncol(two)], as.integer)
  setnames(one, colnames(one), c("gene",samplex))
  setnames(two, colnames(two), c("gene",sampley))
  Dtable <- merge(one, two, by = c("gene"), all = T)
  Dtable[is.na(Dtable)] = 0
  ##
  rownames(Dtable) <- Dtable[,1]
  Dtable <- Dtable[,-1]
  Dtable <- data.matrix(Dtable)
  condition <- factor(c(rep(sampleOne, times = length(samplex)), rep(sampleTwo, times = length(sampley))))
  dds <- DESeqDataSetFromMatrix(countData = Dtable,data.frame(condition), ~ condition)
  dds$condition <- factor(dds$condition, levels = c(sampleOne, sampleTwo))
  dds <- DESeq(dds)
  res <- results(dds)
  res <- res[order(res$padj),]
  return(as.data.frame(res))
}


############################################################################
##.tagCount is slow
.tagCount <- function(cs, tss.raw, samples, useMultiCore, numCores){
  cols <- c("chr","pos","strand", samples)
  tss <- tss.raw[,.SD, .SDcols = cols]
  ##define variable as a NULL value
  chr = strand = start = end = NULL

  if(useMultiCore){
    if(is.null(numCores)){
      numCores <- detectCores()
    }
    print(paste("process is running on",numCores, "cores..."))
    tags <- mclapply(seq_len(cs[,.N]),function(r){
      data <- tss[tss$chr == cs[r,chr] & tss$strand == cs[r,strand] & tss$pos >= cs[r,start] & tss$pos <= cs[r,end],]
      temp <- sapply(as.list(samples), function(s){
        sum(data[,.SD,.SDcols = s])
      })
      return(temp)
    }, mc.cores = numCores)

  }else{
    tags <- lapply(seq_len(cs[,.N]),function(r){
      data <- tss[tss$chr == cs[r,chr] & tss$strand == cs[r,strand] & tss$pos >= cs[r,start] & tss$pos <= cs[r,end],]
      temp <- sapply(as.list(samples), function(s){
        sum(data[,.SD,.SDcols = s])
      })
      return(temp)
    })
  }
  tags <- data.frame(matrix(unlist(tags), nrow=length(tags), byrow=T),stringsAsFactors=FALSE)
  colnames(tags) <- samples
  cs <- cbind(cs,tags)
  return(cs)
}
