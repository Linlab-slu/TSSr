###################################################################################################################
##.deseq2 function calcuates gene differential expression based on Deseq2 package
##.deseq2 function takes two assigned clusters and TSS.all.samples table
##users need to provide which sample they want to compare and 
##run script with the following example command:           
##.deseq2(clustersx.asn,clustersy.asn, TSS.all.samples, 
##                              samplex <- c("ScerBY4741.1","ScerBY4741.2"),
##                              sampley <- c("ScerArrest.1","ScerArrest.2"),
##                              sampleMergedx <- "ScerBY4741",sampleMergedy <- "ScerArrest")
############################################################################
##TSS.all.samples is the raw tss merged tables, before any sums
.deseq2 <- function(clustersx.asn,clustersy.asn, TSS.all.samples, samplex,sampley, sampleMergedx,sampleMergedy){
  setDT(clustersx.asn)
  setDT(clustersy.asn)
  setDT(TSS.all.samples)
  ##clustersx.asn is control 
  ##clustersy.asn is test
  cx <- clustersx.asn[!is.na(clustersx.asn$gene),]
  cy <- clustersy.asn[!is.na(clustersy.asn$gene),]
  TSS.all.samples$chr = as.character(TSS.all.samples$chr)
  TSS.all.samples$strand = as.character(TSS.all.samples$strand)
  cx$chr = as.character(cx$chr)
  cx$strand = as.character(cx$strand)
  cy$chr = as.character(cy$chr)
  cy$strand = as.character(cy$strand)
  ##get raw count tables
  xCounts <-.tagCount(cx, samplex)
  yCounts <-.tagCount(cy, sampley)
  xCounts <- xCounts[,-c(2:11)]
  yCounts <- yCounts[,-c(2:11)]
  ##merge the two raw count tables together
  clusters <- merge(xCounts,yCounts, by = c("cluster"), all = T)
  clusters$gene.x <- as.character(clusters$gene.x)
  clusters$gene.y <- as.character(clusters$gene.y)
  genelist <- lapply(seq_len(clusters[,.N]),function(r){ifelse(is.na(clusters[r,gene.x]), clusters[r,gene.y], clusters[r,gene.x])})
  clusters[,gene:= genelist]
  clusters$gene <- as.character(clusters$gene)
  cols <- c("gene", c(samplex,sampley))
  clusters <- clusters[,.SD,.SDcols = cols]
  clusters[is.na(clusters)]=0
  setkey(clusters,gene)
  ##
  cmp <- lapply(as.list(unique(clusters$gene)), function(my.gene) {
    data <- clusters[.(my.gene)]
    return(c(my.gene,colSums(data[,-1])))
  })
  Dtable <- data.frame(matrix(unlist(cmp), nrow=length(cmp), byrow=T),stringsAsFactors=FALSE)
  colnames(Dtable) <- cols
  rownames(Dtable) <- Dtable[,1]
  Dtable <- Dtable[,-1]
  Dtable <- data.matrix(Dtable)
  condition <- factor(c(rep(sampleMergedx, times = length(samplex)), rep(sampleMergedy, times = length(sampley))))
  dds <- DESeqDataSetFromMatrix(countData = Dtable,DataFrame(condition), ~ condition)
  dds <- DESeq(dds)
  res <- results(dds)
  res <- res[order(res$padj),]
  return(res)
}


############################################################################
##.tagCount is slow
.tagCount <- function(cs, samples){
  cols <- c("chr","pos","strand", samples)
  tss <- TSS.all.samples[,.SD, .SDcols = cols]
  tags <- lapply(seq_len(cs[,.N]),function(r){
    data <- tss[tss$chr == cs[r,chr] & tss$strand == cs[r,strand] & tss$pos >= cs[r,start] & tss$pos <= cs[r,end],]
    temp <- sapply(as.list(samples), function(s){
      sum(data[,.SD,.SDcols = s])
    })
    return(temp)
  })
  tags <- data.frame(matrix(unlist(tags), nrow=length(tags), byrow=T),stringsAsFactors=FALSE)
  colnames(tags) <- samples
  cs <- cbind(cs,tags)
  return(cs)
}
