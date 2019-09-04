#####################################################################################################
#####################################################################################################
##.pss function calculates promoter shape with PSS equation (as described in Lu & Lin, bioRxiv 450429;  doi: https://doi.org/10.1101/450429, and modified)
##.pss function takes two input files, tss.tpm and clusters
##tss.tpm table has 4 columns (chr, pos, strand, tpm)
##clusters table has at least 11 columns (cluster,chr,start,end,strand,dominant_tss,tpm,tpm.dominant_tss,q_0.1,q_0.9,interquantile_width) 
##clusters could be clusters.asn
##run script with the following example command:
##.pss(tss.tpm,clusters,useMultiCore= TRUE)



##########################################################################################################
tss <- read.table("35_39.tpm0.1.txt",header = T)
colnames(tss)[4] = "tags"
cs <- read.table("homo.clusters.assigned-500+500.txt",header = T)
ps <- .pss(tss,cs, useMultiCore= TRUE, numCores = NULL)

.pss <- function(ctss.tpm, clusters, useMultiCore= TRUE, numCores = NULL){
  setDT(ctss.tpm)
  setDT(clusters)
  clusters$chr <- as.character(clusters$chr)
  clusters$strand <- as.character(clusters$strand)
  ctss.tpm$tags <- as.numeric(ctss.tpm$tags)
  setkey(clusters,cluster)
  if (useMultiCore) {
    library(parallel)
    if (is.null(numCores)) {
      numCores <- detectCores()
    }
    print(paste("process is running on", numCores, "cores..."))
    
    cs.pss <- mclapply(seq_len(clusters[,.N]), function(x) {
      data <- clusters[x,]
      ctss <- ctss.tpm[chr == data$chr & strand == data$strand & pos >= data$q_0.1 & pos <= data$q_0.9,]
      temp <- sum(ctss[,tags])
      data$pss <- -sum(sapply(ctss[,tags],function(y){y/temp*log(y/temp,2)}))*log(clusters[x,interquantile_width],2)
      return(data)
      }, mc.cores = numCores)
    }
  else{
    cs.pss <- lapply(seq_len(clusters[,.N]), function(x) {
      data <- clusters[x,]
      ctss <- ctss.tpm[chr == data$chr & strand == data$strand & pos >= data$q_0.1 & pos <= data$q_0.9,]
      temp <- sum(ctss[,tags])
      data$pss <- -sum(sapply(ctss[,tags],function(y){y/temp*log(y/temp,2)}))*log(clusters[x,interquantile_width],2)
      return(data)
    })
  }
  cs.pss <- rbindlist(cs.pss)
  setDF(cs.pss)
  return(cs.pss)
} 


