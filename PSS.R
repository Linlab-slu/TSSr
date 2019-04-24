#####################################################################################################
#####################################################################################################
##.pss function calculates promoter shape with PSS equation (as described in Lu & Lin, bioRxiv 450429;  doi: https://doi.org/10.1101/450429, and modified)
##.pss function takes two input files, tss.tpm and clusters
##tss.tpm table has 4 columns (chr, pos, strand, tpm)
##clusters table has at least 11 columns (cluster,chr,start,end,strand,dominant_tss,tpm,tpm.dominant_tss,q_0.1,q_0.9,interquantile_width) 
##clusters could be clusters.asn
##run script with the following example command:
##.pss(tss.tpm,clusters)



##########################################################################################################


.pss <- function(tss.tpm, clusters){
  setDT(tss.tpm)
  setDT(clusters)
  clusters <- clusters[, chr := as.character(chr)]
  clusters <- clusters[, strand := as.character(strand)]
  tss.tpm$tags <- as.numeric(tss.tpm$tags)
  setkey(clusters,cluster)
  cs.pss <- sapply(seq_len(clusters[,.N]), function(x) {
    tss <- tss.tpm[chr == clusters[x,chr] & strand == clusters[x,strand] & pos >= clusters[x, q_0.1] & pos <= clusters[x,q_0.9],]
    temp <- sum(tss[,tags])
    -sum(sapply(tss[,tags],function(y){y/temp*log(y/temp,2)}))*log(clusters[x,interquantile_width],2)
  })
  clusters[,pss:= cs.pss]
  setDF(clusters)
  return(clusters)
}

