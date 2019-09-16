###################################################################################################################
##.Ds function calcuates gene differential expression based on Deseq2 package
##.Ds function takes two assigned clusters and library sizes of the two samples
##users need to provide which sample they want to compare and 
##run script with the following example command:           
##.Ds(cx,cy, librarySizex, librarySizey, useRawCount = TRUE)

############################################################################################################
.Ds <- function(cx,cy, librarySizex, librarySizey, useRawCount = TRUE){
  ##cx is control
  ##cy is treat
  cx <- cx[,c("cluster","strand","dominant_tss","tpm","gene")]
  cy <- cy[,c("cluster","strand","dominant_tss","tpm","gene")]
  clusters <- merge(cx,cy, by = c("cluster","strand"), all = T)
  clusters[,c("dominant_tss.x","dominant_tss.y","tpm.x","tpm.y")][is.na(clusters[,c("dominant_tss.x","dominant_tss.y","tpm.x","tpm.y")])]=0
  clusters[,dominant_tss := do.call(pmax,clusters[,c("dominant_tss.x","dominant_tss.y")])]
  clusters$gene.x <- as.character(clusters$gene.x)
  clusters$gene.y <- as.character(clusters$gene.y)
  clusters <- clusters[!is.na(clusters$gene.x) | !is.na(clusters$gene.y),]
  genelist <- lapply(seq_len(clusters[,.N]),function(r){ifelse(is.na(clusters[r,gene.x]), clusters[r,gene.y], clusters[r,gene.x])})
  clusters[,gene:= genelist]
  clusters$gene <- as.character(clusters$gene)
  clusters <- clusters[,c("cluster","strand","dominant_tss","gene","tpm.x","tpm.y")]
  setkey(clusters,gene)
  ##
  cmp <- lapply(as.list(unique(clusters$gene)), function(my.gene) {
    data <- clusters[.(my.gene)]
    if(sum(data[,tpm.x]) == 0 | sum(data[,tpm.y])==0){
      Ds <- NA
      pval = NA
    }else if(nrow(data) >= 2){
      setorder(data,-tpm.y)
      ##find the top two clusters which have higher expression level than all the others if present
      ID1 <- data[1,cluster]
      setorder(data,-tpm.x)
      ID2 <- ifelse(ID1==data[1,cluster], data[2,cluster], data[1,cluster])
      data <- data[cluster %in% c(ID1,ID2),]
      if(data$strand[1] == "+"){setorder(data, -dominant_tss)}
      else{setorder(data, dominant_tss)}
      ##now the first row of data is proximal core promoter
      ##the second row of data is the distal cor promoter
      Ds <- log(((data[1,tpm.y]+0.1)/(data[2,tpm.y]+0.1))/((data[1,tpm.x]+0.1)/(data[2,tpm.x]+0.1)),2)
      Xtable <- data[,c("tpm.x","tpm.y")]
      if(useRawCount){
        Xtable$tpm.x <- Xtable$tpm.x * librarySizex/1000000
        Xtable$tpm.y <- Xtable$tpm.y * librarySizey/1000000}
      pval <- suppressWarnings(chisq.test(data[,c("tpm.x","tpm.y")])$p.value)
    }else{
      Ds <- NA
      pval = NA}
    return(c(my.gene,Ds,pval))
  })
  Ds <- data.frame(matrix(unlist(cmp), nrow=length(cmp), byrow=T),stringsAsFactors=FALSE)
  colnames(Ds) <- c("gene","Ds","pval")
  Ds$padj <- p.adjust(Ds[,"pval"],method = "BH")
  return(Ds)
}


