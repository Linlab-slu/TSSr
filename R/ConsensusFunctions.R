##############################################################################################################
##
.getConsensus <- function(gr1,cy, dis){
  colnames(cy)[3:4] <- c("start.c","end.c")
  cy[,start := dominant_tss-round(dis/2)]
  cy[,end := dominant_tss + round(dis/2)]
  gr2 <- makeGRangesFromDataFrame(cy, keep.extra.columns= F)
  hits <- findOverlaps(gr1, gr2)
  new <- c(union(gr1[queryHits(hits)], gr2[subjectHits(hits)]),gr1[-queryHits(hits)],gr2[-subjectHits(hits)])
  return(new)
}

##############################################################################################################
##
.getConsensusQuantile <- function(tc, gr, tss.temp){
  tss_clusters <- lapply(seq_len(gr[,.N]), function(x) {
    data <- gr[x,]
    temp <- tc[chr == gr[x,chr] & strand == gr[x,strand] & dominant_tss >= gr[x, start] & dominant_tss <= gr[x,end],]
    if(nrow(temp) >0){
      s <- tss.temp[chr == temp[1,chr] & strand == temp[1, strand] & pos >= temp[,min(start)] & pos <= temp[,max(end)],]
      q1 <- s[which(cumsum(tags) > 0.1*tags),min(pos)]
      q9 <- s[order(-pos)][which(cumsum(tags) > 0.1*tags),max(pos)]
      list(gr[x, consensusCluster]
           ,gr[x, chr[[1]]]
           ,min(s[,pos])
           ,max(s[,pos])
           ,gr[x, strand[[1]]]
           ,s[which.max(tags),pos]
           ,s[,sum(tags)]
           ,s[,max(tags)]
           ,q1
           ,q9
           ,q9 - q1 + 1)
    }
  })
  tss_clusters <- rbindlist(tss_clusters)
  setnames(tss_clusters, c( "cluster"
                            , "chr", "start", "end", "strand"
                            , "dominant_tss", "tags", "tags.dominant_tss"
                            , "q_0.1", "q_0.9", "interquantile_width" ))
  return(tss_clusters)
}