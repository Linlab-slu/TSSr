#####################################################################################################
#####################################################################################################
##.assign2gene function assign clusters to reference
##.assign2gene function takes two input files, clusters and refernce
##users need to provide organism scientific name like Saccharomyces cerevisiae
##reference could be "gff3" or "gtf" file, or an URL
##clusters table has 11 columns (cluster,chr,start,end,strand,dominant_tss,tpm,tpm.dominant_tss,q_0.1,q_0.9,interquantile_width) 
##run script with the following example command:
##.assign2gene(clusters,reference,organim,upstream=1000, upstreamOverlap = 500)


.assign2gene <- function(clusters,ref,upstream, upstreamOverlap, downstream){
  ref[, "subset" := paste0(seqnames,"_",strand)]
  setkey(ref,subset)
  clusters[, "subset" := paste0(chr,"_",strand)]
  setkey(clusters,subset)
  asn <- lapply(as.list(clusters[,unique(subset)]), function(x) {
    cs <- clusters[.(x)]
    cs[,subset:= NULL]
    colnames(cs)[3:4] <- c("start.c","end.c")
    gr <- makeGRangesFromDataFrame(cs, keep.extra.columns= T, start.field = "dominant_tss", end.field = "dominant_tss")
    ref_sub <- subset(ref[ref$subset == x,]) ## why ref_sub <- ref[.(x)] not working when x == "chrM_-"?
    if(nrow(ref_sub) == 0){
      mcols(gr)[,"gene"] <- NA
    }else{
      if(cs$strand[1] == "+"){
        setorder(ref_sub,start)
        ref_sub[,end.b:=data.table::shift(end, 1, fill = 0,type='lag')]
        ref_sub[,width:=data.table::shift(width, 1, fill = 1000,type='lag')]
        ref_sub[,dis := start - end.b]
        ref_sub[,dis:= ifelse(dis < 0, 0,dis)]
        ref_sub[,up:= ifelse(dis > upstream, upstream, 
                             ifelse(dis + width  <= upstreamOverlap, dis + width -1,
                                    ifelse(dis < upstreamOverlap, upstreamOverlap, dis)))]
        ref_sub[,end:= start + downstream] ##start - 1 -> start April10
        ref_sub[,start := end - downstream - up +1]
      }else{
        setorder(ref_sub,end)
        ref_sub[,end.b:=data.table::shift(start, 1, fill = 0,type='lead')]
        ref_sub[(.N),end.b := end + 1000]
        ref_sub[,width:=data.table::shift(width, 1, fill = 1000,type='lead')]
        ref_sub[,dis := end.b - end] ##start -> end April10
        ref_sub[,dis:= ifelse(dis < 0,0,dis)]
        ref_sub[,up:= ifelse(dis > upstream, upstream, 
                             ifelse(dis + width  <= upstreamOverlap, dis + width -1,
                                    ifelse(dis < upstreamOverlap, upstreamOverlap, dis)))]
        ref_sub[,start:= end- downstream]##end + 1 -> end April10
        ref_sub[,end:= start + downstream + up -1]
      }
      rownames(ref_sub) <- ref_sub$gene_id
      ref_sub <- makeGRangesFromDataFrame(ref_sub, keep.extra.columns= F)
      hits <- findOverlaps(gr,ref_sub)
      hits <- breakTies(hits, method = "first")
      hits <- methods::as(hits, "List")
      hits <- extractList(names(ref_sub), hits)
      hits <- as.character(hits)
      mcols(gr)[,"gene"] <- hits
    }
    gr <- as.data.frame(gr)
    gr$dominant_tss <- gr$start
    colnames(gr)[c(1,7,8)] <- c("chr","start","end")
    gr <- gr[,c(6,1,7,8,5,15,9:14)]
    return(gr)
  })
  setorder(do.call("rbind",asn), cluster)
}



#####################################################################################################

