#####################################################################################################
##.assign2gene function assign cs.temp to reference
.assign2gene <- function(cs.temp,ref,upstream, upstreamOverlap, downstream, filterCluster){
  ref[, "subset" := paste0(seqnames,"_",strand)]
  setkey(ref,subset)
  cs.temp[, "subset" := paste0(chr,"_",strand)]
  setkey(cs.temp,subset)
  asn <- lapply(as.list(cs.temp[,unique(subset)]), function(x) {
    cs <- cs.temp[.(x)]
    cs[,subset:= NULL]
    colnames(cs)[3:4] <- c("start.c","end.c")
    gr <- makeGRangesFromDataFrame(cs, keep.extra.columns= T, start.field = "dominant_tss", end.field = "dominant_tss")
    ref_sub <- ref[.(x)]
    ref_coding <- ref[.(x)]
    if(!(x %in% ref[,unique(subset)])){
      mcols(gr)[,"gene"] <- NA
      mcols(gr)[,"inCoding"] <- NA
    }else{
      if(cs$strand[1] == "+"){
        setorder(ref_sub,start)
        ref_sub[,end.b:=data.table::shift(end, 1, fill = 0,type='lag')]
        ref_sub[,width:=data.table::shift(width, 1, fill = 1000,type='lag')]
        ref_sub[,dis := start - end.b]
        ##        ref_sub[,dis:= ifelse(dis < 0, 0,dis)] ##Oct1
        ref_sub[,up:= ifelse(dis > upstream, upstream, 
                             ifelse(dis + width  <= upstreamOverlap, dis + width -1,
                                    ifelse(dis < upstreamOverlap, upstreamOverlap, dis)))]
        ##add down to solve overlap issue
        ref_sub[,start.a:=data.table::shift(start, 1, fill = 1000,type='lead')]
        ref_sub[, dis.start := start.a - start]
        ref_sub[, down.up := data.table::shift(up, 1, fill = 1000,type='lead')]
        ref_sub[,down:= ifelse(dis.start > downstream + down.up, downstream, 
                             ifelse(dis.start  < down.up, 0, dis.start-down.up))]
        ##
        ref_sub[,end:= start + down] ##start - 1 -> start April10
        ref_sub[,start := end - down - up +1]
      }else{
        setorder(ref_sub,end)
        ref_sub[,end.b:=data.table::shift(start, 1, fill = 0,type='lead')]
        ref_sub[(.N),end.b := end + 1000]
        ref_sub[,width:=data.table::shift(width, 1, fill = 1000,type='lead')]
        ref_sub[,dis := end.b - end] ##start -> end April10
        ##        ref_sub[,dis:= ifelse(dis < 0,0,dis)]
        ref_sub[,up:= ifelse(dis > upstream, upstream, 
                             ifelse(dis + width  <= upstreamOverlap, dis + width -1,
                                    ifelse(dis < upstreamOverlap, upstreamOverlap, dis)))]
        ##add down to solve overlap issue
        ref_sub[,start.a:=data.table::shift(end, 1, fill = 1000,type='lag')]
        ref_sub[, dis.start := end -start.a]
        ref_sub[, down.up := data.table::shift(up, 1, fill = 1000,type='lag')]
        ref_sub[,down:= ifelse(dis.start >= downstream+ down.up, downstream, 
                               ifelse(dis.start  < down.up, 0,dis.start - down.up))]
        ##
        ref_sub[,start:= end- down]##end + 1 -> end April10
        ref_sub[,end:= start + down + up -1]
      }
      rownames(ref_sub) <- ref_sub$gene_id
      ref_sub <- makeGRangesFromDataFrame(ref_sub, keep.extra.columns= F)
      ##find overlap with promoter region
      hits <- findOverlaps(gr,ref_sub)
      hits <- breakTies(hits, method = "first")
      hits <- methods::as(hits, "List")
      hits <- extractList(names(ref_sub), hits)
      hits <- as.character(hits)
      mcols(gr)[,"gene"] <- hits
      if(filterCluster == TRUE){
        ##find overlap with coding regions
        ##coding
        rownames(ref_coding) <- ref_coding$gene_id
        ref_coding <- makeGRangesFromDataFrame(ref_coding, keep.extra.columns = F)
        hits <- findOverlaps(gr,ref_coding)
        hits <- breakTies(hits, method = "first")
        hits <- methods::as(hits, "List")
        hits <- extractList(names(ref_coding), hits)
        hits <- as.character(hits)
        mcols(gr)[,"inCoding"] <- hits
      }
    }
    gr <- as.data.frame(gr)
    gr$dominant_tss <- gr$start
    colnames(gr)[c(1,7,8)] <- c("chr","start","end")
    gr <- gr[,c(6,1,7,8,5,ncol(gr),9:(ncol(gr)-1))]
    setDT(gr)
    return(gr)
  })
  setorder(do.call("rbind",asn), cluster)
}
#####################################################################################################

