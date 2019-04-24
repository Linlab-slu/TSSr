###################################################################################################################
###################################################################################################################
##.plotTSS function plots TSS graph
##.plotTSS function takes three input files, tss.tpm, clusters, and gene.list
##tss.tpm table has at least 4 columns (chr, pos, strand, tpm). tpm value is negative in minus strand
##clusters table has 11 columns (cluster,chr,start,end,strand,dominant_tss,tpm,tpm.dominant_tss,q_0.1,q_0.9,interquantile_width)
##gene.list table has at least 5 columns (gene,chr, start, end, strand)
##run script with the following example command:
##.plotTSS(tss.tpm,clusters.cl, gene.list, up.dis = 500, down.dis=100)

##sampleLables <- c("YPD","Arrest","DD","DSA","Gal","Glc","H2O2","HS","NaCl")

.plotTSS <- function(tss, clusters, gene.list, sampleLables, up.dis = 500, down.dis=100){
  tss.p <- tss[tss$strand == "+",]
  tss.m <- tss[tss$strand == "-",]
  tss.m[,c(4:(3+length(sampleLables)))] = tss.m[,c(4:(3+length(sampleLables)))] * (-1)
  tss <- rbind(tss.p,tss.m)
  
    pdf(file = "TSS graph.pdf", height = 7, width = 11)
  for (my.gene in 1:nrow(gene.list)){
    message("Creating graph ",my.gene," ",gene.list$gene[my.gene],"...")
    df <- gene.list[my.gene,]
    if(df$strand == "+"){
      p <- df$start - up.dis
      q <- df$end + down.dis
    }else{
      p <- df$start - down.dis
      q <- df$end + up.dis
    }
    ##Genome range track
    range.gr <- GRanges(seqnames = as.character(df$chr),ranges = IRanges(start = p,end = q),strand = as.character(df$strand))
    gtrack <- GenomeAxisTrack(range = range.gr,fontcolor = "black")
    ##GeneRegion track
    gene.gr <- makeGRangesFromDataFrame(df,keep.extra.columns=FALSE,ignore.strand=FALSE,seqinfo=NULL,seqnames.field=c("chr"),
                                         start.field="start",end.field=c("end"),strand.field="strand",starts.in.df.are.0based=FALSE)
    atrack.gene <- GeneRegionTrack(gene.gr, name = "gene",col.title="black")
    ##clusters track
    clusters_sub <- clusters[clusters$chr == as.character(df$chr) & clusters$strand == as.character(df$strand) & clusters$q_0.1 >= p & clusters$q_0.9 <= q,]
    clusters.gr <- makeGRangesFromDataFrame(clusters_sub,keep.extra.columns=FALSE,ignore.strand=FALSE,seqinfo=NULL,seqnames.field=c("chr"),
                                             start.field="q_0.1",end.field=c("q_0.9"),strand.field="strand",starts.in.df.are.0based=FALSE)
    atrack.clusters <- AnnotationTrack(clusters.gr, name = "clusters",col.title="black")
    ##Data track
    tss_sub <- tss[tss$chr == as.character(df$chr) & tss$strand == as.character(df$strand) & tss$pos >= p & tss$pos <= q,]
    data_tss_track <- list()
    dtrack <- list()
    for (my.condition in 1:length(sampleLables)){
      data_tss_track[[my.condition]] <- makeGRangesFromDataFrame(tss_sub[,c(1:3,(my.condition+3))],keep.extra.columns=TRUE,ignore.strand=FALSE,seqinfo=NULL,seqnames.field=c("chr"),
                                                       start.field="pos",end.field=c("pos"),strand.field="strand",starts.in.df.are.0based=FALSE)
      dtrack[[my.condition]] <- DataTrack(data_tss_track[[my.condition]], name = sampleLables[my.condition],type = "h",col = rainbow(length(sampleLables))[my.condition],baseline = 0, 
                               col.baseline = "grey",col.title="black",col.axis = "black")
    }
    ##plot Genome range track, gene track, clusters track, TSS track
    plotTracks(c(list(gtrack, atrack.gene,atrack.clusters),dtrack),main = gene.list$gene[my.gene])
  }
  dev.off()
}
###################################################################################################################
