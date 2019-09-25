##########################################################################################################
.plotCorrelation <- function(TSS.all.samples){
  z <- TSS.all.samples[,-c(1:3)]
  # Customize lower panel
  panel.cor <- function(x, y){
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1), xlog = FALSE, ylog = FALSE)
    r <- round(cor(x, y), digits=2)
    txt <- paste0( r)
    cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * r)
  }
  # Customize upper panel
  upper.panel<-function(x, y){
    points(x,y, pch = ".",col = "#00AFBB")
  }
  # Create the plots
  suppressWarnings(pairs(z, lower.panel = upper.panel,upper.panel = panel.cor, log = "xy"))
}

###################################################################################################################
###################################################################################################################
##.plotTSS function plots TSS graph
##.plotTSS function takes three input files, tss.tpm, clusters, and ref
##tss.tpm table has at least 4 columns (chr, pos, strand, tpm). tpm value is negative in minus strand
##clusters table has 11 columns (cluster,chr,start,end,strand,dominant_tss,tpm,tpm.dominant_tss,q_0.1,q_0.9,interquantile_width)
##ref table has at least 5 columns (gene,chr, start, end, strand)
##run script with the following example command:
##.plotTSS(tss.tpm,clusters.cl, ref, up.dis = 500, down.dis=100)
.plotTSS <- function(tss, clusters,df, samples, up.dis, down.dis){
  setnames(df, colnames(df)[c(1,6)], c("chr","gene"))
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
  for (my.condition in 1:length(samples)){
    temp <- tss_sub[,.SD, .SDcols = c("chr","pos","strand",samples[my.condition])]
    data_tss_track[[my.condition]] <- makeGRangesFromDataFrame(temp,keep.extra.columns=TRUE,ignore.strand=FALSE,seqinfo=NULL,seqnames.field=c("chr"),
                                                               start.field="pos",end.field=c("pos"),strand.field="strand",starts.in.df.are.0based=FALSE)
    dtrack[[my.condition]] <- DataTrack(data_tss_track[[my.condition]], name = samples[my.condition],type = "h",col = rainbow(length(samples))[my.condition],baseline = 0, 
                                        col.baseline = "grey",col.title="black",col.axis = "black")
  }
  ##plot Genome range track, gene track, clusters track, TSS track
  plotTracks(c(list(gtrack, atrack.gene,atrack.clusters),dtrack),main = df$gene)
}
###################################################################################################################

.plotInterQuantile <- function(tc, tagsThreshold){
  temp <- tc[tags >= tagsThreshold ]
  
}