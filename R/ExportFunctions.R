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
##.plotTSS function takes three input files, tss.tpm, cs, and ref
##tss.tpm table has at least 4 columns (chr, pos, strand, tpm). tpm value is negative in minus strand
##cs table has 11 columns (cluster,chr,start,end,strand,dominant_tss,tpm,tpm.dominant_tss,q_0.1,q_0.9,interquantile_width)
##ref table has at least 5 columns (gene,chr, start, end, strand)
##run script with the following example command:
##.plotTSS(tss.tpm,cs.cl, ref, up.dis = 500, down.dis=100)
.plotTSS <- function(tss, cs,df, samples, Bidirection, up.dis, down.dis,yFixed){
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
  gene.gr <- makeGRangesFromDataFrame(df,keep.extra.columns=FALSE,ignore.strand=FALSE,
                                      seqinfo=NULL,seqnames.field=c("chr"),
                                      start.field="start",end.field=c("end"),
                                      strand.field="strand",starts.in.df.are.0based=FALSE)
  atrack.gene <- GeneRegionTrack(gene.gr, name = "gene",col.title="black",
                                 transcriptAnnotation = "gene")
  ##cluster and Data track
  if(Bidirection == TRUE){
    tss_sub <- tss[tss$chr == as.character(df$chr)
                   & tss$pos >= p & tss$pos <= q,]
    tss_sub[, strand:= df$strand]
  }else{
    tss_sub <- tss[tss$chr == as.character(df$chr) & tss$strand == as.character(df$strand)
                   & tss$pos >= p & tss$pos <= q,]
  }
  
  # setup y_range
  if(yFixed==TRUE){
    tss_sub1 <- tss_sub[,.SD, .SDcols = c(samples)]
    max(tss_sub1)
    y_range=c(0,max(tss_sub1))
  }else{
    y_range=NULL  
  }
  
  data_tss_track <- list()
  dtrack <- list()
  s <- c()
  for (my.sample in seq(samples)){
    ##clusters track
    cs.temp <- cs[[samples[my.sample]]]
    cs_sub <- cs.temp[cs.temp$chr == as.character(df$chr) & cs.temp$strand == as.character(df$strand) & cs.temp$q_0.1 >= p & cs.temp$q_0.9 <= q,]
    cs.gr <- makeGRangesFromDataFrame(cs_sub,keep.extra.columns=FALSE,ignore.strand=FALSE,seqinfo=NULL,seqnames.field=c("chr"),
                                            start.field="q_0.1",end.field=c("q_0.9"),strand.field="strand",starts.in.df.are.0based=FALSE)
    # atrack.cs <- AnnotationTrack(cs.gr, name = paste(samples[my.sample],"clusters",sep = " "),col.title="black")
    atrack.cs <- AnnotationTrack(cs.gr, name = paste(samples[my.sample],"clusters",sep = " "),
                                 id = cs_sub$cluster, col.title="black")
    ##tss track
    temp <- tss_sub[,.SD, .SDcols = c("chr","pos","strand",samples[my.sample])]
    data_tss_track <- makeGRangesFromDataFrame(temp,keep.extra.columns=TRUE,ignore.strand=FALSE,seqinfo=NULL,seqnames.field=c("chr"),
                                                               start.field="pos",end.field=c("pos"),strand.field="strand",starts.in.df.are.0based=FALSE)
    dtrack <- DataTrack(data_tss_track, name = paste(samples[my.sample],"TSS", sep = " "),type = "h",col = rainbow(length(samples))[my.sample],baseline = 0,
                                        col.baseline = "grey",col.title="black",col.axis = "black",ylim = y_range)
    s <- c(s, atrack.cs, dtrack)
  }
  ##plot Genome range track, gene track, clusters track, TSS track
  plotTracks(c(list(gtrack, atrack.gene),s),main = df$gene,
             featureAnnotation="id",fontcolor.feature="black", cex.feature=0.7, shape = "arrow")
}
###################################################################################################################
.getBed <- function(p){
  ##define variable as a NULL value
  chr = q_0.1 = q_0.9 = dominant_tss = NULL
  pbed <- lapply(as.list(seq_len(p[,.N])), function(i){
    if(p[i,q_0.1] == p[i, q_0.9]){
      nrBlocks = length(unique(c(p[i,start],p[i,end],p[i,q_0.1],p[i,q_0.9])))
    }else{nrBlocks = length(unique(c(p[i,start],p[i,end],p[i,q_0.1],p[i,q_0.9]))) -1}

    if(nrBlocks == 3){
      block.sizes = paste(1, p[i,q_0.9]-p[i,q_0.1]+1, 1, sep = ",")
      block.starts = paste(0, p[i,q_0.1]-p[i,start], p[i,end]-p[i,start], sep = ",")
    }else if(nrBlocks ==2){
      if(p[i,q_0.1] == p[i, start]){
        block.sizes = paste(p[i,q_0.9]-p[i,q_0.1]+1,1, sep = ",")
        block.starts = paste(0, p[i,end]-p[i,start], sep = ",")
      }else if(p[i,end] == p[i, q_0.9]){
        block.sizes = paste(1, p[i,q_0.9]-p[i,q_0.1]+1, sep = ",")
        block.starts = paste(0, p[i,q_0.1]-p[i,start], sep = ",")
      }else{
        message("\nWhat is the additional condition which has two blocks...")
        print(p[i,])
      }
    }else{
      block.sizes = paste(p[i,end]-p[i,start])
      block.starts = paste(0)
    }
    list(p[i,chr]
         ,p[i,start]-1
         ,p[i,end]
         ,"."
         ,0
         ,p[i,strand]
         ,p[i,dominant_tss]-1
         ,p[i,dominant_tss]
         ,0
         ,nrBlocks
         ,block.sizes
         ,block.starts)
  })
  pbed <- rbindlist(pbed)
  return(pbed)
}
