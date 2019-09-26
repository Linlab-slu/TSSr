################################################################################################
setGeneric("plotCorrelation",function(object,...)standardGeneric("plotCorrelation"))
setMethod("plotCorrelation","TSSr", function(object
                                             ,samples = "all"
){
  message("Plotting TSS correlations...")
  tss.raw <- object@TSSrawMatrix
  if(samples == "all"){
    tss <- tss.raw
  }else{
    cols <- c("chr","pos","strand", samples)
    tss <- tss.raw[,.SD, .SDcols = cols]
  }
  pdf(file = paste("TSS_correlation_plot_of_", paste(samples, collapse = "_"), "_samples.pdf", sep = "")
      ,width = 8, height = 8, onefile = T, bg = "transparent", family = "Helvetica", fonts = NULL)
  .plotCorrelation(tss)
  dev.off()
})

################################################################################################
setGeneric("plotPCA",function(object,...)standardGeneric("plotPCA"))
setMethod("plotPCA","TSSr", function(object
                                     ,TSS.threshold=10
){
  message("Plotting TSS PCA...")
  tss <- object@TSSrawMatrix
  tss <- tss[,4:ncol(tss)]
  tss <- tss[apply(tss>TSS.threshold, 1, any),]
  y <- t(tss)
  sampleLabels <- object@sampleLabels
  sampleLabelsMerged <- object@sampleLabelsMerged
  mergeIndex <- object@mergeIndex
  s <- sampleLabelsMerged[mergeIndex]
  pdf(file = "PCA_plot.pdf"
      ,width = 8, height = 8, onefile = T, bg = "transparent", family = "Helvetica", fonts = NULL)
  print(autoplot(prcomp(y), data = data.frame(sample = s)
           ,colour = "sample", size = 3) + theme_minimal()+theme(text = element_text(size=12)))
  dev.off()
})

################################################################################################
setGeneric("plotInterQuantile",function(object,...)standardGeneric("plotInterQuantile"))
setMethod("plotInterQuantile","TSSr", function(object
                                               ,samples = "all"
                                               ,tagsThreshold = 1
){
  message("Plotting interquantile graphs...")
  TCs <- object@clusterShape
  sampleLabels <- object@sampleLabelsMerged
  if(samples == "all"){
    tc <- TCs
    pdf(file = paste("Interquantile_plot_of_ALL_samples.pdf", sep = "")
        ,width = 8, height = 8, onefile = T, bg = "transparent", family = "Helvetica", fonts = NULL)
    for(i in 1:length(sampleLabels)){
      temp <- tc[[sampleLabels[i]]]
      temp <- temp[tags >= tagsThreshold & interquantile_width <= 200,]
      hist(temp$interquantile_width, breaks = 40, col = rainbow(length(sampleLabels))[i]
           , xlab = "TC interquantile width q0.1-q0.9", ylab = "Frequency", main = sampleLabels[i])
    }
  }else{
    tc <- TCs[[samples]]
    pdf(file = paste("Interquantile_plot_of_", paste(samples, collapse = "_"), "_samples.pdf", sep = "")
        ,width = 8, height = 8, onefile = T, bg = "transparent", family = "Helvetica", fonts = NULL)
    for(i in 1:length(samples)){
      temp <- tc[[samples[i]]]
      temp <- temp[tags >= tagsThreshold & interquantile_width <= 200,]
      hist(temp$interquantile_width, breaks = 40, col = rainbow(length(samples))[i]
           , xlab = "TC interquantile width q0.1-q0.9", ylab = "Frequency", main = samples[i])
    }
  }
  dev.off()
})

################################################################################################
setGeneric("plotShape",function(object,...)standardGeneric("plotShape"))
setMethod("plotShape","TSSr", function(object
                                       ,samples = "all"
){
  message("Plotting Shape graphs...")
  TCs <- object@clusterShape
  sampleLabels <- object@sampleLabelsMerged
  if(samples == "all"){
    tc <- TCs
    pdf(file = paste("Shape_plot_of_ALL_samples.pdf", sep = "")
        ,width = 8, height = 8, onefile = T, bg = "transparent", family = "Helvetica", fonts = NULL)
    for(i in 1:length(sampleLabels)){
      temp <- tc[[sampleLabels[i]]]
      hist(temp$shape.score, breaks = 40, col = rainbow(length(sampleLabels))[i]
           , xlab = "shape score", ylab = "Frequency", main = sampleLabels[i])
    }
  }else{
    tc <- TCs[[samples]]
    pdf(file = paste("Shape_plot_of_", paste(samples, collapse = "_"), "_samples.pdf", sep = "")
        ,width = 8, height = 8, onefile = T, bg = "transparent", family = "Helvetica", fonts = NULL)
    for(i in 1:length(samples)){
      temp <- tc[[samples[i]]]
      hist(temp$shape.score, breaks = 40, col = rainbow(length(samples))[i]
           , xlab = "shape score", ylab = "Frequency", main = samples[i])
    }
  }
  dev.off()
})

################################################################################################
setGeneric("plotDE",function(object,...)standardGeneric("plotDE"))
setMethod("plotDE","TSSr", function(object
                                    ,withGeneName = "TRUE"
                                    ,xlim=c(-2.5, 2.5)
                                    ,ylim=c(0,10)
){
  message("Plotting DE graphs...")
  pdf(file = paste("Volcano_plot.pdf", sep = ""),width = 8, height = 8,bg = "transparent"
      , family = "Helvetica", fonts = NULL)
  D.names <- names(object@DEtables)
  for(i in 1:length(D.names)){
    res <- object@DEtables[[D.names[i]]]$DEtable
    plot(res$log2FoldChange,-log10(res$pvalue), pch = 20, xlim = xlim, ylim = ylim
         ,main = "Volcano plot", xlab = "log2FoldChange", ylab = "-log10(pvalue)")
    # Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
    with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
    with(subset(res, abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="orange"))
    with(subset(res, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="green"))
    if(withGeneName == "TRUE"){
      with(subset(res, padj<.05 & abs(log2FoldChange)>1), textxy(log2FoldChange, -log10(pvalue), labs=gene, cex=.8))
    }
    
  }
  ##
  dev.off()
})

################################################################################################
setGeneric("plotTSS",function(object,...)standardGeneric("plotTSS"))
setMethod("plotTSS","TSSr", function(object
                                     ,samples
                                     ,genelist
                                     ,up.dis =500
                                     ,down.dis = 500
){
  message("Plotting TSS graphs...")
  ##initialize data
  tss.filtered <- object@TSSfilteredMatrix
  clusters <- object@consensusClusters[[1]]
  tss.raw <- object@TSSfilteredMatrix
  refGFF <- object@refSource
  organismName <- object@organismName
  sampleLabelsMerged <- object@sampleLabelsMerged
  
  ##prepare tss table
  if(FALSE %in% unique(samples %in% sampleLabelsMerged)){
    stop("No data for one or more given samples! ")
  }else{
    cols <- c("chr","pos","strand", samples)
    tss <- tss.filtered[,.SD, .SDcols = cols]
  }
  tss.p <- tss[tss$strand == "+",]
  tss.m <- tss[tss$strand == "-",]
  tss.m[,(samples) := lapply(.SD, "*",-1), .SDcols = samples]
  tss <- rbind(tss.p,tss.m)
  
  ##prepare annotation file
  txdb <- suppressWarnings(makeTxDbFromGFF(refGFF, organismName, format = "auto"))
  ref <- setDT(as.data.frame(genes(txdb)))
  ref <- ref[gene_id %in% genelist,]
  
  pdf(file = paste("TSS_graphs.pdf", sep = "")
      ,width = 10, height = 8, onefile = T, bg = "transparent", family = "Helvetica", fonts = NULL)
  for (i in 1:nrow(ref)){
    df <- ref[i,]
    .plotTSS(tss, clusters,df, samples, up.dis, down.dis)
  }
  dev.off()
})

################################################################################################
setGeneric("exportTSStable",function(object,...)standardGeneric("exportTSStable"))
setMethod("exportTSStable","TSSr", function(object
                                            ,data = "raw"
                                            ,merged = "TRUE"
){
  message("Exporting TSS table...")
  if(data == "raw"){
    if(merged == "TRUE"){
      tss <- object@TSSmergedMatrix
      write.table(tss, file = paste("ALL.samples.TSS",data,"txt", sep = "."), sep = "\t", quote = F, row.names = F)
    }else{
      tss <- object@TSSrawMatrix
      write.table(tss, file = paste("ALL.samples.TSS",data,"txt", sep = "."), sep = "\t", quote = F, row.names = F)
    }
  }else if(data == "filtered"){
    tss <- object@TSSfilteredMatrix
    write.table(tss, file = paste("ALL.samples.TSS",data,"txt", sep = "."), sep = "\t", quote = F, row.names = F)
  }else{
    stop("No data for the given TSS data type!")
  }
})

################################################################################################
setGeneric("exportTagClustersTable",function(object,...)standardGeneric("exportTagClustersTable"))
setMethod("exportTagClustersTable","TSSr", function(object
                                                    ,data = "assigned"
){
  message("Exporting tagClusters table...")
  if(data == "assigned"){
    tc <- object@assignedClusters
    samples <- object@sampleLabelsMerged
    for(i in 1:length(samples)){
      temp <- tc[[samples[i]]]
      write.table(temp, file = paste(samples[i],"tagClusters.assigned","txt", sep = "."), sep = "\t", quote = F, row.names = F)
    }
  }else if(data == "unassigned"){
    tc <- object@unassignedClusters
    samples <- object@sampleLabelsMerged
    for(i in 1:length(samples)){
      temp <- tc[[samples[i]]]
      write.table(temp, file = paste(samples[i],"tagClusters.unassigned","txt", sep = "."), sep = "\t", quote = F, row.names = F)
    }
  }else{
    stop("No data for the given tag cluster data type!")
  }
})

################################################################################################
setGeneric("exportShapeTable",function(object,...)standardGeneric("exportShapeTable"))
setMethod("exportShapeTable","TSSr", function(object
){
  message("Exporting promoter shape table...")
  s <- object@clusterShape
  if(!is.null(s)){
    samples <- object@sampleLabelsMerged
    for(i in 1:length(samples)){
      temp <- s[[samples[i]]]
      write.table(temp, file = paste(samples[i],"promoter.shape","txt", sep = "."), sep = "\t", quote = F, row.names = F)
    }
  }else{
    stop("No data for the promoter shape!")
  }
})

################################################################################################
setGeneric("exportDETable",function(object,...)standardGeneric("exportDETable"))
setMethod("exportDETable","TSSr", function(object
                                           ,data = "all"
){
  message("Exporting differential expression table...")
  D.names <- names(object@DEtables)
  if(data == "all"){
    for(i in 1:length(D.names)){
      temp <- object@DEtables[[D.names[i]]]$DEtable
      write.table(temp, file = paste(D.names[i],"DE.table.aLL.txt", sep = "."), sep = "\t", quote = F, row.names = F)
    }
  }else if(data == "sig"){
    for(i in 1:length(D.names)){
      temp <- object@DEtables[[D.names[i]]]$DEsig
      write.table(temp, file = paste(D.names[i],"DE.table.sig.txt", sep = "."), sep = "\t", quote = F, row.names = F)
    }
  }else{
    stop("No data for the differential expression!")
  }
})

################################################################################################
setGeneric("exportShiftTable",function(object,...)standardGeneric("exportShiftTable"))
setMethod("exportShiftTable","TSSr", function(object
){
  message("Exporting promoter shift table...")
  D.names <- names(object@PromoterShift)
  for(i in 1:length(D.names)){
    temp <- object@PromoterShift[[D.names[i]]]
    write.table(temp, file = paste(D.names[i],"promoter.shift.table.txt", sep = "."), sep = "\t", quote = F, row.names = F)
  }
})

################################################################################################
