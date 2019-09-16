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
                                    ,ylim=c()
){
  message("Plotting DE graphs...")
  res <- object@DEtables$DEtable
  pdf(file = paste("Volcano_plot.pdf", sep = ""),width = 8, height = 8,bg = "transparent"
      , family = "Helvetica", fonts = NULL)
  ##
  if(is.null(xlim)){
    with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot"))
  }else{
    with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=xlim,ylim = ylim))
  }
  
  # Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
  with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
  with(subset(res, abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="orange"))
  with(subset(res, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="green"))
  if(withGeneName == "TRUE"){
    with(subset(res, padj<.05 & abs(log2FoldChange)>1), textxy(log2FoldChange, -log10(pvalue), labs=gene, cex=.8))
  }
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
  clusters <- object@tagClusters[[1]]
  tss.raw <- object@TSSfilteredMatrix
  refGFF <- object@refSource
  organismName <- object@organismName
  
  ##prepare tss table
  if(FALSE %in% unique(samples %in% sampleLabels)){
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
      write.table(tss, file = "TSS_table_ALL_samples_merged.txt", sep = "\t", quote = F, row.nams = F)
    }else{
      tss <- object@TSSrawMatrix
      write.table(tss, file = "TSS_table_ALL_samples_unmerged.txt", sep = "\t", quote = F, row.nams = F)
    }
  }else if(data == "filtered"){
    tss <- object@TSSfilteredMatrix
    write.table(tss, file = "TSS_table_ALL_samples_filtered.txt", sep = "\t", quote = F, row.nams = F)
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
      write.table(temp, file = paste("TagClusters_assigned_",samples[i],".txt", sep = ""), sep = "\t", quote = F, row.nams = F)
    }
  }else if(data == "unassigned"){
    tc <- object@unassignedClusters
    samples <- object@sampleLabelsMerged
    for(i in 1:length(samples)){
      temp <- tc[[samples[i]]]
      write.table(temp, file = paste("TagClusters_unassigned_",samples[i],".txt", sep = ""), sep = "\t", quote = F, row.nams = F)
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
      write.table(temp, file = paste("Peomoter_shape_",samples[i],".txt", sep = ""), sep = "\t", quote = F, row.nams = F)
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
  if(data == "all"){
    temp <- object@DEtables$DEtable
    write.table(temp, file = paste("DE_table_ALL.txt", sep = ""), sep = "\t", quote = F, row.nams = F)
  }else if(data == "sig"){
    temp <- object@DEtables$DEsig
    write.table(temp, file = paste("DE_table_sig.txt", sep = ""), sep = "\t", quote = F, row.nams = F)
  }else{
    stop("No data for the differential expression!")
  }
})

################################################################################################
