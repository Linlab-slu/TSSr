################################################################################################
#' Pairwise scatter plots and correlations of TSS signal
#'
#' @description Calculates the pairwise correlation coefficients between samples and
#' creates a matix showing pairwise scatter plots and correlation coefficients.
#'
#' @usage plotCorrelation(object, samples = "all")
#' @param object A TSSr object.
#' @param samples Specify samples to be plotted. Can be either "all" to plot all samples in the object
#' or a subset of samples in the object. Default is "all".
#'
#' @return
#' @export
#'
#' @examples
#' \donttest{
#' plotCorrelation(exampleTSSr, samples = "all")
#' }
#'
setGeneric("plotCorrelation",function(object, samples = "all")standardGeneric("plotCorrelation"))
#' @rdname plotCorrelation
#' @export
setMethod("plotCorrelation",signature(object = "TSSr"), function(object, samples){
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
#' Plotting principle component analysis (PCA)
#'
#' @description Calculates principle component analysis (PCA) of all samples and creates a biplot
#' which includes the position of each sample in terms of PC1 and PC2.
#'
#' @usage plotTssPCA(object, TSS.threshold =10)
#'
#' @param object A TSSr object.
#' @param TSS.threshold Only TSSs with raw signal >= TSS.threshold will be included in PCA analysis
#'
#' @return
#' @export
#'
#' @examples
#' \donttest{
#' plotTssPCA(exampleTSSr)
#' }
setGeneric("plotTssPCA",function(object, TSS.threshold =10)standardGeneric("plotTssPCA"))
#' @rdname plotTssPCA
#' @export
setMethod("plotTssPCA",signature(object = "TSSr"), function(object, TSS.threshold){
  message("Plotting TSS PCA...")
  tss <- object@TSSrawMatrix
  tss <- tss[,4:ncol(tss)]
  tss <- tss[apply(tss >= TSS.threshold, 1, any),]
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
#' Plot core promoter interquantile width
#'
#' @description Plots histograms of the interquantile width of processed clusters.
#' @usage plotInterQuantile(object, samples ="all", tagsThreshold = 1)
#'
#' @param object A TSSr object.
#' @param samples Specify samples to be plotted. Default is "all".
#' @param tagsThreshold Excludes clusters with tags < tagsThreshold.
#'
#' @return
#' @export
#'
#' @examples
#' \donttest{
#' plotInterQuantile(exampleTSSr, samples = "all")
#' }
setGeneric("plotInterQuantile",function(object, samples = "all", tagsThreshold = 1)standardGeneric("plotInterQuantile"))
#' @rdname plotInterQuantile
#' @export
setMethod("plotInterQuantile",signature(object = "TSSr"), function(object, samples, tagsThreshold){
  message("Plotting interquantile graphs...")
  TCs <- object@clusterShape
  sampleLabels <- object@sampleLabelsMerged
  ##define variable as a NULL value
  interquantile_width = NULL

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
#' Plot core promoter shape
#'
#' @description Plots histograms of core promoter shape scores.
#'
#' @usage plotShape(object, samples = "all")
#'
#' @param object A TSSr object.
#' @param samples Specify samples to be plotted. Default is "all".
#'
#' @return
#' @export
#'
#' @examples
#' \donttest{
#' plotShape(exampleTSSr)
#' }
setGeneric("plotShape",function(object, samples = "all")standardGeneric("plotShape"))
#' @rdname plotShape
#' @export
setMethod("plotShape",signature(object = "TSSr"), function(object ,samples){
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
#' Plot gene differential expressions
#'
#' @description Vocano plots of gene differential expression (with DESeq2 method) results.
#'
#' @usage plotDE(object, withGeneName = "TRUE",xlim, ylim)
#'
#' @param object A TSSr object.
#' @param withGeneName Specify whether to display names for genes which are differentially expressed. Default is "TRUE".
#' @param xlim Only enes of which log2FoldChange value within the xlim range are plotted. Default xlim = c(-2.5, 2.5).
#' @param ylim Only genes of which -log10(pvalue) within the ylim range are plotted. Default ylim = c(0, 10).
#'
#' @return
#' @export
#'
#' @examples
#' \donttest{
#' plotDE(exampleTSSr, withGeneName = "TRUE")
#' plotDE(exampleTSSr, withGeneName = "FALSE")
#' }
setGeneric("plotDE",function(object
                             ,withGeneName = "TRUE"
                             ,xlim=c(-2.5, 2.5)
                             ,ylim=c(0,10))standardGeneric("plotDE"))
#' @rdname plotDE
#' @export
setMethod("plotDE",signature(object = "TSSr"), function(object, withGeneName, xlim, ylim){
  message("Plotting DE graphs...")
  pdf(file = paste("Volcano_plot.pdf", sep = ""),width = 8, height = 8,bg = "transparent"
      , family = "Helvetica", fonts = NULL)
  D.names <- names(object@DEtables)
  ##define variable as a NULL value
  padj = log2FoldChange = NULL

  for(i in 1:length(D.names)){
    res <- object@DEtables[[D.names[i]]]$DEtable
    plot(res$log2FoldChange,-log10(res$pvalue), pch = 20, xlim = xlim, ylim = ylim
         ,main = D.names[i], xlab = "log2FoldChange", ylab = "-log10(pvalue)")
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
#' Plot TSSs and clusters
#'
#' @description Plots Gviz-track of TSSs, clusters, and genes.
#' @usage plotTSS(object,samples,tssData = "processed",clusters = "filtered",
#' clusterThreshold = 0.02,genelist,Bidirection = TRUE,up.dis =500,down.dis = 500)
#'
#' @param object A TSSr object.
#' @param samples Specify samples to be included for plotting.
#' @param tssData Specify which TSS data to be included for plotting: "raw" or "process".
#' @param clusters Specify which cluster data to be included for plotting: "all", "assigned", or "processed".
#' @param clusterThreshold Ignore downstream clusters if signal < filterClusterThreshold*the strongest
#' clusters within the same gene promoter region. Default value = 0.02.
#' @param genelist List of gene names used for plotting.
#' @param Bidirection Specify whether to display bidirectional TSS signals within defined region. Default is TRUE.
#' @param up.dis Distance upstream of genes to specify plotting range. Default value = 500.
#' @param down.dis Distance downstream of genes to specify plotting range. Default value = 500.
#'
#' @return
#' @export
#'
#' @examples
#' \donttest{
#' plotTSS(exampleTSSr, samples=c("control","treat"), genelist=c("YBL017C","YBL067C")
#' ,up.dis =500, down.dis = 500)
#' }
setGeneric("plotTSS",function(object,samples
                              ,tssData = "processed"
                              ,clusters = "filtered"
                              ,clusterThreshold = 0.02
                              ,genelist
                              ,Bidirection= TRUE
                              ,up.dis =500
                              ,down.dis = 500)standardGeneric("plotTSS"))
#' @rdname plotTSS
#' @export
setMethod("plotTSS",signature(object = "TSSr"), function(object, samples, tssData, clusters, clusterThreshold
                                                         , genelist, Bidirection, up.dis, down.dis){
  message("Plotting TSS graphs...")
  ##define variable as a NULL value
  gene_id = NULL

  ##initialize data
  if(clusters == "all"){
    cs <- object@consensusClusters
  }else if(clusters == "assigned"){
    cs <- object@assignedClusters
  }else if(clusters == "filtered"){
    cs <- object@filteredClusters
  }else{
    stop("No cluster data for the given clusters option! ")
  }
  if(tssData == "processed"){
    tss <- object@TSSprocessedMatrix
  }else if(tssData == "raw"){
    tss <- object@TSSmergedMatrix
  }
  refGFF <- object@refSource
  organismName <- object@organismName
  sampleLabelsMerged <- object@sampleLabelsMerged

  ##prepare tss table
  if(FALSE %in% unique(samples %in% sampleLabelsMerged)){
    stop("No data for one or more given samples! ")
  }else{
    cols <- c("chr","pos","strand", samples)
    tss <- tss[,.SD, .SDcols = cols]
    cs <- cs[samples]
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
    .plotTSS(tss, cs,df, samples, Bidirection, up.dis, down.dis)
  }
  dev.off()
})

################################################################################################
#' Export TSS tables
#'
#' @description Exports TSS tables to text file.
#' @usage exportTSStable(object, data = "raw", merged = "TRUE")
#'
#' @param object A TSSr object.
#' @param data Specify which data will be exported: "raw" or "processed". Default is "raw".
#' @param merged Specify whether to export merged TSS table. Used only if data = "raw".
#' @return
#' @export
#'
#' @examples
#' \donttest{
#' exportTSStable(exampleTSSr)
#' exportTSStable(exampleTSSr, data="raw")
#' }
setGeneric("exportTSStable",function(object, data = "raw", merged = "TRUE")standardGeneric("exportTSStable"))
#' @rdname exportTSStable
#' @export
setMethod("exportTSStable",signature(object = "TSSr"), function(object, data, merged){
  message("Exporting TSS table...")
  if(data == "raw"){
    if(merged == "TRUE"){
      tss <- object@TSSprocessedMatrix
      write.table(tss, file = paste("ALL.samples.TSS",data,"txt", sep = "."), sep = "\t", quote = F, row.names = F)
    }else{
      tss <- object@TSSrawMatrix
      write.table(tss, file = paste("ALL.samples.TSS",data,"txt", sep = "."), sep = "\t", quote = F, row.names = F)
    }
  }else if(data == "processed"){
    tss <- object@TSSprocessedMatrix
    write.table(tss, file = paste("ALL.samples.TSS",data,"txt", sep = "."), sep = "\t", quote = F, row.names = F)
  }else{
    stop("No data for the given TSS data type!")
  }
})

################################################################################################
#' Export cluster tables
#'
#' @description Export cluster tables to text files.
#' @usage exportClustersTable(object, data = "filtered")
#'
#' @param object A TSSr object.
#' @param data Specify which cluster data will be exported: "tagClusters", "consensusClusters",
#' "assigned", "unassigned", "filtered". Default is "filtered".
#'
#' @return
#' @export
#'
#' @examples
#' \donttest{
#' 	exportClustersTable(exampleTSSr, data = "tagClusters")
#' 	exportClustersTable(exampleTSSr, data = "consensusClusters")
#' 	exportClustersTable(exampleTSSr, data = "assigned")
#' 	exportClustersTable(exampleTSSr, data = "unassigned")
#' 	exportClustersTable(exampleTSSr, data = "filtered")
#' }
setGeneric("exportClustersTable",function(object, data = "filtered")standardGeneric("exportClustersTable"))
#' @rdname exportClustersTable
#' @export
setMethod("exportClustersTable",signature(object = "TSSr"), function(object, data){
  if(data == "tagClusters"){
    message("Exporting tagClusters table...")
    tc <- object@tagClusters
    samples <- object@sampleLabelsMerged
    for(i in 1:length(samples)){
      temp <- tc[[samples[i]]]
      write.table(temp, file = paste(samples[i],"tagClusters","txt", sep = "."), sep = "\t", quote = F, row.names = F)
    }
  }else if(data == "consensusClusters"){
    message("Exporting consensusClusters table...")
    tc <- object@consensusClusters
    samples <- object@sampleLabelsMerged
    for(i in 1:length(samples)){
      temp <- tc[[samples[i]]]
      write.table(temp, file = paste(samples[i],"consensusClusters","txt", sep = "."), sep = "\t", quote = F, row.names = F)
    }
  }else if(data == "assigned"){
    message("Exporting assignedClusters table...")
    tc <- object@assignedClusters
    samples <- object@sampleLabelsMerged
    for(i in 1:length(samples)){
      temp <- tc[[samples[i]]]
      write.table(temp, file = paste(samples[i],"assignedClusters","txt", sep = "."), sep = "\t", quote = F, row.names = F)
    }
  }else if(data == "unassigned"){
    message("Exporting unassignedClusters table...")
    tc <- object@unassignedClusters
    samples <- object@sampleLabelsMerged
    for(i in 1:length(samples)){
      temp <- tc[[samples[i]]]
      write.table(temp, file = paste(samples[i],"unassignedClusters","txt", sep = "."), sep = "\t", quote = F, row.names = F)
    }
  }else if(data == "filtered"){
    message("Exporting filteredClusters table...")
    tc <- object@filteredClusters
    samples <- object@sampleLabelsMerged
    for(i in 1:length(samples)){
      temp <- tc[[samples[i]]]
      write.table(temp, file = paste(samples[i],"filteredClusters","txt", sep = "."), sep = "\t", quote = F, row.names = F)
    }
  }else{
    stop("No data for the given tag cluster data type!")
  }
})

################################################################################################
#' Export core promoter shape score tables
#'
#' @description Exports core promoter shape score tables to text files. Shape score is calculated with
#' shapeCluster(object) method.
#' @usage exportShapeTable(object)
#'
#' @param object A TSSr object.
#'
#' @return
#' @export
#'
#' @examples
#' \donttest{
#' exportShapeTable(exampleTSSr)
#' }
setGeneric("exportShapeTable",function(object)standardGeneric("exportShapeTable"))
#' @rdname exportShapeTable
#' @export
setMethod("exportShapeTable",signature(object = "TSSr"), function(object
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
#' Export gene differential expression results table
#'
#' @description Exports gene differential expression results table to text files.
#' @usage exportDETable(object, data = "sig")
#'
#' @param object A TSSr object.
#' @param data Specify which data will be exported: "all" or "sig".
#'
#' @return
#' @export
#'
#' @examples
#' \donttest{
#' exportDETable(exampleTSSr, data="sig")
#' }
setGeneric("exportDETable",function(object, data = "sig")standardGeneric("exportDETable"))
#' @rdname exportDETable
#' @export
setMethod("exportDETable",signature(object = "TSSr"), function(object, data){
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
#' Export core promoter shift table
#'
#' @description Export core promoter shift tables to text files.
#' @usage exportShiftTable(object)
#'
#' @param object A TSSr object.
#'
#' @return
#' @export
#'
#' @examples
#' \donttest{
#' exportShiftTable(exampleTSSr)
#' }
setGeneric("exportShiftTable",function(object)standardGeneric("exportShiftTable"))
#' @rdname exportShiftTable
#' @export
setMethod("exportShiftTable",signature(object = "TSSr"), function(object
){
  message("Exporting promoter shift table...")
  D.names <- names(object@PromoterShift)
  for(i in 1:length(D.names)){
    temp <- object@PromoterShift[[D.names[i]]]
    write.table(temp, file = paste(D.names[i],"promoter.shift.table.txt", sep = "."), sep = "\t", quote = F, row.names = F)
  }
})

################################################################################################
#' Creating Bedgraph/BigWig tracks of TSSs
#' @description Creates bedGraph/BigWig files of TSSs that can be visualized in the UCSC Genome Browser
#'  and Integrative Genomics Viewer (IGV).
#'
#' @usage exportTSStoBedgraph(object, data = "processed", format = "bedGraph",oneFile = FALSE)
#'
#' @param object A TSSr object.
#' @param data Specify which data will be exported: "raw" or "processed". Default is "processed".
#' @param format The format of output files: "bedGraph" or "BigWig". Default is "bedGraph".
#' @param oneFile Logical, specify whether to export individual TSS tracks into the one bedGraph
#' file (TRUE) of in separate bedGraph files (FALSE).
#' @return
#' @export
#'
#' @examples
#' \donttest{
#' exportTSStoBedgraph(exampleTSSr, data = "processed", format = "bedGraph")
#' }
setGeneric("exportTSStoBedgraph",function(object,data = "processed"
                                          ,format = "bedGraph"
                                          ,oneFile = FALSE)standardGeneric("exportTSStoBedgraph"))
#' @rdname exportTSStoBedgraph
#' @export
setMethod("exportTSStoBedgraph",signature(object = "TSSr"), function(object, data, format, oneFile){
  Genome <- .getGenome(object@genomeName)
  sampleLabelsMerged <- object@sampleLabelsMerged
  ##define variable as a NULL value
  score = strand = NULL

  if(data == "processed"){
    tss.dt <- object@TSSprocessedMatrix
  }else{tss.dt <- object@TSSrawMatrix}
  for (i in 1:length(sampleLabelsMerged)){
    temp <- tss.dt[,.SD, .SDcols = c("chr","pos","strand",sampleLabelsMerged[i])]
    setnames(temp, colnames(temp)[[4]], "score")
    temp <- temp[score >0,]
    if(oneFile == TRUE){
      message("Exporting TSS to bedgraph...")
      temp[, score := ifelse(strand == "+", score, score*(-1))]
      temp <- makeGRangesFromDataFrame(temp, start.field = "pos", end.field = "pos", keep.extra.columns = TRUE)
      export(temp,paste(sampleLabelsMerged[i], "TSS", data, "bedGraph", sep = "."), format = "bedGraph")
    }else{
      temp.p <- temp[strand == "+",]
      temp.m <- temp[strand == "-",]
      temp.m[, score := score*(-1)]
      temp.p <- makeGRangesFromDataFrame(temp.p, start.field = "pos", end.field = "pos", keep.extra.columns = TRUE)
      temp.m <- makeGRangesFromDataFrame(temp.m, start.field = "pos", end.field = "pos", keep.extra.columns = TRUE)
      if(format == "bedGraph"){
        message("Exporting TSS to bedgraph...")
        export(temp.p,paste(sampleLabelsMerged[i], "TSS", data, "plus.bedGraph", sep = "."), format = "bedGraph")
        export(temp.m,paste(sampleLabelsMerged[i], "TSS", data, "minus.bedGraph", sep = "."), format = "bedGraph")
      }else if(format == "BigWig"){
        message("Exporting TSS to BigWig...")
        seqlengths(temp.p) <- seqlengths(Genome)[seqnames(Genome) %in% as.character(temp.p@seqnames)]
        seqlengths(temp.m) <- seqlengths(Genome)[seqnames(Genome) %in% as.character(temp.m@seqnames)]
        export(temp.p,paste(sampleLabelsMerged[i], "TSS", data, "plus.BigWig", sep = "."), format = "BigWig")
        export(temp.m,paste(sampleLabelsMerged[i], "TSS", data, "minus.BigWig", sep = "."), format = "BigWig")
      }
    }
  }
})
################################################################################################
#' Creating bed files of clusters
#'
#' @description Creates bed files of clusters.
#' @usage exportClustersToBed(object, data = "consensusClusters", filtered = TRUE)
#'
#' @param object A TSSr object.
#' @param data Specify which data will be exported: "tagClusters" or "consensusClusters". Default is "consensusClusters".
#' @param filtered Specify which consensus clusters will be exported. Used only if data = "consensusClusters. Default is TRUE.
#'
#' @return
#' @export
#'
#' @examples
#' \donttest{
#' exportTSStoBedgraph(exampleTSSr, data = "tagClusters")
#' exportTSStoBedgraph(exampleTSSr, dataexampleTSSr = "consensusClusters")
#' }
setGeneric("exportClustersToBed",function(object,data = "consensusClusters", filtered = TRUE)
  standardGeneric("exportClustersToBed"))
#' @rdname exportClustersToBed
#' @export
setMethod("exportClustersToBed",signature(object = "TSSr"), function(object, data, filtered){
  message("Exporting clusters to bed...")
  sampleLabelsMerged <- object@sampleLabelsMerged
  if(data == "tagClusters"){
    cs <- object@tagClusters
  }else if(data == "consensusClusters"){
    if(filtered == TRUE){
      cs <- object@filteredClusters
    }else{
      cs <- object@consensusClusters
    }
  }
  for (i in 1:length(sampleLabelsMerged)){
    temp <- cs[[sampleLabelsMerged[i]]]
    temp <- .getBed(temp)
    df <- file(paste(sampleLabelsMerged[i],data,"bed", sep = "."), open = "wt")
    writeLines(paste('track name="',sampleLabelsMerged[i]
                     ,"(",data
                     ,'(TC) q0.1-q0.9)" description="'
                     ,sampleLabelsMerged[i],"(",data
                     ,'(TC) q0.1-q0.9)" '
                     ,'visibility="pack" color=0,255,255', sep = ""), df)
    write.table(temp, df, sep = "\t", quote = F, row.names = F, col.names = F)
    close(df)
  }
})

################################################################################################
