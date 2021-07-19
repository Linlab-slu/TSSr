################################################################################################
#' Analysis of gene differential expression.
#'
#' @description Analyzes gene-level differential expression using DESeq2 method (Love et al., 2014).
#' @usage deGene(object,comparePairs=list(c("control","treat")), pval = 0.01,
#'  useMultiCore=FALSE, numCores = NULL)
#'
#' @param object A TSSr object.
#' @param comparePairs Specified list of sample pairs for comparison with DESeq2 method.
#' @param pval Genes with adjusted p value >= pVal will be returned. Default value = 0.01.
#' @param useMultiCore Logical indicating whether multiple cores are used (TRUE) or not (FALSE). Default is FALSE.
#' @param numCores Number of cores are used in clustering step. Used only if useMultiCore = TRUE. Default is NULL.
#'
#'
#' @export
#'
#' @examples
#' \donttest{
#' deGene(exampleTSSr,comparePairs=list(c("control","treat")), pval = 0.01)
#' }
setGeneric("deGene",function(object, comparePairs=list(c("control","treat")), pval=0.01)standardGeneric("deGene"))
#' @rdname deGene
#' @export
setMethod("deGene",signature(object = "TSSr"), function(object, comparePairs, pval){
  ##initialize data
  message("\nCalculating gene differential expression...")
  objName <- deparse(substitute(object))
  sampleLabels <- object@sampleLabels
  sampleLabelsMerged <- object@sampleLabelsMerged
  gene_raw_count <- object@TAGtables
  ##define variable as a NULL value
  padj = NULL
  D_updated <- lapply(as.list(seq(comparePairs)), function(i){
  sampleOne <- comparePairs[[i]][1]
  samplex1 <- sampleLabels[which(mergeIndex ==which(sampleLabelsMerged == sampleOne))]
  sampleTwo <- comparePairs[[i]][2]
  sampley1 <- sampleLabels[which(mergeIndex ==which(sampleLabelsMerged == sampleTwo))]
  #different with expression by Deseq
  ##get raw count tables
  one <- gene_raw_count[[sampleOne]]
  two <- gene_raw_count[[sampleTwo]]
  ##merge the two raw count tables together by genes
  one[,2:ncol(one)] <- sapply(one[,2:ncol(one)], as.integer)
  two[,2:ncol(two)] <- sapply(two[,2:ncol(two)], as.integer)
  #setnames(one, colnames(one), c("gene",samplex))
  #setnames(two, colnames(two), c("gene",sampley))
  Dtable <- merge(one, two, by = c("gene"), all = T)
  Dtable[is.na(Dtable)] = 0
  rownames(Dtable) <- Dtable[,1]
  Dtable1 <- Dtable[,-1]
  Dtable1 <- data.matrix(Dtable1)
  #default two replicates
  condition <- factor(c(rep(sampleOne, times = length(samplex1)), rep(sampleTwo, times = length(sampley1))))
  dds <- DESeqDataSetFromMatrix(countData = Dtable1,data.frame(condition), ~ condition)
  dds$condition <- factor(dds$condition, levels = c(sampleOne, sampleTwo))
  dds <- DESeq(dds)
  res <- results(dds)
  res1 <- as.data.frame(res)
  res1$gene <- rownames(res1)
  res2 <- merge(Dtable,res1,by=c("gene"),all=T)
  res2 <- res2[order(res2$padj),]
  DE.dt <- as.data.frame(res2)
  DE.sig <- subset(DE.dt, padj < pval)
  #DE.dt <- DE.dt[,c(ncol(DE.dt), 1:(ncol(DE.dt)-1))]
  #DE.sig <- DE.sig[,c(ncol(DE.sig), 1:(ncol(DE.sig)-1))]
  setDT(DE.dt)
  setDT(DE.sig)
  DE <- list("DEtable" = DE.dt, "DEsig" = DE.sig)
  return(DE)
})
D.names <- sapply(as.list(seq(comparePairs)), function(i){
  paste0(comparePairs[[i]][1],"_VS_",comparePairs[[i]][2], sep ="")
})
names(D_updated) <- D.names
return(D_updated)
object@DEtables <- D_updated
assign(objName, object, envir = parent.frame())
})
