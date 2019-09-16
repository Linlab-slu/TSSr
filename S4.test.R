################################################################################################
##S4 test
################################################################################################
setwd("C:/Users/zhaolianlu/Google Drive/CAGE/TSSr_package/data")
library(Rsamtools)
library(stringr)
library(data.table)
library(rtracklayer)
library(GenomicRanges)
library(GenomicFeatures)
library(GenomicAlignments)
library(Gviz)
library(Deseq2)
library(BSgenome.Scerevisiae.UCSC.sacCer3)
################################################################################################
inputFiles <- list.files()[1:4]
myTSSr <- new("TSSr", genomeName = "BSgenome.Scerevisiae.UCSC.sacCer3"
              ,inputFiles = inputFiles
              ,inputFilesType= "bam"
              ,sampleLabels = c("SL01","SL02","SL03","SL04")
              ,sampleLabelsMerged = c("control","treat")
              ,refSource = "saccharomyces_cerevisiae.SGD.gff"
              ,organismName = "saccharomyces cerevisiae"
)
myTSSr <- getTSS(myTSSr)
myTSSr <- mergeSamples(myTSSr, mergeIndex = c(1,1,2,2))
myTSSr <- normalizeTSS(myTSSr)
myTSSr <- filterTSS(myTSSr, method = "TPM")
myTSSr <- clusterTSS(myTSSr, data = "filtered", method = "peakclu",clusterThreshold = 1, useMultiCore=FALSE)
myTSSr <- shapeCluster(myTSSr, method = "PSS",useMultiCore= FALSE, numCores = NULL)
myTSSr <- annotateCluster(myTSSr,clusters = "tagClusters",annotationType = "genes",upstream=1000,upstreamOverlap = 500,downstream = 0)
myTSSr <- deGene(myTSSr,"control","treat", pval = 0.01)
myTSSr <- shiftPromoter(myTSSr,"control","treat", pval = 0.1) ##must be based on consensus, not done yet, Sep16
##plot
plotCorrelation(myTSSr, samples = "all")
plotInterQuantile(myTSSr)
plotShape(myTSSr)
plotDE(myTSSr, withGeneName = "TRUE")
plotTSS(myTSSr,samples=c("control","treat"),genelist=c("YBR298C"),up.dis =500,down.dis = 500)
##exporting data table
exportTSStable(myTSSr, data = "filtered")
exportTagClustersTable(myTSSr, data = "assigned")
exportShapeTable(myTSSr)
exportDETable(myTSSr, data = "sig")
