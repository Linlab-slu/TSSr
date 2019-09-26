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
library(DESeq2)
library(BSgenome.Scerevisiae.UCSC.sacCer3)
library(calibrate)
library(ggfortify)
################################################################################################
inputFiles <- c("S01.sorted.bam","S02.sorted.bam","S03.sorted.bam","S04.sorted.bam")
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
myTSSr <- clusterTSS(myTSSr, data = "filtered", method = "peakclu",clusterThreshold = 1, useMultiCore=TRUE, numCores = NULL)
myTSSr <- consensusCluster(myTSSr, data = "filtered", dis = 50)
myTSSr <- shapeCluster(myTSSr,data = "consensusClusters" , method = "PSS",useMultiCore= TRUE, numCores = NULL)
myTSSr <- annotateCluster(myTSSr,clusters = "consensusClusters",annotationType = "genes",upstream=1000,upstreamOverlap = 500,downstream = 0)
myTSSr <- deGene(myTSSr,comparePairs=list(c("control","treat")), pval = 0.01)
myTSSr <- shiftPromoter(myTSSr,comparePairs=list(c("control","treat")), pval = 0.01)
##plot
plotCorrelation(myTSSr, samples = "all")
plotPCA(myTSSr)
plotInterQuantile(myTSSr)
plotShape(myTSSr)
plotDE(myTSSr, withGeneName = "TRUE")
plotTSS(myTSSr,samples=c("control","treat"),genelist=c("YBL017C","YBL067C"),up.dis =500,down.dis = 500)
##exporting data table
exportTSStable(myTSSr, data = "filtered")
exportTagClustersTable(myTSSr, data = "assigned")
exportShapeTable(myTSSr)
exportDETable(myTSSr, data = "sig")
exportShiftTable(myTSSr)
##export to bedGraph/BigWig
exportTSStoBedgraph(myTSSr, data = "filtered", format = "bedGraph")
exportTSStoBedgraph(myTSSr, data = "filtered", format = "BigWig")
exportClustersToBed(myTSSr, data = "consensusClusters")

##check run time
ptm <- proc.time()
myTSSr <- clusterTSS(myTSSr, data = "filtered", method = "peakclu",clusterThreshold = 1, useMultiCore=TRUE, numCores = NULL)
proc.time() - ptm
