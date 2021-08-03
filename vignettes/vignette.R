## ----eval=TRUE----------------------------------------------------------------
library(TSSr)

## ----eval=TRUE----------------------------------------------------------------
inputFiles <- c("S01.sorted.bam", "S02.sorted.bam", "S03.sorted.bam", "S04.sorted.bam")
myTSSr <- new("TSSr", genomeName = "BSgenome.Scerevisiae.UCSC.sacCer3"
              ,inputFiles = inputFiles
              ,inputFilesType= "bam"
              ,sampleLabels = c("SL01","SL02","SL03","SL04")
              ,sampleLabelsMerged = c("control","treat")
              ,refSource = "saccharomyces_cerevisiae.SGD.gff"
              ,organismName = "saccharomyces cerevisiae")

## ----eval=TRUE----------------------------------------------------------------
getTSS(myTSSr)
plotCorrelation(myTSSr, samples = "all")
plotTssPCA(myTSSr, TSS.threshold=10)
mergeSamples(myTSSr)
normalizeTSS(myTSSr)
filterTSS(myTSSr, method = "TPM",tpmLow=0.1)

## ----eval=TRUE----------------------------------------------------------------
exportTSStable(myTSSr, data = "raw", merged = "TRUE")
exportTSStoBedgraph(myTSSr, data = "processed", format = "bedGraph")
exportTSStoBedgraph(myTSSr, data = "processed", format = "BigWig")

## ----eval=TRUE----------------------------------------------------------------
clusterTSS(myTSSr, method = "peakclu",peakDistance=100,extensionDistance=30
           ,localThreshold = 0.02,clusterThreshold = 1
	         ,useMultiCore=FALSE, numCores=NULL)
exportClustersTable(myTSSr, data = "filtered")
exportClustersToBed(myTSSr, data = "tagClusters")

## ----eval=TRUE----------------------------------------------------------------
consensusCluster(myTSSr, dis = 50, useMultiCore = FALSE)
exportClustersToBed(myTSSr, data = "consensusClusters")

## ----eval=TRUE----------------------------------------------------------------
shapeCluster(myTSSr,clusters = "consensusClusters", method = "PSS",useMultiCore= FALSE, numCores = NULL)
plotInterQuantile(myTSSr,samples = "all",tagsThreshold = 1)

## ----eval=TRUE----------------------------------------------------------------
exportShapeTable(myTSSr)

## ----eval=TRUE----------------------------------------------------------------
annotateCluster(myTSSr,clusters = "consensusClusters",filterCluster = TRUE,
                filterClusterThreshold = 0.02, annotationType = "genes"
                ,upstream=1000, upstreamOverlap = 500, downstream = 0)

## ----eval=TRUE----------------------------------------------------------------
deGene(myTSSr,comparePairs=list(c("control","treat")), pval = 0.01,useMultiCore=FALSE, numCores=NULL)
exportShiftTable(myTSSr)
plotDE(myTSSr, withGeneName = "TRUE",xlim=c(-2.5, 2.5),ylim=c(0,10))
exportDETable(myTSSr, data = "sig")

## ----eval=TRUE----------------------------------------------------------------
shiftPromoter(myTSSr,comparePairs=list(c("control","treat")), pval = 0.01)
plotTSS(myTSSr,samples=c("control","treat"),tssData = "processed",clusters = "filtered",clusterThreshold = 0.02
        ,genelist=c("YBL017C","YBL067C"),up.dis =500,down.dis = 100)

