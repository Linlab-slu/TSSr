## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----github, eval=FALSE-------------------------------------------------------
#  devtools::install_github("Linlab-slu/TSSr")

## ----library, results='hide', message=FALSE-----------------------------------
library(TSSr)

## ----citation, eval=TRUE------------------------------------------------------
citation("TSSr")

## ----exampleTSSr, results="hide", tidy=FALSE----------------------------------
# Load the example data
data("exampleTSSr")
myTSSr <- exampleTSSr

## ----newTSSr, eval=FALSE------------------------------------------------------
#  # Provide bam files
#  inputFiles <- c("S01.sorted.bam", "S02.sorted.bam", "S03.sorted.bam", "S04.sorted.bam")
#  myTSSr <- new("TSSr", genomeName = "BSgenome.Scerevisiae.UCSC.sacCer3"
#            ,inputFiles = inputFiles
#            ,inputFilesType= "bam"
#            ,sampleLabels = c("SL01","SL02","SL03","SL04")
#            ,sampleLabelsMerged = c("control","treat")
#            ,mergeIndex = c(1,1,2,2)
#            ,refSource = "saccharomyces_cerevisiae.SGD.gff"
#            ,organismName = "saccharomyces cerevisiae")
#  myTSSr

## ----getTSS, eval=FALSE-------------------------------------------------------
#  # Get TSS
#  getTSS(myTSSr)

## ----TSSprocessing, tidy=FALSE------------------------------------------------
# Merge replicates
mergeSamples(myTSSr)
# Normalization
normalizeTSS(myTSSr)
# TSS filtering
filterTSS(myTSSr, method = "TPM", tpmLow = 0.1)

## ----TSSclustering, eval=FALSE------------------------------------------------
#  # TSS clustering
#  clusterTSS(myTSSr, method = "peakclu",peakDistance=100,extensionDistance=30
#           ,localThreshold = 0.02,clusterThreshold = 1
#           ,useMultiCore=FALSE, numCores=NULL)
#  
#  # Aggregating consensus clusters
#  consensusCluster(myTSSr, dis = 50, useMultiCore = FALSE)

## ----shapeCluster, eval=FALSE-------------------------------------------------
#  # Calculating core promoter shape score
#  shapeCluster(myTSSr,clusters = "consensusClusters", method = "PSS",
#               useMultiCore= FALSE, numCores = NULL)

## ----annotateCluster, tidy=FALSE----------------------------------------------
# Assign clusters to the annotated features
annotateCluster(myTSSr,clusters = "consensusClusters",filterCluster = TRUE,
              filterClusterThreshold = 0.02, annotationType = "genes"
              ,upstream=1000, upstreamOverlap = 500, downstream = 0)

## ----deGEne, tidy=FALSE-------------------------------------------------------
# Assign clusters to the annotated features
deGene(myTSSr,comparePairs=list(c("control","treat")), 
       pval = 0.01,useMultiCore=FALSE, numCores=NULL)

## ----shiftPromoter, tidy=FALSE------------------------------------------------
# Calcuate core promoter shifts
shiftPromoter(myTSSr,comparePairs=list(c("control","treat")), pval = 0.01)

