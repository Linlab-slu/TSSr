##########################################################################################################
##TSSr pacakge main page
##all packages need to be loaded before run the following codes
library(Rsamtools)
library(stringr)
library(data.table)
library(rtracklayer)
library(GenomicFeatures)
library(Gviz)
library(GenomicRanges)
library(Deseq2)
library(BSgenome.Scerevisiae.UCSC.sacCer3)
##########################################################################################################
samplex <- c("ScerBY4741.1","ScerBY4741.2")
sampley <- c("ScerArrest.1","ScerArrest.2")
sampleMergedx <- "ScerBY4741"
sampleMergedy <- "ScerArrest"
sample.labels <- c("ScerBY4741.1","ScerBY4741.2","ScerArrest.1","ScerArrest.2")
organism = "Saccharomyces cerevisiae"
reference <- "~/saccharomyces_cerevisiae.SGD.gff"
genome <- BSgenome.Scerevisiae.UCSC.sacCer3
##If input files are bam files
bam.files <- c("ScerBY4741.1.sorted.bam","ScerBY4741.2.sorted.bam","ScerArrest.1.sorted.bam","ScerArrest.2.sorted.bam")

TSS.all.samples <- .getTSS_from_bam(bam.files, genome, sequencingQualityThreshold = 10, mappingQualityThreshold = 20,removeNewG = TRUE,correctG = TRUE)

.plotCorrelation(TSS.all.samples)


library.sizes <- .library.sizes(TSS.all.samples)

clusters <- clusterTSS2(tss,genome, pVal=0.01, normalizeData=TRUE, 
                          peakDistance=100, extensionDistance=30, nonOverLapping=TRUE,
                          useMultiCore=FALSE, numCores=NULL)
pss <- .pss(tss.tpm,clusters)
##write.table(clusters, file = paste(samplesLablesMerged[i],"cluster.txt", sep = "."), sep = "\t", row.names = F, quote = F)
##provide gff reference to be assigned
clusters.asn <- .assign2gene(clusters,reference,organims,upstream=500, downstream = 0)
##write.table(clusters.asn, file = paste(samplesLablesMerged[i],"cluster.txt", sep = "."), sep = "\t", row.names = F, quote = F)
Ds <- .Ds(clustersx.asn,clustersy.asn,librarySizex, librarySizey)

De <- .deseq2(clustersx.asn,clustersy.asn, tss.all, samplex,sampley, sampleMergedx,sampleMergedy)
DeSig <- subset(De, padj <= 0.05)
##write.table(as.data.frame(DeSig),file = "Deseq2_sig", sep = "\t", row.names = T)

##provide gene.list which you want to plot
.plotTSS(tss.tpm,clusters.cl, gene.list, up.dis = 500, down.dis=100)

