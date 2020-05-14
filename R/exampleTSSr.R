
exampleTSSr <- new("TSSr", genomeName = "BSgenome.Scerevisiae.UCSC.sacCer3"
              ,inputFiles = inputFiles
              ,inputFilesType= "bam"
              ,sampleLabels = c("SL01","SL02","SL03","SL04")
              ,sampleLabelsMerged = c("control","treat")
              ,mergeIndex = c(1,1,2,2)
              ,refSource = "saccharomyces_cerevisiae.SGD.gff"
              ,organismName = "saccharomyces cerevisiae")
