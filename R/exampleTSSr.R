
exampleFiles <- list.files(system.file("extdata", package = "TSSr")
                           , "bam$", full.names = TRUE)
exampleGFF <- list.files(system.file("extdata", package = "TSSr")
                         , "gff$", full.names = TRUE)
exampleTSSr <- new("TSSr", genomeName = "BSgenome.Scerevisiae.UCSC.sacCer3"
              ,inputFiles = exampleFiles
              ,inputFilesType= "bam"
              ,sampleLabels = c("SL01","SL02","SL03","SL04")
              ,sampleLabelsMerged = c("control","treat")
              ,mergeIndex = c(1,1,2,2)
              ,refSource = "saccharomyces_cerevisiae.SGD.gff"
              ,organismName = "saccharomyces cerevisiae")

