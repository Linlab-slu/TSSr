
###################################################################################################################
##.getTSS_from_TSStable function calls TSS from one TSStable file
##.getTSS_from_tss function takes one TSS table file
##run script with the following example command:
##.getTSS_from_TSStable(TSStable.file)

.getTSS_from_TSStable <- function(TSStable.file){
  if(length(TSStable.file) > 1){
    stop("Only one file should be provided when inputFilesType = \"TSStable\"!")
  }
  if(file.exists(TSStable.file) == FALSE){
    stop("Could not locate input file ", TSStable.file)
  }			
  TSS.all.samples <- read.table(file = TSStable.file, header = F, stringsAsFactors = FALSE)
  if(ncol(TSS.all.samples) != (length(sample.labels) + 3)){
    stop("Number of provided sample labels must match the number of samples in the TSS table!")
  }
  return(TSS.all.samples)
}


