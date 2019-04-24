
###################################################################################################################
##.getTSS_from_tss function calls TSS from bam file
##.getTSS_from_tss function takes tss files
##run script with the following example command:
##.getTSS_from_tss(tss.files)

.getTSS_from_tss <- function(tss.files){
  first <- TRUE
  for(i in 1:length(tss.files)) {
    message("\nReading in file: ", tss.files[i], "...")
    TSS <- read.table(file = tss.files[i], header = F, sep = "\t", colClasses = c("character", "integer", "character", "integer"), col.names = c("chr", "pos", "strand", sample.labels[i]))
    TSS <- data.table(TSS)
    setkeyv(TSS, cols = c("chr", "pos", "strand"))
    if(first == TRUE) {
      TSS.all.samples <- TSS
    }else{
      TSS.all.samples <- merge(TSS.all.samples, TSS, all = TRUE)
    }			
    first <- FALSE
  }
  TSS.all.samples <- data.frame(TSS.all.samples)
  return(TSS.all.samples)
}
