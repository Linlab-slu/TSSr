


.library.sizes <- function(TSS.all.samples)
  {as.integer(colSums(TSS.all.samples[,c(4:ncol(TSS.all.samples)), drop = F], na.rm = T))}
