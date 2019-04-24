#####################################################################################################
#####################################################################################################
##.plotCorrelations plots pair-wise pearson correlations 
##.plotCorrelations function takes one input file TSS.all.samples
##run script with the following example command:
##.plotCorrelations(TSS.all.samples)



##########################################################################################################
.plotCorrelations <- function(TSS.all.samples){
  z <- TSS.all.samples[,-c(1:3)]
  # Customize lower panel
  panel.cor <- function(x, y){
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1), xlog = FALSE, ylog = FALSE)
    r <- round(cor(x, y), digits=2)
    txt <- paste0( r)
    cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * r)
  }
  # Customize upper panel
  upper.panel<-function(x, y){
    points(x,y, pch = ".",col = "#00AFBB")
  }
  # Create the plots
  suppressWarnings(pairs(z, lower.panel = upper.panel,upper.panel = panel.cor, log = "xy"))
}




