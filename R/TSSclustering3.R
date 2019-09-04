# This is a script to cluster cage sequencing data
# Input is a TSS table for one condition (4 columns: chromosome, position, strand, and tag counts)
# This script consists of three functions:
# main --> clusterTSS
# filter --> filterWithPoisson
# cluster --> clusterByPeak

# After loading all functions, run script with the following example command:
# tss.clusters <- .clusterTSS2(data, Scerevisiae, useMultiCore = TRUE)
# where "data" is tss table, "Scerevisiae" is the loaded genome, and "useMultiCore" allows for faster parallelization of clustering

###############################################################################

.clusterTSS <- function(tss.dt, genome, filter = "TPM",pVal=0.01, tpmLow = 0.1, normalizeData=TRUE, 
                         peakDistance=100, localThreshold = 0.02, extensionDistance=30, nonOverLapping=TRUE,
                         useMultiCore=FALSE, numCores=NULL)  
{
  # force columns 1 to 3 to characters
  
  # convert dataframe to datatable
  setDT(tss.dt)
  setnames(tss.dt, colnames(tss.dt)[[4]], "tags")
  
  # calculate size of genome
  genomeSize <- 0
  for (chrom in 1:length(genome)) {
    genomeSize <- genomeSize + length(genome[[chrom]])
  }
  
  # find sequence coverage depth
  coverageDepth <- tss.dt[,sum(tags)]
  if(filter = "poisson"){
    tss.dt <- .filterWithPoisson(tss.dt, coverageDepth, genomeSize, pVal)
    if (normalizeData) {
      sizePerMillion <- coverageDepth / 1e6
      tss.dt[, tags := round(tags / sizePerMillion, 6)]
    }
  }else if(filter = "TPM"){
    sizePerMillion <- coverageDepth / 1e6
    tss.dt[, tags := round(tags / sizePerMillion, 6)]
    tss.dt <- tss.dt[tags >= tpmLow,]
  }
  
  # filter tss data
  
  # normalize tags (tags per million)
  # if false assume data is already normalized as tpm
  
  # add column to subset data
  tss.dt[, "subset" := paste0(chr,"_",strand)]
  setkey(tss.dt,subset)
  
  # pass sub datatables to peak-caller and clustering functions
  if (useMultiCore) {
    library(parallel)
    if (is.null(numCores)) {
      numCores <- detectCores()
    }
    print(paste("process is running on", numCores, "cores..."))
    clusters <- mclapply(as.list(tss.dt[,unique(subset)]), function(x) {
      data <- tss.dt[.(x)]
      setkey(data, NULL)
      setorder(data, pos)
      cluster.data <- .clusterByPeak(data, peakDistance, localThreshold, extensionDistance, nonOverLapping)
    }, mc.cores = numCores)
  }else{
    clusters <- lapply(as.list(tss.dt[,unique(subset)]), function(x) {
      data <- tss.dt[.(x)]
      setkey(data, NULL)
      setorder(data, pos)
      cluster.data <- .clusterByPeak(data, peakDistance, localThreshold, extensionDistance, nonOverLapping)
    })
  }
  
  tss.clusters <- rbindlist(clusters, use.names=TRUE, fill=TRUE)
  na.omit(tss.clusters)
  # update cluster IDs
  tss.clusters[, cluster := .I]
  # covert list of datatables to dataframe object and return
  return(setDF(tss.clusters))
} 

###############################################################################

.filterWithPoisson <- function(data, coverageDepth, genomeSize, pVal) {
  # calculate lambda value (average)
  lambda <- coverageDepth / (genomeSize * 2)
  # get cutoff value
  cutoff <- qpois(pVal, lambda, lower.tail=FALSE, log.p=FALSE)
  # filter tss table
  data <- data[tags >= cutoff,]
  return(data)
}

###############################################################################

.clusterByPeak <- function(tss.dt, peakDistance, localThreshold, extensionDistance, nonOverLapping) {
  # create copy for reference later
  copied.dt <- copy(tss.dt)
  setkey(tss.dt, pos)
  
  # get peakID
  # TODO could potentially by optimized more
  peakID <- vapply(seq_len(tss.dt[,.N]), function(x) {
    id <- 0
    temp <- tss.dt[x,pos]
    ##ZL
    if(tss.dt[x,pos] == tss.dt[pos>temp-peakDistance & pos<temp+peakDistance,][which(tags == max(tags)),pos][1]){
      id <- x
    }else{id <- 0}
    return(id)
  }, numeric(1))
  
  # manipulate data.table to collapse clustered rows
  tss.dt[, peak := peakID]
  tss.dt[, ID := .I]
  #######################################################################################################################
  ##local filtering
  #######################################################################################################################
  if(unique(tss.dt$strand)== "+"){
    localF <- sapply(peakID[peakID >0],function(i){
      temp <- tss.dt[pos %between% c(tss.dt$pos[i],tss.dt$pos[i]+100),]
      temp$ID[which(temp$tag < tss.dt$tags[i] * localThreshold)]
    })}
  else{
    localF <- sapply(peakID[peakID >0],function(i){
      temp <- tss.dt[pos %between% c(tss.dt$pos[i]-100,tss.dt$pos[i]),]
      temp$ID[which(temp$tag < tss.dt$tags[i] * localThreshold)]
    })}
  if(length(unlist(localF)) >0){tss.dt <- tss.dt[-unlist(localF),]}
  #######################################################################################################################
  #######################################################################################################################
  tss.dt[, forward := ifelse(data.table::shift(pos,1,type="lead") < pos + extensionDistance, 1, 0)] #  
  tss.dt[, reverse := ifelse(data.table::shift(pos,1,type="lag") > pos - extensionDistance, 1, 0)]
  tss.dt <- tss.dt[,list(peak=max(peak),start=min(pos),end=max(pos),tags=sum(tags)),by=.(rleid(peak, forward, reverse))]##ZL?
  
  # get start and end boundaries for clusters
  # TODO revisit this code for better optimization
  clusters <- lapply(as.list(tss.dt[peak>0,rleid]), function(x) {
    start <- tss.dt[x,start]
    end <- tss.dt[x,end]
    
    if (x-1>0 && tss.dt[x-1,!peak>0] && tss.dt[x-1,end] > start - extensionDistance) {
      start <- tss.dt[x-1,start]
      if (x-2>0 && tss.dt[x-2,!peak>0] && tss.dt[x-2,end] > start - extensionDistance) {
        start <- tss.dt[x-2,start]
      }
    }
    if (x+1<tss.dt[,.N] && tss.dt[x+1,!peak>0] && tss.dt[x+1,start] < end + extensionDistance) {
      end <- tss.dt[x+1,end]
      if (x+2<tss.dt[,.N] && tss.dt[x+2,!peak>0] && tss.dt[x+2,start] < end + extensionDistance) {
        end <- tss.dt[x+2,end]
      }
    }
    list(start, end)
  })
  
  clusters <- rbindlist(clusters)
  
  # deal with overlapping clusters here
  # TODO this section needs some more optimization/work
  if (nonOverLapping) {
    rowVec <- which(clusters$V2 >= data.table::shift(clusters$V1,1,type="lead"))
    if (length(rowVec)>0) {
      #######################################################################################################################
      #######################################################################################################################
      for(i in 1:length(rowVec)){clusters$V1[rowVec[i]+1] = clusters$V1[rowVec[i]]}
      clusters <- clusters[-rowVec,]
      #######################################################################################################################
      #######################################################################################################################
#      newbounds <- mapply(function(x,y) {
#        temp <- copied.dt[between(pos, x, y, incbounds=FALSE),]
#        ##ZL 
#        
#        nb <- 0
#       if (length(temp[tags==min(tags),pos]) > 1) {
#          nb <- temp[tags==min(tags),pos][[1]]  # where should we split??
#        }else { nb <- temp[tags==min(tags),pos] }
#        nb
#      }, x=clusters$V1[rowVec+1],y=clusters$V2[rowVec])##mapply end
#      #######################################################################################################################
#      #######################################################################################################################
#      clusters$V2[rowVec] <- newbounds
#      clusters$V1[rowVec+1] <- newbounds + c(1)
#      
#      tooSmall <- which(clusters$V1[rowVec]+1 >= clusters$V2[rowVec] | clusters$V1[rowVec+1]+1 >= clusters$V2[rowVec+1])
#      if (length(tooSmall)>0) {
#        clusters$V2[rowVec][tooSmall] <- clusters$V2[rowVec+1][tooSmall]
#        clusters$V1[rowVec+1][tooSmall] <- clusters$V1[rowVec][tooSmall]
#      }
    }##
  }
  
#  clusters <- unique(clusters)
  
  # get full clustering data
  # core promoter boundaries are calculated here (i.e. cumsum distribution)
  tss_clusters <- lapply(as.list(seq_len(clusters[,.N])), function(i) {
    start <- clusters[i,V1]
    end <- clusters[i,V2]
    #copied.dt[, ID := .I]##NEW April18
    cluster.data <- copied.dt[pos %between% c(start,end), ]
    tpm <- cluster.data[,sum(tags)]
    q1 <- cluster.data[which(cumsum(tags) > 0.1*tpm),min(pos)]
    q9 <- cluster.data[order(-pos)][which(cumsum(tags) > 0.1*tpm),max(pos)]
    list(i
         ,cluster.data[,chr[[1]]]
         ,start
         ,end
         ,cluster.data[,strand[[1]]]
         #,cluster.data[which(ID %in% peakID), pos]  # NEW - use id column to find the intersection for the peakID vector (should hopefully only be 1!)
         ,cluster.data[which.max(tags),pos]
         ,tpm
         ,cluster.data[,max(tags)]
         ,q1
         ,q9
         ,q9 - q1 + 1)
  })
  
  # set names
  tss_clusters <- rbindlist(tss_clusters)
  setnames(tss_clusters, c( "cluster"
                             , "chr", "start", "end", "strand"
                             , "dominant_tss", "tpm", "tpm.dominant_tss"
                             , "q_0.1", "q_0.9", "interquantile_width" ))
  return(tss_clusters)
}