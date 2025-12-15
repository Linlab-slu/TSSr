################################################################################################
.filterWithPoisson <- function(data, coverageDepth, genomeSize, pVal) {
    # calculate lambda value (average)
    lambda <- coverageDepth / (genomeSize * 2)
    # get cutoff value
    cutoff <- qpois(pVal, lambda, lower.tail = FALSE, log.p = FALSE)
    # fiter tss table
    data[data < cutoff, ] <- 0
    return(data)
}
