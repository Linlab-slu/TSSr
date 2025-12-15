# TSSr 0.99.6

* Preparing for Bioconductor submission
* Updated LICENSE to standard MIT
* Updated CITATION format to bibentry()
* Moved BSgenome.Scerevisiae.UCSC.sacCer3 to Suggests

# TSSr 0.1.0

* Initial release (2021-03-08)
* Code has been tested and is stable
* Core functionality complete:
  - TSS identification from BAM, bed, BigWig, ctss, and TSS table files
  - TSS clustering using peakclu algorithm
  - Consensus cluster generation across samples
  - Core promoter shape analysis (PSS and SI metrics)
  - Cluster annotation to downstream genes
  - Differential expression analysis via DESeq2
  - Promoter shift detection between conditions
  - Export to bedGraph, BigWig, BED, and table formats
