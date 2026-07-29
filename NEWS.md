# TSSr 0.99.16 (2026-06-16)

* Mutating workflow functions now return a modified `TSSr` object without
  changing the caller's object. Assign the result to retain each operation:

      Before: mergeSamples(myTSSr)
      After:  myTSSr <- mergeSamples(myTSSr)

* Standardized all 25 exported S4 generics and restricted method dispatch to
  the `object` argument; configuration arguments no longer participate in
  dispatch.

* Removed Bioconductor-build-specific test skips. Core workflow tests, the
  seven workflow examples, and the vignette now execute their documented
  functions on the reduced chromosome I example data.

* Updated the vignette installation instructions to use
  `BiocManager::install("TSSr")`.

* Fixed BAM CIGAR width handling so aligned read intervals are calculated in reference coordinates. Soft clips, hard clips, and insertions no longer inflate the reference span used to determine minus-strand TSS coordinates.
* Restored the BAM terminal correction to the original G-only biological rule while correcting the implementation. With `softclippingAllowed = FALSE`, plus-strand reads trim only leading mismatched `G` bases, and minus-strand reads trim trailing `C` bases as transcript-sense 5' `G` bases. Non-G mismatches are not trimmed, and consecutive mismatched terminal G bases are removed until the first matched G or non-G read base.
* Clarified that `softclippingAllowed = TRUE` uses the aligner's aligned 5' boundary directly and skips uncoded G correction.

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
