# TSSr 0.99.19 (2026-07-30)

* Declared `GenomeInfoDbData` as a suggested dependency and added an
  actionable guard for GFF-based annotation. Existing `refTable` workflows do
  not require it; building a TxDb from GFF with organism metadata now reports
  the exact Bioconductor installation command when the taxonomy data package
  is unavailable.

* Added 13 documented, read-only accessors for TSS matrices, reference data,
  library sizes, and analysis results. Table accessors return independent base
  `data.frame` copies rather than exposing the package's internal `data.table`
  representation. User-facing examples, the README, and the vignette now use
  these accessors instead of direct slot access.

* `DEtables()` explicitly selects either all tested genes or significant genes
  with `result = "all"` or `result = "significant"`.

* Updated the `getTSS()` man-page example to import the bundled TSStable
  fixture. All reviewer-requested user-facing examples now execute their
  documented functions during `R CMD check`.

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

* Added a compact `show()` method for `TSSr` objects. Displaying an object now
  reports its data dimensions and analysis result counts without printing full
  slot contents.

* Fixed `deGene()` so TAG count tables are recomputed from the current assigned
  clusters and accumulated in memory across multiple comparisons. This removes
  stale cached results and the previous disk-backed return channel.

* Updated the reduced example workflow to use an effective TPM filtering
  threshold and to verify that filtering strictly reduces the processed TSS
  matrix.

* Added explicit normalization state to `TSSr` objects. Filtering no longer
  guesses state from column sums, so normalized data remain valid after rows
  are removed and repeated normalization is rejected reliably.

* Promoter-shift `Ds`, `pval`, and `padj` columns are now consistently numeric
  so downstream sorting, plotting, and arithmetic use their intended types.

* Fixed promoter-shift chi-squared tests to use raw counts reconstructed from
  TPM values and sample library sizes. The previous implementation calculated
  the raw-count table but mistakenly tested the TPM table instead. In the
  bundled chromosome I workflow, all 65 testable genes retain the same effect
  sizes, while the genes significant at `pval = 0.01` decrease from 25 to 24;
  `YAR019W-A` is no longer significant. Expected sparse-table approximation
  warnings increase from 35 to 42 and are not suppressed.

* Fixed restoration of graphical parameters after drawing correlation panels.
  Correlation documentation now states how non-positive values are handled on
  logarithmic scatterplot axes.

* Corrected `exportTSStable()` raw export semantics. Raw exports now always
  originate from `TSSrawMatrix`; `merged = TRUE` combines raw sample counts
  according to `mergeIndex` and writes `ALL.samples.TSS.raw.merged.txt`.
  Unmerged raw and processed exports remain distinct as
  `ALL.samples.TSS.raw.txt` and `ALL.samples.TSS.processed.txt`.

* Made TSS table import/export round trips robust to exponent notation.
  Imported genomic positions are validated and stored as integers, exported
  coordinates use integer notation, and TSStable sample values are written
  without scientific notation without changing global R options.

* Allowed one multi-sample TSStable file to be paired with all of its sample
  column labels. Constructor and class validation now compare `mergeIndex`
  with `sampleLabels`, while retaining one-file-per-sample validation for the
  other input formats.

* Moved feature-specific heavy dependencies DESeq2, Gviz, ggfortify, and
  calibrate from Imports to Suggests. Core package loading no longer imports
  them; functions that need them load their namespace on demand and report an
  actionable installation command when a package is unavailable.

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
