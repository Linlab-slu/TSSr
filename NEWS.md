# TSSr 0.99.20 (2026-08-04)

* Added `clusterTSS(method = "peakcluMax")`, a greedy maximum-signal
  alternative to the published `peakclu` method. It repeatedly retains the
  strongest remaining TSS, resolves ties toward the lower genomic position,
  and suppresses candidates in the retained peak's closed `peakDistance`
  window. Peak selection is implemented in O(n log n) time and O(n) memory;
  the default remains `method = "peakclu"`.

* Made the 10th- and 90th-percentile TSS cluster boundaries deterministic
  across platforms. Cumulative signal that equals a percentile threshold
  within a scale-aware floating-point tolerance is now treated as having
  reached that threshold, consistently with the documented definition that
  the interval contains at least 80% of the cluster signal. The same rule is
  used for tag clusters and consensus clusters in single- and multicore code.

* Added public-workflow regression tests for exact-threshold behavior and
  invariance under positive signal scaling. Cluster coordinates and other
  discrete fields remain exact in cross-platform comparisons; only aggregate
  floating-point tag sums use a relative tolerance. On the bundled chromosome
  I data, the deterministic rule changes one control cluster's lower boundary
  from 123531 to 123502 and its interquantile width from 16 to 45. Its promoter
  shape score changes accordingly; assigned clusters, enhancers, differential
  expression, TAG tables, and promoter-shift results are unchanged.

* Reduced repeated work in tests and examples without skipping any code paths.
  Workflow examples now run on the bundled 100-row, four-sample, two-strand
  TSStable input rather than recomputing against the full example object. On
  the Linux Bioconductor devel container used for validation, no example takes
  more than five seconds and `R CMD check --no-vignettes --no-manual` decreased
  from 566 to 421 seconds while all 1,000 test assertions still run.

# TSSr 0.99.19 (2026-07-30)

* `getTSS()` now streams BAM and paired-end BAM inputs in bounded chunks.
  The new `bamYieldSize` argument controls the number of alignments read per
  chunk and defaults to one million. Chunk processing stops only at end of
  file, so exact chunk multiples and fully filtered inputs cannot create an
  extra empty record. Chunked imports retain the current reference-CIGAR,
  strand, terminal-G, quality, and proper-first-mate behavior.

* `exportTSStable()` now accepts `outputPrefix` and `outputSuffix` while
  retaining the established `TSS.raw.txt`, `TSS.raw.merged.txt`, and
  `TSS.processed.txt` suffixes by default. Existing output files are no longer
  silently replaced: callers must explicitly set `overwrite = TRUE` to
  replace one.

* Added `combineTSSTables()` and `splitTSSTable()` for assembling any number
  of independent raw-count TSS tables and extracting explicit sample columns.
  The functions align the union of `chr`, `pos`, and `strand` coordinates,
  fill missing sample values with zero, remove all-zero rows when requested,
  validate raw counts and chromosome compatibility, and preserve versioned
  data type, genome, and chromosome metadata. `exportTSStable()` now writes
  the same metadata so processed tables cannot be mistaken for raw counts.

* `sampleLabelsMerged` and `mergeIndex` are now optional when no samples need
  to be combined. Omitting both creates a one-to-one mapping from each sample
  to its own analysis group, so `mergeSamples()` can be skipped. Supplying
  only one of the two grouping arguments now gives an immediate error.

* GFF-based annotation no longer requires `GenomeInfoDbData` or a taxonomic
  lookup. `organismName` remains available as optional descriptive metadata
  for backward compatibility, but cluster-to-gene assignment is determined
  from the coordinates and strands in `refSource` or `refTable`. Earlier
  versions passed this value only to metadata on a temporary TxDb object; it
  did not participate in the annotation calculation.

* Added a reproducible GFF3 fixture containing the first 12 complete genes and
  associated features from the original SGD R64-2-1 chromosome-I annotation.
  Public `annotateCluster()` tests compare its results with an equivalent
  `refTable` workflow on both strands.

* Added public-API import coverage for BED, per-sample TSS tables, BigWig, and
  paired-end BAM. The tests verify format-specific coordinate, strand,
  aggregation, missing-site, signed-score, and proper-first-mate behavior. A
  reproducible 346-byte paired BAM fixture is generated from real sacCer3
  sequence and requires no BAM index.

* Added 13 documented, read-only accessors for TSS matrices, reference data,
  library sizes, and analysis results. Table accessors return independent base
  `data.frame` copies rather than exposing the package's internal `data.table`
  representation. User-facing examples, the README, and the vignette now use
  these accessors instead of direct slot access.

* `DEtables()` explicitly selects either all tested genes or significant genes
  with `result = "all"` or `result = "significant"`.

* Updated the `getTSS()` man-page example to import a reproducible 100-row,
  four-sample TSStable fixture containing 50 TSSs from each strand. All
  reviewer-requested user-facing examples now execute their documented
  functions during `R CMD check`.

* Added a reproducible, unindexed micro-BAM fixture built from real sacCer3
  chromosome-I sequence. Its public `getTSS()` regression test covers both
  strands, insertion/deletion/soft-clipping CIGAR operations, terminal
  matched and mismatched G correction, quality filtering, and count
  aggregation. The expected minus-strand coordinates distinguish reference
  width from the former query-width implementation.

* Removed inert legacy QA assertions that did not exercise TSSr. Single-sample
  consensus, DESeq2 output columns, and plot/export prerequisites now have
  behavior-level assertions, so missing results cannot produce an empty PASS.

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
  warnings increase from 35 to 42 and are not suppressed. The p-value direction
  is scale-dependent: reconstructed counts reduce effective sample size when
  `librarySize / 1e6` is below one and increase it when the factor is above
  one; effect-size calculations are unchanged.

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

* Optimized `clusterTSS(method = "peakclu")` peak detection by replacing
  repeated per-row window scans with indexed range-maximum lookups. Permanent
  regression tests cover open-window boundaries, tied signals, public return
  semantics, and exact agreement with the previously generated tag clusters.

* Removed an unused genome load from `clusterTSS()`, reducing fixed overhead
  without changing clustering output.

* Historical note: the untagged development commit `09fa236` temporarily used
  version 0.99.17 for the peakclu optimization. The published `v0.99.17` tag
  identifies a different Bioconductor-review checkpoint; both histories are
  reconciled in 0.99.19.

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
