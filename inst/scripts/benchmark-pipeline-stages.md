# TSSr pipeline performance benchmarks

This note records the benchmark procedure and the performance decisions made
on the `post-review-improvements` branch in August 2026. The purpose is to
preserve both successful optimizations and negative results so that future
maintainers can reproduce the measurements instead of repeating exploratory
work.

## Benchmark driver

`benchmark-pipeline-stages.R` measures any subset of the public TSSr workflow:

```text
Rscript benchmark-pipeline-stages.R CONFIG OUTDIR \
    [STAGES] [REPETITIONS] [BASELINE_DIR]
```

`STAGES` is either `all` or a comma-separated subset of:

```text
getTSS,merge,normalize,filter,cluster,consensus,shape,annotate,enhancer,de,shift
```

When `BASELINE_DIR` is supplied, the driver restores the snapshot immediately
before each selected stage. This makes before/after comparisons use exactly
the same input object. Each run writes:

- `timings.tsv`, updated after every completed stage;
- one RDS snapshot for each measured stage;
- `warnings.rds` with warnings separated by stage and repetition; and
- `metadata.rds` with package, configuration, input, and session details.

If the original BAM files are not mounted during a snapshot-only benchmark,
their configured paths are retained in the metadata and their MD5 values are
recorded as `NA`. Full-pipeline runs should mount the original inputs so their
checksums can be captured.

## Validation datasets

All retained algorithm changes were tested on two independent four-sample
datasets:

- the yeast BAM dataset used by the original TSSr workflow; and
- an Arabidopsis CAGE BAM dataset.

Candidate stages were run from the same serialized predecessor object as the
baseline. Output tables were compared in full, including values, row order,
column order, and attributes. Floating-point comparisons used a tolerance of
`1e-12` only where the stage already produced floating-point values.

## Retained optimizations

| Stage | Dataset | Baseline | Optimized | Speedup | Output |
|---|---:|---:|---:|---:|---|
| `shapeCluster()` | yeast | 26.863 s | 0.100-0.144 s | 187-269x | equivalent |
| `shapeCluster()` | Arabidopsis | 365.806 s | 0.323-0.376 s | 973-1,133x | equivalent |
| `consensusCluster()` | yeast | 79.832 s | 2.277 s | 35.1x | identical |
| `consensusCluster()` | Arabidopsis | 750.626 s | 8.827 s | 85.0x | identical |
| `deGene()` | yeast | 14.246 s | 2.906 s | 4.9x | equivalent |
| `deGene()` | Arabidopsis | 97.218 s | 3.555 s | 27.3x | equivalent |

The shape and differential-expression changes replace repeated scans of a
complete TSS table with a single point-to-interval mapping followed by grouped
aggregation. The consensus change computes cluster summaries in grouped
passes rather than repeatedly subsetting the sample table. These changes keep
the public API and calculated results unchanged.

## Rejected optimization

A `getTSS()` candidate requested fewer BAM fields and filtered reads by
sequencing quality and chromosome before constructing `GRanges`. On the yeast
four-BAM benchmark it changed runtime from 175.907 seconds to 175.428 seconds,
an improvement of only 0.27%. The TSS matrix and library sizes were identical,
but the gain was below the predeclared 10% retention threshold. The candidate
implementation and its implementation-specific test were therefore removed.

This result indicates that BAM field materialization and temporary `GRanges`
metadata are not important bottlenecks for the current terminal-G correction
workflow. Future `getTSS()` work should be driven by a new profile rather than
retrying this field-pruning approach.

## Verification gate

After the retained changes, the Linux test suite reported:

```text
FAIL 0 | WARN 75 | SKIP 0 | PASS 1412
```

The warnings are the package's already documented statistical and plotting
warnings; none were introduced by these performance changes. A clean Linux
source build, examples, tests, and vignette rebuild completed with `R CMD
check --no-manual` status `OK`. BiocCheck 1.48.0 reported 0 errors, 0 warnings,
and 5 advisory notes. Full package build, check, and BiocCheck results should
be refreshed whenever additional optimizations are added.
