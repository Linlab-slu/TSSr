# Test TSSr S4 class creation and validation

test_that("TSSr class can be instantiated with valid parameters", {
  # Create a minimal TSSr object

  obj <- new("TSSr",
             genomeName = "BSgenome.Scerevisiae.UCSC.sacCer3",
             inputFilesType = "bam",
             sampleLabels = c("sample1", "sample2"))


  expect_s4_class(obj, "TSSr")
  expect_equal(obj@genomeName, "BSgenome.Scerevisiae.UCSC.sacCer3")
  expect_equal(obj@inputFilesType, "bam")
  expect_equal(obj@sampleLabels, c("sample1", "sample2"))
})

test_that("TSSr class validates inputFilesType", {
  # Test that invalid inputFilesType is rejected
  expect_error(
    new("TSSr",
        genomeName = "BSgenome.Scerevisiae.UCSC.sacCer3",
        inputFilesType = "invalid_type",
        sampleLabels = c("sample1")),
    regexp = "inputFilesType.*must be one of"
  )
})

test_that("TSSr class accepts all valid inputFilesType values", {
  valid_types <- c("bam", "bamPairedEnd", "bed", "tss", "TSStable", "BigWig")

  for (type in valid_types) {
    obj <- new("TSSr",
               genomeName = "test_genome",
               inputFilesType = type,
               sampleLabels = c("sample1"))
    expect_s4_class(obj, "TSSr")
    expect_equal(obj@inputFilesType, type)
  }
})

test_that("TSSr class has correct default slot values", {
  obj <- new("TSSr",
             genomeName = "test",
             inputFilesType = "bam",
             sampleLabels = "s1")


  # Check default empty values

  expect_equal(length(obj@TSSrawMatrix), 0)
  expect_equal(length(obj@tagClusters), 0)
  expect_equal(length(obj@consensusClusters), 0)
  expect_equal(length(obj@DEtables), 0)
  expect_type(obj@tagClusters, "list")
  expect_type(obj@consensusClusters, "list")
})
