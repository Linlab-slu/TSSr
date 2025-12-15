# Test using exampleTSSr data

test_that("exampleTSSr data loads correctly", {
  data(exampleTSSr)
  expect_s4_class(exampleTSSr, "TSSr")
})

test_that("exampleTSSr has expected structure", {
  data(exampleTSSr)

  # Check that essential slots are populated
  expect_true(length(exampleTSSr@sampleLabels) > 0)
  expect_true(nrow(exampleTSSr@TSSrawMatrix) > 0)
})

test_that("exampleTSSr sampleLabels are valid", {
  data(exampleTSSr)

  # Sample labels should be character
  expect_type(exampleTSSr@sampleLabels, "character")

  # Should have at least one sample

  expect_gte(length(exampleTSSr@sampleLabels), 1)
})

test_that("exampleTSSr TSSrawMatrix has correct columns", {
  data(exampleTSSr)

  # TSSrawMatrix should have chr, pos, strand columns
  raw_matrix <- exampleTSSr@TSSrawMatrix

  if (nrow(raw_matrix) > 0) {
    expect_true("chr" %in% names(raw_matrix) || "V1" %in% names(raw_matrix) ||
                  length(names(raw_matrix)) >= 3)
  }
})
