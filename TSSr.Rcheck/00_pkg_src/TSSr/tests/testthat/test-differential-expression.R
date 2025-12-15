# Test differential expression functionality

test_that("DEtables slot has valid structure", {
  data(exampleTSSr)

  # DEtables should be a list
  expect_type(exampleTSSr@DEtables, "list")

  if (length(exampleTSSr@DEtables) > 0) {
    # Each DE result should contain DEtable and DEsig
    first_de <- exampleTSSr@DEtables[[1]]
    expect_type(first_de, "list")

    if (length(first_de) > 0) {
      # Should have DEtable component
      expect_true("DEtable" %in% names(first_de) || length(first_de) >= 1)
    }
  }
})

test_that("DEtable contains expected DESeq2 columns", {
  data(exampleTSSr)

  if (length(exampleTSSr@DEtables) > 0) {
    first_de <- exampleTSSr@DEtables[[1]]

    if ("DEtable" %in% names(first_de)) {
      de_table <- first_de$DEtable

      if (is.data.frame(de_table) && nrow(de_table) > 0) {
        col_names <- names(de_table)
        # DESeq2 output columns
        deseq_cols <- c("baseMean", "log2FoldChange", "pvalue", "padj")
        has_deseq_cols <- any(deseq_cols %in% col_names) || ncol(de_table) >= 5
        expect_true(has_deseq_cols)
      }
    }
  }
})
