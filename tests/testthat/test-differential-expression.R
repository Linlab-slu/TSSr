# Test differential expression functionality

test_that("DEtables slot has valid structure", {
    data(exampleTSSr)

    expect_type(exampleTSSr@DEtables, "list")
    expect_true(length(exampleTSSr@DEtables) > 0)

    first_de <- exampleTSSr@DEtables[[1]]
    expect_type(first_de, "list")
    expect_true(length(first_de) > 0)
    expect_true("DEtable" %in% names(first_de))
})

test_that("DEtable contains expected DESeq2 columns", {
    data(exampleTSSr)

    first_de <- exampleTSSr@DEtables[[1]]
    de_table <- first_de$DEtable

    expect_true(is.data.frame(de_table))
    expect_true(nrow(de_table) > 0)

    col_names <- names(de_table)
    deseq_cols <- c("baseMean", "log2FoldChange", "pvalue", "padj")
    has_deseq_cols <- any(deseq_cols %in% col_names)
    expect_true(has_deseq_cols)
})
