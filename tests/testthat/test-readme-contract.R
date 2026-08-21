test_that("README uses the validated public constructor", {
    readme_file <- system.file(
        "README.md",
        package = "TSSr",
        mustWork = TRUE
    )
    readme <- readLines(readme_file, warn = FALSE)

    expect_false(any(grepl(
        "new\\s*\\(\\s*[\"']TSSr[\"']",
        readme
    )))
    expect_true(any(grepl(
        "myTSSr\\s*<-\\s*TSSr\\s*\\(",
        readme
    )))
})
