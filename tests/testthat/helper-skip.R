# testthat::skip_on_bioc() only checks IS_BIOC_BUILD_MACHINE, which is
# not set on the Bioconductor Single Package Builder (verified empirically
# on builds 0.99.10–0.99.14). This helper adds redundant signals that DO
# fire on SPB workers, while remaining inert on devtools::test, CI, and
# user installs so coverage there is unchanged.
.skip_on_bioc_spb <- function() {
  on_bioc <- identical(Sys.info()[["user"]], "biocbuild") ||
    grepl("bbs-", R.home(), fixed = TRUE) ||
    nzchar(Sys.getenv("BBS_HOME")) ||
    identical(Sys.getenv("IS_BIOC_BUILD_MACHINE"), "true")
  if (on_bioc) {
    testthat::skip("Skipping slow test on Bioconductor build machine")
  }
}
