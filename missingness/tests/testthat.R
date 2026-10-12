# Run from the missingness/ directory:   Rscript tests/testthat.R
library(testthat)
source("load_missingness.R"); load_missingness(".")
test_dir("tests/testthat", reporter = "summary", stop_on_failure = TRUE)
