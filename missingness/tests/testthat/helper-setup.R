# Shared fixtures. Functions are loaded by tests/testthat.R (or here, if a test file
# is run on its own from the missingness/ directory).
if (!exists("make_config")) {
  root <- if (file.exists("load_missingness.R")) "." else if (file.exists("../../load_missingness.R")) "../.." else stop("run from missingness/")
  source(file.path(root, "load_missingness.R")); load_missingness(root)
}
suppressMessages(library(data.table))

cfg_test <- function(...) make_config(list(boot = list(B = 30L), ...))

# Brute-force gap detection on an explicit calendar grid (the slow, obviously-correct way)
naive_gaps <- function(dates) {
  dates <- sort(unique(dates))
  full <- seq(min(dates), max(dates), by = "day")
  miss <- !(full %in% dates)
  r <- rle(miss); ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1L
  idx <- which(r$values)
  data.table(gap_start = full[starts[idx]], gap_length = r$lengths[idx])
}

# Random day table: n persons, each with a random window and random observed days
random_days <- function(n = 25, seed = 1) {
  set.seed(seed)
  rbindlist(lapply(seq_len(n), function(i) {
    W <- sample(30:120, 1); p <- runif(1, 0.6, 0.97)
    d0 <- as.Date("2023-01-01") + sample(0:200, 1)
    obs <- c(1L, which(runif(W - 1L) < p) + 1L)
    data.table(person_id = 100L + i, date = d0 + obs - 1L)
  }))
}
