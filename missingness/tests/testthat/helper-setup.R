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

# Shared by the onset, collinearity and wear-consistency tests
small_study <- function(beta, seed = 7) {
  cfg <- cfg_test()
  inp <- simulate_inputs(cfg, n_per_cohort = c(Insomnia = 80, Diabetes = 80, Healthy = 160, Prostate = 40),
                         beta = beta, seed = seed)
  list(cfg = cfg, inp = inp,
       b = suppressMessages(prepare_stream_bundle("activity", inp$day_tables, inp$persons, cfg, trailing = TRUE)))
}
primary_spec <- function(cfg, which_stream = "activity") {
  g <- make_onset_grid(cfg)
  g[g$stream == which_stream & g$role == "primary"]
}


# Synthetic device-table rows for a person table: cohort-specific missing shares and device mixes, some scale-only people,
# some people with two devices, some rows without a date.
sim_device_rows <- function(persons, cfg, seed = 11) {
  set.seed(seed)
  mix <- list(Healthy  = c("Charge 5" = .30, "Versa 4" = .20, "MobileTrack" = .15, "Versa 3" = .10, "Inspire 2" = .08, "Sense" = .07, "Alta HR" = .05, "Luxe" = .05),
              Insomnia = c("Charge 5" = .25, "Versa 4" = .20, "MobileTrack" = .25, "Versa 3" = .08, "Inspire 2" = .08, "Sense" = .06, "Alta HR" = .04, "Zip" = .04),
              Diabetes = c("Charge 5" = .25, "Versa 4" = .20, "MobileTrack" = .22, "Versa 3" = .10, "Inspire 2" = .08, "Sense" = .06, "Alta HR" = .05, "Flex" = .04),
              Prostate = c("Charge 5" = .35, "Versa 4" = .25, "Inspire 2" = .10, "Versa 3" = .10, "Sense" = .10, "Charge 4" = .10))
  miss <- c(Healthy = .06, Insomnia = .05, Diabetes = .045, Prostate = 0)
  scale_only <- c(Healthy = .01, Insomnia = .0, Diabetes = .0, Prostate = 0)
  rows <- lapply(seq_len(nrow(persons)), function(i) {
    co <- as.character(persons$cohort[i]); pid <- persons$person_id[i]; u <- runif(1)
    if (u < miss[[co]]) return(NULL)
    if (u < miss[[co]] + scale_only[[co]]) return(data.table(person_id = pid, device_date = as.Date("2022-01-01") + 0:2, device_version = "Aria Air", device_type = "SCALE"))
    v <- sample(names(mix[[co]]), 1, prob = mix[[co]])
    n <- sample(3:12, 1)
    d <- data.table(person_id = pid, device_date = as.Date("2021-01-01") + sort(sample(0:700, n)), device_version = v, device_type = "TRACKER")
    if (runif(1) < .12) d <- rbind(d, data.table(person_id = pid, device_date = as.Date("2022-09-01") + 0:4, device_version = "Versa 4", device_type = "TRACKER"))
    if (runif(1) < .1) d[sample.int(nrow(d), 1L), device_date := as.Date(NA)]
    d
  })
  rbindlist(rows)
}
