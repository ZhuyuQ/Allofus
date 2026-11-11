# Stream model: risk set, outcome, exposure without leakage, recording history, tests, the joint sandwich (R/stream_model.R, R/plots_stream_model.R).
for (f in c("R/stream_states.R", "R/plots_stream_states.R", "R/stream_model.R", "R/plots_stream_model.R")) source(file.path(if (file.exists(f)) "." else "../..", f))

d0 <- as.Date("2024-03-04")

# One person, 22 days, built by hand (see the expectations below): days 1-14 fully recorded (day 5 with 18 hours of heart rate), day 15 activity only,
# days 16-20 fully recorded, day 21 excluded (device without heart rate), day 22 fully recorded.
make_days <- function(steps_fun = function(j) 1000 * j) {
  j <- 1:22
  D <- data.table(person_id = 1L, date = d0 + j - 1L, reason = "analysed", act_after = TRUE, hr_present = TRUE, sleep_main = TRUE, slp_row = TRUE, steps = steps_fun(j),
                  dev_version = "Versa 3", hr_min = 1300)
  D[j == 5L, hr_min := 1080]; D[j == 15L, `:=`(hr_present = FALSE, sleep_main = FALSE, hr_min = NA_real_)]; D[j == 21L, reason := "device_without_hr"]
  D[]
}

test_that("the outcome categories: every other pattern is category 4", {
  expect_equal(stm_outcome(c("111", "000", "100", "110", "101", "001", "011", "010", NA)), c(0L, 1L, 2L, 3L, 4L, 4L, 4L, 4L, NA))
})

test_that("the risk set: fully recorded d with 20 hours, seven earlier baseline days, outcome from d + 1", {
  sc <- stm_config(); X <- stm_days(make_days(), sc)
  persons <- data.table(person_id = 1L, cohort = "Healthy", sex = "Female", age_first = 50, cci_score = 1, cci_group = "1-2")
  rs <- stm_risk_set(X, persons, cfg_test(), sc); R <- rs$R
  expect_equal(as.integer(R$date - d0) + 1L, c(9:14, 16:19))                          # day 5 (18 h) is no baseline day; days 9+ have seven earlier ones; day 15 (activity only) is no risk-set day; day 20 precedes an excluded day; day 22 is last
  expect_equal(R$outcome, c(0L, 0L, 0L, 0L, 0L, 2L, 0L, 0L, 0L, 0L))                 # day 14 is followed by the activity-only day 15
  expect_equal(rs$flow$n_rows, c(20L, 19L, 12L, 11L, 10L)); expect_equal(rs$flow$n_dropped[-1L], c(1L, 7L, 1L, 1L))
})

test_that("the usual level uses only EARLIER qualifying days: nothing from day d on enters it", {
  sc <- stm_config(); X <- stm_days(make_days(), sc)
  row14 <- X[date == d0 + 13L]                                                       # d = day 14: baseline = qualifying days 1-13 without day 5
  expect_equal(row14$m_before, mean(log1p(1000 * setdiff(1:13, 5)))); expect_equal(row14$x, log1p(14000)); expect_equal(row14$n_before, 12L)
  X2 <- stm_days(make_days(function(j) ifelse(j >= 14, 99999, 1000 * j)), sc)        # change the steps of every day from day 14 on
  expect_equal(X2[date == d0 + 13L, m_before], row14$m_before)                       # the baseline of day 14 is unchanged
})

test_that("recording history on these day statuses", {
  X <- stm_days(make_days(), stm_config())
  r <- X[date == d0 + 16L]                                                           # day 17: day 15 was a spell of one day, returned on day 16
  expect_equal(c(r$any_prior, r$prev_spell, r$dsr), c(TRUE, 1L, 1L))
  r0 <- X[date == d0 + 11L]; expect_equal(c(r0$any_prior, r0$dsr), c(FALSE, 0L))
})

test_that("event counts, thresholds, variants and tests", {
  sc <- stm_config()
  R <- data.table(person_id = rep(1:2, each = 6), cohort = factor("Healthy"), outcome = c(0, 0, 1, 0, 3, 0, 0, 0, 0, 2, 0, 0), study_day = rep(c(5L, 10L, 15L, 20L, 25L, 30L), 2))
  ec <- stm_event_counts(R, sc)[group == "all"]
  expect_equal(ec$events, c(9L, 1L, 1L, 1L, 0L)); expect_false(any(ec$effect_claim_allowed))                           # far below the thresholds
  expect_equal(nrow(stm_variant_rows(R, "first_transition", sc)), 3L + 4L)                                              # up to and including each person's first outcome other than 0
  expect_equal(nrow(stm_variant_rows(R, "first14", sc)), 8L)
  coef <- data.table(model = c("c1", "c2", "c3"), term = "w", estimate = log(c(0.8, 1.001, 0.9)), se = c(0.01, 0.005, 0.02), p = 0.5, or = c(0.8, 1.001, 0.9),
                     or_low = 0, or_high = 0, or90_low = exp(log(c(0.8, 1.001, 0.9)) - 1.645 * c(0.01, 0.005, 0.02)), or90_high = exp(log(c(0.8, 1.001, 0.9)) + 1.645 * c(0.01, 0.005, 0.02)))
  info <- data.table(model = c("c1", "c2", "c3"), n_events = c(2000L, 7000L, 7000L))
  t <- stm_tests(coef, info, sc)
  expect_equal(t$equivalent, c(FALSE, TRUE, FALSE)); expect_equal(t$effect_claim_allowed, c(TRUE, TRUE, TRUE)); expect_equal(t$equivalence_claim_allowed, c(FALSE, TRUE, TRUE))
  expect_lt(t[model == "c2", p_tost], 0.05); expect_gt(t[model == "c3", p_tost], 0.05)
})

test_that("the multinomial model equals nnet::multinom and its cluster covariance equals the hand-made sandwich", {
  set.seed(8); n <- 3000; g <- 150
  dat <- data.table(person_id = sample(seq_len(g), n, TRUE), w = rnorm(n), z = rnorm(n))
  dat[, outcome := { eta <- cbind(-1.5 - 0.4 * w, -1.8 + 0.2 * w, -1.2 - 0.1 * z, -2 + 0.1 * w + 0.3 * z); pr <- cbind(1, exp(eta)); pr <- pr / rowSums(pr)
                     apply(pr, 1L, function(p) sample(0:4, 1L, prob = p)) }]
  fit <- stm_multinom(dat, c("w", "z")); ref <- nnet::multinom(factor(outcome) ~ w + z, dat, trace = FALSE, reltol = 1e-12, maxit = 500)
  expect_true(fit$converged); expect_equal(unname(t(fit$beta)), unname(coef(ref)), tolerance = 1e-4)
  X <- cbind(1, dat$w, dat$z); P <- predict(ref, type = "probs")[, -1L]; Y <- vapply(1:4, function(k) as.numeric(dat$outcome == k), numeric(n))
  S <- do.call(cbind, lapply(1:4, function(k) rowsum(X * (Y[, k] - P[, k]), dat$person_id)))
  bread <- unname(vcov(ref)[as.vector(outer(1:3, 0:3 * 3, `+`)), as.vector(outer(1:3, 0:3 * 3, `+`))])
  expect_equal(unname(fit$V), unname((g / (g - 1)) * bread %*% crossprod(S) %*% bread), tolerance = 1e-3)
  ct <- stm_contrasts(fit, list(c(1L, 2L), c(1L, 3L))); expect_equal(ct$difference[1L], fit$beta["w", 1L] - fit$beta["w", 2L]); expect_equal(ct$p_holm, stats::p.adjust(ct$p, "holm"))
})

test_that("a simulated study end to end: units resume, tables and the three views are written", {
  cfg <- cfg_test(); st <- st_config(); sc <- stm_config(list(n_profile_rows = 200L, n_beta_draws = 100L, n_re_grid = 11L))
  inp <- suppressMessages(simulate_inputs(cfg, n_per_cohort = c(Insomnia = 60, Diabetes = 60, Healthy = 120, Prostate = 40), seed = 5))
  p <- copy(inp$persons)[, device_family := "Charge"]
  D <- st_build_days(as_day_store(inp$day_tables, inp$persons), p, st_trackers(sim_device_rows(inp$persons, cfg), p, cfg), st, cfg)
  rs <- stm_risk_set(stm_days(D, sc), p, cfg, sc); ck <- tempfile()
  out <- suppressWarnings(suppressMessages(stm_run_models(rs$R, "main", cfg, sc, ck, list(k = 1), secondary = TRUE, multinomial = TRUE)))
  expect_equal(length(list.files(ck, pattern = "^streammodel__main__.*[.]rds$")), 9L)
  t0 <- Sys.time(); again <- suppressMessages(stm_run_models(rs$R, "main", cfg, sc, ck, list(k = 1), secondary = TRUE, multinomial = TRUE)); expect_lt(as.numeric(difftime(Sys.time(), t0, units = "secs")), 5)
  res <- stm_collect(again, sc, "main"); expect_true(all(c("c1", "c2", "c3", "c4") %in% res$info$model)); expect_true(all(is.finite(res$risk_curve$prob_pct)))
  expect_true(nrow(res$contrasts) == 2L && all(c("c0", "c4") %in% res$risk_curve[source == "multinomial", model]))
  figs <- stm_figures(res, stm_deciles(rs$R, sc), sc, cbind(stm_event_counts(rs$R, sc), variant = "main")); expect_equal(length(figs), 5L)
  d <- tempfile(); stm_write(c(res, list(event_counts = stm_event_counts(rs$R, sc))), figs, d)
  expect_true(all(file.exists(file.path(d, "figures", paste0(names(figs), ".png"))))); expect_true(file.exists(file.path(d, "tables", "stream_model_event_counts.tsv")))
})
