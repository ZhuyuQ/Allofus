test_that("FULL-PERIOD standardization (sensitivity): within-person z has pooled SD 1 and between-person z has mean 0, SD 1", {
  cfg <- cfg_test(); set.seed(1)
  days <- data.table(person_id = rep(1:30, each = 40), date = rep(as.Date("2024-01-01") + 0:39, 30),
                     steps = rpois(1200, rep(runif(30, 3000, 12000), each = 40)), sedentary_minutes = 800, mvpa_minutes = 10)
  setkey(days, person_id, date)
  pt <- predictor_table(days, "act_steps", cfg, standardization = "full_period")
  w <- pt$w_act_steps
  expect_equal(sqrt(sum(w^2) / (length(w) - 30)), 1, tolerance = 1e-9)
  expect_equal(max(abs(pt[, mean(w_act_steps), by = person_id]$V1)), 0, tolerance = 1e-9)
  b <- unique(pt[, .(person_id, b_act_steps)])$b_act_steps
  expect_equal(mean(b), 0, tolerance = 1e-9); expect_equal(sd(b), 1, tolerance = 1e-9)
})

test_that("trailing baseline uses only previous observations and needs >= 7 of them", {
  cfg <- cfg_test()
  days <- data.table(person_id = 1, date = as.Date("2024-01-01") + 0:19, steps = 1:20, sedentary_minutes = 1, mvpa_minutes = 1)
  setkey(days, person_id, date)
  pt <- predictor_table(days, "act_steps", cfg)
  sd_w <- attr(pt, "scales")$sd_within
  x <- log1p(1:20); tb <- x - pt$tw_act_steps * sd_w
  expect_true(all(is.na(tb[1:7])))                                # fewer than 7 earlier days
  expect_equal(tb[8], mean(x[1:7]), tolerance = 1e-9)
  expect_equal(tb[10], mean(x[1:9]), tolerance = 1e-9)
  expect_equal(tb[20], mean(x[6:19]), tolerance = 1e-9)           # window of 14, excluding today
})

test_that("partial-wear rule evaluates and tolerates a missing column", {
  d <- data.table(steps = c(100, 5000), sedentary_minutes = c(1400, 600))
  expect_equal(partial_wear_flag(d, "steps < 500 | sedentary_minutes >= 1380"), c(TRUE, FALSE))
  expect_error(partial_wear_flag(d, "no_such_col < 1"), "cannot be evaluated")      # a bad rule must not silently flag nothing
})

test_that("trailing_mean_excl_current ignores NAs, excludes today and respects the minimum count", {
  x <- c(1, 2, NA, 4, 5, 6, 7, 8, 9, 10)
  got <- trailing_mean_excl_current(x, n = 4, min_n = 3)
  expect_true(all(is.na(got[1:3])))                       # fewer than 3 non-missing values so far
  expect_equal(got[5], mean(c(1, 2, 4)))                  # positions 1-4, the NA skipped
  expect_equal(got[8], mean(x[4:7]))                      # the 4 positions before today
  expect_equal(got[10], mean(x[6:9]))
  expect_equal(trailing_mean_excl_current(1:20, 14, 7)[20], mean(6:19))
})

test_that("risk-set rows pair day d's own predictor with the outcome on day d + 1 (no look-ahead in the pairing)", {
  cfg <- cfg_test()
  inp <- suppressMessages(simulate_inputs(cfg, n_per_cohort = c(Insomnia = 20, Diabetes = 20, Healthy = 40, Prostate = 10), seed = 8))
  b <- suppressMessages(prepare_stream_bundle("activity", inp$day_tables, inp$persons, cfg, trailing = TRUE))
  days <- inp$day_tables$activity
  # independent construction on an explicit calendar: is tomorrow observed?
  obs <- days[, .(person_id, date, obs = TRUE)]
  chk <- merge(b$rs[, .(person_id, date, y_ge1, w_act_steps)], obs[, .(person_id, tomorrow = date - 1L, next_seen = TRUE)],
               by.x = c("person_id", "date"), by.y = c("person_id", "tomorrow"), all.x = TRUE)
  expect_equal(chk$y_ge1, as.integer(is.na(chk$next_seen)))                    # gap starts tomorrow <=> tomorrow not observed
  # the predictor on row (p, d) is day d's own value, standardised with person-level constants
  # (primary = pre-gap standardization: the day's own value against the mean of the person's EARLIER observed days)
  x <- log1p(pmax(days$steps, 0)); d2 <- copy(days)[, x := x][, m := prior_mean(x, cfg$standardization$baseline_min_n), by = person_id]
  sdw <- b$scales[predictor == "act_steps", sd_within]
  d2[, w := (x - m) / sdw]
  j <- merge(b$rs[, .(person_id, date, w_act_steps)], d2[, .(person_id, date, w)], by = c("person_id", "date"))
  expect_equal(j$w_act_steps, j$w)
  # (full-period standardization: person-level constants over the whole period)
  bf <- suppressMessages(prepare_stream_bundle("activity", inp$day_tables, inp$persons, cfg, standardization = "full_period"))
  d3 <- copy(days)[, x := x][, m := mean(x, na.rm = TRUE), by = person_id]
  d3[, w := (x - m) / bf$scales[predictor == "act_steps", sd_within]]
  j3 <- merge(bf$rs[, .(person_id, date, w_act_steps)], d3[, .(person_id, date, w)], by = c("person_id", "date"))
  expect_equal(j3$w_act_steps, j3$w)
})

# ---- pre-gap standardization ---------------------------------------------------------------------------------------

pg_days <- function(n_person = 25, n_day = 60, seed = 3) {
  set.seed(seed)
  d <- data.table(person_id = rep(seq_len(n_person), each = n_day), date = rep(as.Date("2024-01-01") + seq_len(n_day) - 1L, n_person))
  d[, steps := pmax(100, round(rnorm(.N, rep(runif(n_person, 3000, 12000), each = n_day), 1500)))]
  d[, `:=`(sedentary_minutes = 800, mvpa_minutes = 10)]
  setkey(d, person_id, date); d
}

test_that("prior_mean uses only EARLIER observations, skips NAs and needs min_n of them; through_mean includes today", {
  x <- c(1, 2, NA, 4, 10, 6)
  p <- prior_mean(x, 2)
  expect_true(all(is.na(p[1:2])))                           # fewer than 2 earlier observations
  expect_equal(p[3], mean(c(1, 2))); expect_equal(p[4], mean(c(1, 2))); expect_equal(p[5], mean(c(1, 2, 4))); expect_equal(p[6], mean(c(1, 2, 4, 10)))
  th <- through_mean(x, 2)
  expect_true(is.na(th[1])); expect_equal(th[2], 1.5); expect_equal(th[4], mean(c(1, 2, 4))); expect_equal(th[6], mean(c(1, 2, 4, 10, 6)))
})

test_that("pre-gap standardization: x = m + sd_w * w and m = mu_b + sd_b * b, the constants come from the rows that have a baseline", {
  cfg <- cfg_test(); days <- pg_days(); n_min <- cfg$standardization$baseline_min_n
  pt <- predictor_table(days, "act_steps", cfg, trailing = FALSE, keep_x = TRUE, keep_baseline = FALSE)
  sc <- attr(pt, "scales"); expect_equal(sc$standardization, "pre_gap")
  x <- pt$x_act_steps; m <- days[, prior_mean(log1p(steps), n_min), by = person_id]$V1
  expect_equal(x - sc$sd_within * pt$w_act_steps, m)                         # the within part is the deviation from the prior usual
  ok <- !is.na(m)
  expect_equal((m - sc$sd_between * pt$b_act_steps)[ok], rep(sc$mean_of_person_means, sum(ok)))   # and the usual is mu_b + sd_b * b
  expect_equal(sc$sd_within, sqrt(mean((x[ok] - m[ok])^2)))                  # reference distributions: rows with a baseline
  expect_equal(c(mean(pt$b_act_steps[ok]), sd(pt$b_act_steps[ok])), c(0, 1), tolerance = 1e-9)
  expect_true(all(is.na(pt$w_act_steps[days[, seq_len(.N) <= n_min, by = person_id]$V1])))   # the first days have no baseline
  expect_true(all(is.na(pt$b_act_steps) == is.na(pt$w_act_steps)))
})

test_that("NO LOOK-AHEAD: changing what a person does LATER never changes the person's earlier deviation from their usual", {
  cfg <- cfg_test(); days <- pg_days(); n_min <- cfg$standardization$baseline_min_n
  dev_of <- function(d) { pt <- predictor_table(d, "act_steps", cfg, trailing = FALSE, keep_x = FALSE); pt$w_act_steps * attr(pt, "scales")$sd_within }
  d1 <- copy(days); d2 <- copy(days)
  later <- d2$person_id == 7L & d2$date >= as.Date("2024-02-01")
  d2[later, steps := steps * 3]                                              # person 7 changes from Feb 1 on (a "gap-period" or recovery)
  a <- dev_of(d1); b <- dev_of(d2)
  before <- days$person_id == 7L & days$date < as.Date("2024-02-01")
  expect_equal(a[before], b[before], tolerance = 1e-12)                       # unchanged before the change (as a deviation in log-step units)
  expect_true(any(abs(a[days$person_id == 7L & !before] - b[days$person_id == 7L & !before]) > 1e-6, na.rm = TRUE))   # and it does change afterwards
  # the full-period standardization DOES leak the future (documented; it is the sensitivity analysis)
  fp <- function(d) { pt <- predictor_table(d, "act_steps", cfg, trailing = FALSE, keep_x = FALSE, standardization = "full_period"); pt$w_act_steps * attr(pt, "scales")$sd_within }
  expect_true(any(abs(fp(d1)[before] - fp(d2)[before]) > 1e-6, na.rm = TRUE))
})

test_that("the trailing baseline uses the pre-gap scale, and effect sizes mean the same thing in every standardization", {
  cfg <- cfg_test(); days <- pg_days()
  pt <- predictor_table(days, "act_steps", cfg, trailing = TRUE, keep_x = TRUE)
  sdw <- attr(pt, "scales")$sd_within
  tb <- days[, trailing_mean_excl_current(log1p(steps), cfg$onset$trailing_n, cfg$onset$trailing_min_n), by = person_id]$V1
  expect_equal(pt$tw_act_steps, (pt$x_act_steps - tb) / sdw)
  # one SD of w is one within-person SD (about the usual) in both; the values differ because the usual and the SD differ
  pf <- predictor_table(days, "act_steps", cfg, trailing = FALSE, standardization = "full_period")
  expect_gt(cor(pt$w_act_steps, pf$w_act_steps, use = "complete.obs"), 0.9)
  expect_error(predictor_table(days, "act_steps", cfg, standardization = "nope"))
})
