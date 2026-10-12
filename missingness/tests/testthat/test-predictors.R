test_that("within-person z has pooled SD 1 and between-person z has mean 0, SD 1", {
  cfg <- cfg_test(); set.seed(1)
  days <- data.table(person_id = rep(1:30, each = 40), date = rep(as.Date("2024-01-01") + 0:39, 30),
                     steps = rpois(1200, rep(runif(30, 3000, 12000), each = 40)), sedentary_minutes = 800, mvpa_minutes = 10)
  setkey(days, person_id, date)
  pt <- predictor_table(days, "act_steps", cfg)
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
  x <- log1p(pmax(days$steps, 0)); d2 <- copy(days)[, x := x][, m := mean(x, na.rm = TRUE), by = person_id]
  sdw <- b$scales[predictor == "act_steps", sd_within]
  d2[, w := (x - m) / sdw]
  j <- merge(b$rs[, .(person_id, date, w_act_steps)], d2[, .(person_id, date, w)], by = c("person_id", "date"))
  expect_equal(j$w_act_steps, j$w)
})
