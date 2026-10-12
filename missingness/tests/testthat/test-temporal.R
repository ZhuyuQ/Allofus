test_that("weekday counts equal a brute-force count", {
  set.seed(2)
  first <- as.Date("2023-03-01") + sample(0:60, 40, TRUE); W <- sample(1:100, 40, TRUE)
  wc <- window_weekday_counts(first, W)
  for (i in seq_along(first)) {
    truth <- table(factor(DOW_LABELS[wday(seq(first[i], by = "day", length.out = W[i]))], levels = DOW_LABELS))
    got <- wc[row == i]
    expect_equal(as.integer(got$den), as.integer(truth))
  }
})

make_pg <- function(days, persons_cohort = "Healthy") {
  g <- get_gap_runs(days); pg <- summarise_person_gaps(days, g); pg[, `:=`(cohort = factor(persons_cohort), stream = "activity")]
  list(pg = pg, g = g)
}

test_that("day-of-week, study-week and calendar-year counts add up to the window and the gaps", {
  days <- random_days(20, seed = 3); x <- make_pg(days)
  dw <- dow_counts(x$pg, x$g); wk <- study_week_counts(x$pg, x$g); yr <- calendar_year_counts(x$pg, x$g)
  tot <- x$pg[, .(win = sum(window_days), mis = sum(n_missing))]
  expect_equal(sum(dw$den), tot$win); expect_equal(sum(dw$num), tot$mis)
  expect_equal(sum(wk$den), tot$win); expect_equal(sum(wk$num), tot$mis)
  expect_equal(sum(yr$den), tot$win); expect_equal(sum(yr$num), tot$mis)
})

test_that("transition counts: continuation days + onsets = missing days", {
  days <- random_days(20, seed = 4); x <- make_pg(days)
  tc <- transition_person_counts(x$pg)
  expect_true(all(tc$continue_days + tc$onsets == tc$miss_days))
  tr <- transition_table(tc, cfg_test())
  expect_true(all(tr$p_onset >= 0 & tr$p_onset <= 1))
  expect_equal(tr$p_onset, sum(tc$onsets) / sum(tc$at_risk))
})

# ---- attained age -------------------------------------------------------------------------------------------

age_fixture <- function(seed = 3, n = c(Insomnia = 15, Diabetes = 15, Healthy = 30, Prostate = 10)) {
  cfg <- cfg_test(); inp <- suppressMessages(simulate_inputs(cfg, n_per_cohort = n, seed = seed))
  d <- inp$day_tables$activity; gaps <- get_gap_runs(d)
  pg <- with_cohort(summarise_person_gaps(d, gaps), inp$persons, "activity")
  list(cfg = cfg, persons = inp$persons, days = d, gaps = gaps, pg = pg)
}

test_that("age of a day is exact at the year boundaries and the same thresholds are used for both counts", {
  age0 <- 60.5                                    # = 60 years at first record + 0.5
  e <- c(0L, 179L, 180L, 181L, 364L, 365L, 546L, 547L)
  expect_equal(age_year_of_day(age0, e), floor(age0 + e / 365.25))
  expect_equal(age_year_of_day(age0, 0:3000), floor(age0 + (0:3000) / 365.25))
})

test_that("window days and missing days are split over ages without loss, and agree with a day-by-day count", {
  fx <- age_fixture()
  ac <- age_counts(fx$pg, fx$gaps, fx$persons)
  tot <- ac[, .(den = sum(den), num = sum(num)), by = person_id]
  chk <- merge(tot, fx$pg[, .(person_id, window_days, n_missing)], by = "person_id")
  expect_equal(chk$den, chk$window_days); expect_equal(chk$num, chk$n_missing)
  expect_true(all(ac$num <= ac$den))
  # brute force for one person: walk every calendar day of the window
  pid <- fx$pg$person_id[fx$pg$n_missing > 0][1]; p <- fx$pg[person_id == pid]
  age_first <- fx$persons[person_id == pid, age_first]
  day <- seq(p$first_date, p$last_date, by = "day"); e <- as.integer(day - p$first_date)
  observed <- day %in% fx$days[person_id == pid, date]
  a <- floor(age_first + 0.5 + e / 365.25)
  brute <- data.table(age = a, miss = !observed)[, .(den = .N, num = sum(miss)), by = age][order(age)]
  got <- ac[person_id == pid][order(age)]
  expect_equal(got$age, brute$age); expect_equal(got$den, brute$den); expect_equal(got$num, brute$num)
})

test_that("the moving window pools 2h+1 neighbouring ages and keeps the totals consistent", {
  fx <- age_fixture(); h <- 2L
  a0 <- age_counts(fx$pg, fx$gaps, fx$persons, 0L); ah <- age_counts(fx$pg, fx$gaps, fx$persons, h)
  expect_equal(sum(ah$den), (2 * h + 1) * sum(a0$den)); expect_equal(sum(ah$num), (2 * h + 1) * sum(a0$num))
  pid <- a0$person_id[1]; x <- a0[person_id == pid]; g <- ah[person_id == pid]
  centre <- x$age[1]
  expect_equal(g[age == centre, den], sum(x[age >= centre - h & age <= centre + h, den]))
  expect_true(all(g$num <= g$den))
})

test_that("rate tables by age follow the same layout as the other rate tables and need >= 30 people per age", {
  fx <- age_fixture(n = c(Insomnia = 120, Diabetes = 120, Healthy = 240, Prostate = 60))
  r <- missing_rate_tables(fx$pg, list(activity = fx$gaps), fx$persons, "activity", fx$cfg)
  expect_setequal(names(r), c("study_week_rates", "calendar_year_rates", "weekday_rates", "age_rates", "missing_trend_tests",
                              "missing_trend_heterogeneity", "missing_weekday_tests", "missing_weekday_odds_ratios"))
  expect_true(all(c("stream", "cohort", "age", "estimate", "lo", "hi", "n_persons") %in% names(r$age_rates)))
  expect_true(all(r$age_rates$n_persons >= 30L)); expect_true(all(r$age_rates$lo <= r$age_rates$estimate & r$age_rates$estimate <= r$age_rates$hi))
})

# ---- statistics behind the "missing by ..." figures ------------------------------------------------------------------

planted_counts <- function(slopes, n = 150, units = 1:30, seed = 1, person_sd = 0.4, daily = 7L) {
  set.seed(seed)
  rbindlist(lapply(seq_along(slopes), function(k) rbindlist(lapply(seq_len(n), function(i) {
    u <- rnorm(1, 0, person_sd)
    data.table(person_id = k * 1000L + i, cohort = names(slopes)[k], week = units,
               den = daily, num = rbinom(length(units), daily, plogis(-2.5 + u + slopes[[k]] * (units - mean(units)))))
  }))))
}

test_that("Wald chi-square test of a linear hypothesis", {
  w <- wald_chisq(c(1, 2), diag(c(1, 4)), diag(2)); expect_equal(w$chisq, 2); expect_equal(w$df, 2); expect_equal(w$p, pchisq(2, 2, lower.tail = FALSE))
  expect_true(is.na(wald_chisq(c(1, 2), matrix(0, 2, 2), diag(2))$chisq))                            # singular covariance -> NA, not an error
})

test_that("clustered logistic regression gives the binomial estimates with person-clustered standard errors", {
  d <- planted_counts(list(A = 0.03), n = 200, person_sd = 0.8); d[, xs := week - mean(week)]
  m <- cluster_logit(d, "xs"); fit <- glm(cbind(num, den - num) ~ xs, data = d, family = binomial())
  expect_equal(unname(m$b), unname(coef(fit))); expect_equal(m$V, sandwich::vcovCL(fit, cluster = d$person_id), ignore_attr = TRUE)
  dc <- copy(d)[, person_level := (person_id %% 7) / 7]                                              # a covariate constant within person: clustering must widen the SE
  mc <- cluster_logit(dc, "person_level"); fc <- glm(cbind(num, den - num) ~ person_level, data = dc, family = binomial())
  expect_gt(sqrt(mc$V["person_level", "person_level"]), sqrt(vcov(fc)["person_level", "person_level"]))
})

test_that("trend statistics recover planted slopes, and the heterogeneity test separates different from equal slopes", {
  # slope 0.05 per week = 0.2 per 4 weeks; the second cohort is flat
  wk <- planted_counts(list(Insomnia = 0.05, Diabetes = 0.0, Healthy = 0.05, Prostate = 0.05))[, stream := "activity"]
  cn <- list(week = wk, year = wk[, .(stream, person_id, cohort, year = week, num, den)][0], age = wk[0], dow = wk[0])
  cn$year <- copy(wk)[, year := week][, week := NULL]; cn$age <- copy(wk)[, age := week][, week := NULL]
  r <- missing_trend_statistics(cn)
  tr <- r$missing_trend_tests[variable == "study_week"]
  expect_equal(tr[cohort == "Insomnia", odds_ratio], exp(0.2), tolerance = 0.1); expect_lt(tr[cohort == "Insomnia", p], 1e-6)
  expect_gt(tr[cohort == "Diabetes", p], 0.01); expect_lt(abs(log(tr[cohort == "Diabetes", odds_ratio])), 0.15)
  expect_true(all(tr[!grepl("^All", cohort), p_holm] >= tr[!grepl("^All", cohort), p]))
  expect_lt(r$missing_trend_heterogeneity[variable == "study_week", p], 1e-6)                          # the flat cohort differs from the others
  eq <- planted_counts(list(Insomnia = 0.05, Diabetes = 0.05, Healthy = 0.05, Prostate = 0.05), seed = 4)[, stream := "activity"]
  cn2 <- list(week = eq, year = copy(eq)[, year := week][, week := NULL], age = copy(eq)[, age := week][, week := NULL], dow = eq[0])
  expect_gt(missing_trend_statistics(cn2)$missing_trend_heterogeneity[variable == "study_week", p], 0.01)
  # too few distinct units or people -> reported as not estimable, not as a result
  tiny <- list(week = eq[person_id %in% (1000L + 1:5)], year = cn2$year[0], age = cn2$age[0], dow = eq[0])
  expect_true(all(missing_trend_statistics(tiny)$missing_trend_tests[variable == "study_week" & !grepl("^All", cohort), status] == "not_estimable"))
})

test_that("weekday statistics find a planted weekday effect and report odds ratios versus Sunday", {
  set.seed(3)
  dw <- rbindlist(lapply(c("Insomnia", "Healthy"), function(co) rbindlist(lapply(1:200, function(i) {
    u <- rnorm(1, 0, 0.3)
    data.table(person_id = (co == "Healthy") * 1000L + i, cohort = co, dow = DOW_LABELS, den = 20L,
               num = rbinom(7, 20, plogis(-2.5 + u + c(0, 0, 0, 0, 0, 0, if (co == "Insomnia") 0.5 else 0))))   # only Insomnia is worse on Saturdays
  }))))[, stream := "activity"]
  r <- missing_weekday_statistics(list(dow = dw))
  t <- r$missing_weekday_tests
  expect_lt(t[cohort == "Insomnia" & test == "weekday effect", p], 1e-6); expect_gt(t[cohort == "Healthy" & test == "weekday effect", p], 0.01)
  expect_lt(t[test == "weekday x cohort", p], 0.001)
  o <- r$missing_weekday_odds_ratios[cohort == "Insomnia"]
  expect_equal(o[weekday == "Sat", odds_ratio_vs_sunday], exp(0.5), tolerance = 0.15); expect_setequal(o$weekday, DOW_LABELS[-1L])
  expect_true(all(o$or_low < o$odds_ratio_vs_sunday & o$odds_ratio_vs_sunday < o$or_high))
})
