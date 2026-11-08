# Weekly rhythm of gap starts and gap ends (R/weekday_rhythm.R, R/plots_weekday_rhythm.R).

wr_cfg <- function() {
  cfg <- cfg_test()
  cfg$weekday_rhythm[c("n_profile_rows", "n_beta_draws", "n_re_draws", "n_curve_draws", "n_person_draws")] <- list(200L, 30L, 15L, 15L, 1000L)
  cfg
}

test_that("the gap-start weekday equals the onset model's dow_t; the sleep night-start shift moves it by exactly one day and nothing else", {
  cfg <- wr_cfg(); sim <- simulate_gap_process(60, 120, seed = 2)
  a <- wr_gap_start_rows(sim$day_tables$activity, sim$persons, "activity", cfg)
  rs <- build_risk_set(sim$day_tables$activity, cfg)
  m <- merge(a[, .(person_id, date, wd)], rs[, .(person_id, date, dow_t)], by = c("person_id", "date"))
  expect_equal(m$wd, as.integer(m$dow_t) - 1L)
  s <- wr_gap_start_rows(sim$day_tables$sleep, sim$persons, "sleep", cfg)
  expect_equal(s$wd, (a$wd - 1L) %% 7L)                                         # the night started the day before the stored date
  same <- c("person_id", "date", "start", "dsr", "study_day", "prev_gap_len_l", "any_prior_gap")
  expect_equal(s[, ..same], a[, ..same])
  cfg$weekday_rhythm$sleep_alignment <- "as_stored"
  expect_equal(wr_gap_start_rows(sim$day_tables$sleep, sim$persons, "sleep", cfg)$wd, a$wd)
  expect_equal(wr_weekday(as.Date("2024-03-04"), "activity", cfg), 0L); expect_equal(wr_weekday(as.Date("2024-03-10"), "activity", cfg), 6L)   # Monday, Sunday
})

test_that("days since return agrees with hand counts, including before the first gap", {
  cfg <- wr_cfg(); d0 <- as.Date("2024-01-01")
  days <- data.table(person_id = 1L, date = d0 + c(0:4, 7:11, 13:15) - 0L)       # observed 1-5, gap 6-7, observed 8-12, gap 13, observed 14-16
  persons <- data.table(person_id = 1L, cohort = "Healthy", sex = "Female", age_first = 50, cci_score = 1, cci_group = "1-2")
  r <- wr_gap_start_rows(days, persons, "activity", cfg)
  expect_equal(r$dsr, c(1:5, 1:5, 1:2)); expect_equal(r$start, c(0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0))
  expect_equal(as.character(r$dsr_f[c(5, 6, 12)]), c("5", "1", "2"))
  expect_equal(nrow(r), 12L)                                                    # the person's last observed day has no next day
})

test_that("gap-end rows: one row per gap day, the outcome only on the return day, 14 days so far dropped and counted", {
  cfg <- wr_cfg(); d0 <- as.Date("2024-01-01")
  days <- data.table(person_id = 1L, date = d0 + c(0:2, 6:9, 24:26))             # gaps of 3 days (indices 3-5) and 14 days (10-23)
  persons <- data.table(person_id = 1L, cohort = "Healthy", sex = "Female", age_first = 50, cci_score = 1, cci_group = "1-2")
  g <- get_gap_runs(days); expect_equal(g$gap_length, c(3L, 14L))
  e <- wr_gap_end_rows(days, persons, "activity", cfg)
  expect_equal(e[, .N, by = gap_length]$N, c(3L, 13L))                          # 3 rows, and 13 of the 14 (14 days so far dropped)
  expect_equal(e$end, c(0, 0, 1, rep(0, 13)))                                   # the 14-day gap's only possible end (14 days so far) is the dropped row
  expect_equal(e$dsf[1:3], 1:3)
  expect_equal(e$out_date[3], d0 + 6)                                           # the return day
  fl <- attr(e, "flow"); expect_equal(fl$n_rows, c(17L, 1L))
})

test_that("US federal holidays are built by rule: known dates 2019-2024 including observed dates", {
  h <- wr_us_holidays(2019:2024)
  known <- as.Date(c("2019-11-28", "2020-11-26", "2021-11-25", "2022-11-24", "2023-11-23", "2024-11-28",      # Thanksgiving
                     "2019-05-27", "2020-05-25", "2021-05-31", "2022-05-30", "2023-05-29", "2024-05-27",      # Memorial Day
                     "2019-09-02", "2020-09-07", "2021-09-06", "2022-09-05", "2023-09-04", "2024-09-02",      # Labor Day
                     "2019-01-21", "2020-01-20", "2021-01-18", "2022-01-17", "2023-01-16", "2024-01-15",      # Martin Luther King Jr. Day
                     "2019-02-18", "2020-02-17", "2021-02-15", "2022-02-21", "2023-02-20", "2024-02-19",      # Presidents Day
                     "2019-10-14", "2020-10-12", "2021-10-11", "2022-10-10", "2023-10-09", "2024-10-14"))     # Columbus Day
  expect_true(all(known %in% h$date))
  obs <- as.Date(c("2020-07-03", "2021-07-05", "2021-12-24", "2022-12-26", "2021-12-31", "2023-01-02", "2021-06-18", "2023-11-10"))   # Saturday -> Friday, Sunday -> Monday
  expect_true(all(obs %in% h$observed))
  expect_false(any(h[format(date, "%Y") == "2020"]$holiday == "Juneteenth"))    # from 2021
  expect_equal(wr_near_holiday(as.Date(c("2024-11-27", "2024-11-29", "2024-11-30", "2024-12-02")), as.Date("2024-11-28"), 1L), c(TRUE, TRUE, FALSE, FALSE))
})

test_that("the first-harmonic projection recovers a planted rhythm; ratios to the average day and peak-to-trough match hand computation", {
  k <- 0:6; a <- 0.4; b <- -0.2; e <- a * cos(2 * pi * k / 7) + b * sin(2 * pi * k / 7)
  beta <- setNames(e[1:6], paste0("weekday", 1:6))
  h <- wr_first_harmonic(beta); expect_equal(unname(h), c(a, b), tolerance = 1e-12)
  expect_equal(wr_peak_day(a, b), (atan2(b, a) * 7 / (2 * pi)) %% 7)
  P <- cbind(c(0.02, 0.03, 0.025, 0.02, 0.01, 0.04, 0.03), c(0.021, 0.031, 0.024, 0.02, 0.011, 0.039, 0.03), c(0.019, 0.029, 0.026, 0.02, 0.009, 0.041, 0.03))
  r <- wr_profile_from_matrix(P)
  expect_equal(r$profile$ratio_to_average_day, P[, 1] / mean(P[, 1])); expect_equal(r$summary$peak_to_trough, max(P[, 1]) / min(P[, 1]))
  expect_equal(c(r$summary$peak_weekday, r$summary$trough_weekday), c(5L, 4L)); expect_equal(r$profile$pct, 100 * P[, 1])
})

test_that("Holm is applied within each family only, H2r is not adjusted, and tests with different random-effects structures are not made", {
  mk <- function(s, ll) rbindlist(lapply(c("S0", "S1", "S2", "S3", "S4", "S5", "E0", "E1", "E4"), function(f)
    data.table(stream = s, fit_id = f, loglik = ll[[f]], df = c(S0 = 10L, S1 = 16L, S2 = 31L, S3 = 37L, S4 = 34L, S5 = 12L, E0 = 20L, E1 = 26L, E4 = 44L)[[f]], re_structure = "independent_harmonics")))
  base <- c(S0 = -1000, S1 = -997, S2 = -990, S3 = -986, S4 = -990, S5 = -1000, E0 = -500, E1 = -497, E4 = -490)
  fits <- rbind(mk("activity", base), mk("sleep", base - c(0, 4, 0, 1, 7, 0, 0, 5, 12)), mk("hr", base - 1))
  t <- wr_tests(fits)
  for (fam in c("H1", "H2", "H3", "H4")) expect_equal(t[family == fam, p_holm], stats::p.adjust(t[family == fam, p], "holm"))
  expect_equal(t[family == "H1", .N], 6L); expect_equal(t[family == "H4", .N], 6L); expect_equal(t[family == "H2", .N], 3L)
  expect_true(all(is.na(t[id == "H2r", p_holm])))
  expect_equal(t[id == "H1s" & stream == "activity", df], 6L); expect_equal(t[id == "H4s" & stream == "activity", df], 18L)
  fits[stream == "activity" & fit_id == "S0", re_structure := "correlated"]
  t2 <- wr_tests(fits); expect_true(is.na(t2[id == "H1s" & stream == "activity", p])); expect_match(t2[id == "H1s" & stream == "activity", note], "differ")
})

test_that("cells with 1-19 people are suppressed, cells with few days are flagged", {
  cfg <- wr_cfg()
  tab <- data.table(n_people = c(5L, 25L, 0L), n_days = c(1000L, 150L, 0L), rate_pct = c(1, 2, NA), n_events = c(3L, 4L, 0L))
  s <- wr_suppress(tab, c("rate_pct", "n_events"), cfg, "n_days")
  expect_equal(s$people_label, c("<20", "25", "0")); expect_true(is.na(s$rate_pct[1]) && is.na(s$n_events[1]) && is.na(s$n_people[1]))
  expect_equal(s$suppressed, c(TRUE, FALSE, FALSE)); expect_equal(s$low_exposure, c(FALSE, TRUE, TRUE))
  expect_equal(s$rate_pct[2], 2)
})

test_that("the stage: fits share rows and random-effects structure within a test, checkpoints resume and change with the settings, every figure builds", {
  cfg <- wr_cfg(); inp <- suppressMessages(simulate_inputs(cfg, n_per_cohort = c(Insomnia = 40, Diabetes = 40, Healthy = 80, Prostate = 30), seed = 5))
  store <- as_day_store(inp$day_tables, inp$persons); ck <- tempfile(); dk <- data_key(store, "activity")
  run <- function(cfg) {
    days <- lazy(function() store$get("activity", character()))
    sr <- lazy(function() wr_gap_start_rows(days$get(), inp$persons, "activity", cfg)); er <- lazy(function() wr_gap_end_rows(days$get(), inp$persons, "activity", cfg))
    b <- lazy(function() suppressMessages(prepare_stream_bundle("activity", stream_inputs(store, "activity", cfg), inp$persons, cfg)))
    suppressWarnings(suppressMessages(run_weekday_stage("activity", sr, er, b, make_onset_grid(cfg)[role == "primary" & stream == "activity"], cfg, ck, dk)))
  }
  st <- run(cfg); expect_true(all(st %in% c("computed", "resumed")) && length(st) == 17L)       # 16 fits + the descriptions
  wr <- collect_weekday_results(ck, cfg, list(activity = dk))
  expect_equal(nrow(wr$fits), 16L); expect_equal(nrow(wr$failures), 0L)
  for (i in seq_len(nrow(wr_test_pairs()))) {
    p <- wr_test_pairs()[i]; big <- wr$fits[fit_id == p$big]; small <- wr$fits[fit_id == p$small]
    expect_equal(big$n_rows, small$n_rows)                                                   # the same rows
    if (p$id != "H3") expect_equal(big$re_structure, small$re_structure)                     # the same random-effects structure
  }
  expect_true(all(wr$fits[fit_id %in% c("S0", "S1", "S2", "S3", "S4", "S1_holidays")]$re_structure == wr$fits[fit_id == "S1"]$re_structure))
  expect_lt(wr$fits[fit_id == "S1_holidays", n_rows], wr$fits[fit_id == "S1", n_rows]); expect_lt(wr$fits[fit_id == "S1_before_first_gap", n_rows], wr$fits[fit_id == "S1", n_rows])
  expect_equal(wr$fits[fit_id == "S1_activity", n_rows], wr$fits[fit_id == "S1_activity_same_rows", n_rows])
  fl <- wr$flow[step == "fit S1_holidays"]                                                    # rows left out of a fit, with the reason
  expect_equal(fl$n_left_out, wr$flow[step == "fit S0", n_left_out] + wr$fits[fit_id == "S0", n_rows] - wr$fits[fit_id == "S1_holidays", n_rows])
  expect_match(fl$reason, "holiday")
  expect_equal(unname(unlist(lapply(c(0, 3, 6), function(k) wr$profile[fit_id == "S1" & weekday == k, ratio_to_average_day] > 0))), rep(TRUE, 3))
  expect_true(all(run(cfg) == "resumed"))                                                    # nothing refitted
  cfg2 <- cfg; cfg2$weekday_rhythm$holiday_window <- 2L
  key <- function(c_) ckpt_key(c_, "weekday", dk, list(id = "S1", weekday_rhythm = c_$weekday_rhythm))
  expect_match(ckpt_mismatch(key(cfg), key(cfg2)), "spec")                                   # the settings are part of every key: a change makes the unit stale
  expect_equal(ckpt_mismatch(key(cfg), key(cfg)), "")
  figs <- wr_figures(wr, "activity", cfg)
  expect_true(all(c("weekday_rhythm_profile", "weekday_rhythm_calendar_vs_cycle_activity", "weekday_rhythm_personal_activity", "weekday_rhythm_cohort", "weekday_rhythm_cohort_level",
                    "weekday_rhythm_sensitivity_start", "weekday_rhythm_sensitivity_end", "weekday_rhythm_gap_length", "weekday_rhythm_dsr_counts_activity") %in% names(figs)))
  expect_true(all(vapply(figs, function(p) inherits(p, c("ggplot", "gtable", "grob")), logical(1))))
  expect_equal(attr(figs$weekday_rhythm_calendar_vs_cycle_activity, "size"), c(11, 11))
  d <- tempfile(); dir.create(d); write_weekday_summary(wr, cfg, d)
  expect_true(file.exists(file.path(d, "results_summary.md")))
})
