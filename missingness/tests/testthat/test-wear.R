# Wear time of the last recorded day: the description (wear.R) and the row filter of the analyses wear_ge<h>h (onset_model.R).

test_that("wear_days marks the last day before a gap and joins the heart-rate minutes; a day with no heart-rate record is NA", {
  d0 <- as.Date("2024-03-01")
  days <- data.table(person_id = c(1L, 1L, 1L, 1L, 1L, 2L, 2L), date = d0 + c(0, 1, 2, 5, 6, 0, 1))
  wear <- data.table(person_id = c(1L, 1L, 1L, 2L), date = d0 + c(0, 1, 5, 0), wear_min = c(1440, 1000, 700, 1300))      # no heart-rate record on days 2 (person 1) and 1 (person 2)
  persons <- data.table(person_id = 1:2, cohort = c("Healthy", "Diabetes"))
  x <- wear_days(days, wear, persons, "activity")
  expect_equal(nrow(x), 5L)                                              # each person's last record has no next day
  expect_equal(x[order(person_id, date)]$last_day, c(FALSE, FALSE, TRUE, FALSE, FALSE))      # only day 2 of person 1 is followed by a gap (days 3-4 missing)
  expect_equal(x[order(person_id, date)]$wear_min, c(1440, 1000, NA, 700, 1300))
  expect_equal(x[person_id == 2L, unique(cohort)], "Diabetes")
})

test_that("the description reproduces hand-computed shares, the contrast, and the distribution is a proper ECDF", {
  s <- small_study(-0.5, seed = 3); cfg <- s$cfg
  x <- wear_days(s$inp$day_tables$activity, sim_wear(s$inp$day_tables), s$inp$persons, "activity")
  r <- wear_describe_stream(x, cfg)
  expect_true(all(abs(r$distribution[, .(m = max(cum_share)), by = .(cohort, group)]$m - 1) < 1e-12))     # every ECDF ends at 1
  expect_true(r$distribution[, all(diff(cum_share) >= 0), by = .(cohort, group)]$V1 |> all())
  expect_true(all(r$distribution$wear_min >= 0 & r$distribution$wear_min <= cfg$wear$max_min))
  for (co in cfg$cohorts) for (h in cfg$wear$hours) {
    k <- x[cohort == co & !is.na(wear_min)]
    sh <- function(g) 100 * mean(k[last_day == g]$wear_min < h * 60)
    sm <- r$summary[cohort == co & threshold_h == h]
    expect_equal(sm[group == "last", pct_below], sh(TRUE)); expect_equal(sm[group == "ordinary", pct_below], sh(FALSE))
    ct <- r$contrast[cohort == co & quantity == sprintf("%% of days under %g h", h)]
    expect_equal(ct$difference, sh(TRUE) - sh(FALSE)); expect_equal(ct$last, sh(TRUE))
    expect_lt(ct$low, ct$difference); expect_gt(ct$high, ct$difference)
  }
  # days with no heart-rate record are counted, not drawn
  s1 <- r$summary[cohort == cfg$cohorts[1] & threshold_h == cfg$wear$hours[1]]
  expect_equal(s1$n_no_hr_record, c(x[cohort == cfg$cohorts[1] & is.na(wear_min), sum(last_day)], x[cohort == cfg$cohorts[1] & is.na(wear_min), sum(!last_day)])[match(s1$group, c("last", "ordinary"))])
  expect_equal(s1$n_days, s1$n_known + s1$n_no_hr_record)
  expect_true(all(r$contrast[, .N, by = .(cohort, quantity)]$N == 1L))
  expect_true(all(r$contrast$p_holm >= r$contrast$p - 1e-12))
  expect_s3_class(plot_wear_last_day(r$distribution, cfg), "ggplot")
})

test_that("the hr stream carries its own heart-rate minutes; the other streams read day_tables$wear", {
  s <- small_study(-0.5, seed = 3); cfg <- s$cfg
  bh <- suppressMessages(prepare_stream_bundle("hr", s$inp$day_tables, s$inp$persons, cfg))
  hr <- s$inp$day_tables$hr
  expect_equal(bh$rs$wear_min, hr[bh$rs[, .(person_id, date)], on = c("person_id", "date"), total_hr_min])
  ba <- suppressMessages(prepare_stream_bundle("activity", s$inp$day_tables, s$inp$persons, cfg))        # no wear table given
  expect_false("wear_min" %in% names(ba$rs))
  expect_error(build_onset_data(ba, primary_spec(cfg), cfg, list(day = "day_before_gap", wear_min = 960)), "no wear_min")
  w <- sim_wear(s$inp$day_tables)
  miss <- s$b$rs[!w, on = c("person_id", "date")]                                         # risk-set days with no heart-rate record
  expect_gt(nrow(miss), 0L); expect_true(all(is.na(miss$wear_min)))
  expect_false(anyNA(s$b$rs[w, on = c("person_id", "date"), nomatch = 0L]$wear_min))
})

test_that("wear_min keeps exactly the rows whose own day d has enough heart-rate minutes, with the predictors of day d and the usual unchanged", {
  s <- small_study(-0.5, seed = 3); cfg <- s$cfg; b <- s$b; spec <- primary_spec(cfg)
  rs <- b$rs; all_preds <- unique(c(cfg$own_predictors[[b$stream]], cfg$activity_predictors))
  cols <- c(paste0("w_", all_preds), paste0("b_", all_preds)); adj <- c("age10", "sex2", "cohort")
  cc <- complete.cases(rs[, c(cols, adj), with = FALSE])
  dd <- build_onset_data(b, spec, cfg, list(day = "day_before_gap"))
  for (h in c(16, 20)) {
    w <- build_onset_data(b, spec, cfg, list(day = "day_before_gap", wear_min = h * 60, keep_wear = TRUE))
    ex <- rs[cc & (wear_min >= h * 60) %in% TRUE]
    expect_equal(nrow(w), nrow(ex)); expect_equal(w$w_act_steps, ex$w_act_steps)       # day d's predictors, the usual of the primary model
    expect_true(all(w$wear_min >= h * 60)); expect_false(anyNA(w$wear_min))               # no heart-rate record = not enough wear
    expect_lt(nrow(w), nrow(dd))
  }
  expect_lt(nrow(build_onset_data(b, spec, cfg, list(day = "day_before_gap", wear_min = 1200))), nrow(build_onset_data(b, spec, cfg, list(day = "day_before_gap", wear_min = 960))))
  expect_false("wear_min" %in% names(dd))                                                 # not kept unless asked
})

test_that("wear_filter_flow adds up: unknown + classes = all rows of day d; kept rows = the rows of each analysis", {
  s <- small_study(-0.5, seed = 3); cfg <- s$cfg; b <- s$b; spec <- primary_spec(cfg)
  f <- wear_filter_flow(b, spec, cfg)
  all_d <- f[rows == "Day d: every row with complete predictors"]
  unk <- f[grepl("no heart-rate record", rows)]
  cls <- f[grepl("^  wear on d ", rows)]
  expect_equal(nrow(cls), length(cfg$wear$hours) + 1L)
  expect_equal(unk$n_rows + sum(cls$n_rows), all_d$n_rows); expect_equal(unk$n_events + sum(cls$n_events), all_d$n_events)
  expect_equal(f[rows == "Primary model rows (predictors of the last full day d - 1)", n_rows], nrow(build_onset_data(b, spec, cfg)))
  for (h in cfg$wear$hours) {
    k <- f[grepl(sprintf("^Kept by wear_ge%gh", h), rows)]
    w <- build_onset_data(b, spec, cfg, list(day = "day_before_gap", wear_min = h * 60))
    expect_equal(c(k$n_rows, k$n_events), c(nrow(w), sum(w$y_ge1)))
  }
  expect_true(all(f$event_pct[f$n_rows > 0] >= 0))
})

test_that("the wear stage writes the distribution, summary, contrast tables and the figure, filed under data quality", {
  cfg <- cfg_test()
  inp <- suppressMessages(simulate_inputs(cfg, n_per_cohort = c(Insomnia = 40, Diabetes = 40, Healthy = 80, Prostate = 20), seed = 5))
  out <- tempfile(); r <- NULL
  invisible(suppressWarnings(suppressMessages(capture.output(r <- run_missingness_pipeline(inp$day_tables, inp$persons, cfg, out, stages = "wear")))))
  expect_true(all(file.exists(file.path(out, "tables", c("wear_last_day_distribution.tsv", "wear_last_day.tsv", "wear_last_day_contrast.tsv")))))
  expect_true(file.exists(file.path(out, "figures", "wear_last_day.png")))
  expect_setequal(unique(r$wear_last_day$stream), c("activity", "sleep", "hr"))
  expect_setequal(unique(r$wear_last_day_contrast$quantity), c("% of days under 16 h", "% of days under 20 h", "mean wear (hours)"))
  ix <- fread(file.path(out, "figures", "figure_index.tsv")); expect_equal(ix[figure == "wear_last_day.png", step], 3L)
  expect_false(is.na(caption_for("wear_last_day")))
})
