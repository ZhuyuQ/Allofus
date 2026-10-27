# Observed steps on the days before a gap (onset_pre_gap_trajectory): the onset model's counterpart of the observed post-gap activity of model B.

test_that("the difference of two cluster-sampled means agrees with a regression with a cluster-robust variance", {
  set.seed(8); id1 <- sample(1:300, 400, TRUE); id0 <- sample(1:300, 5000, TRUE)
  u <- rnorm(300); x1 <- rnorm(400, 0.4) + u[id1]; x0 <- rnorm(5000) + u[id0]
  r <- cluster_diff_ci(x1, id1, x0, id0)
  d <- data.table(x = c(x1, x0), g = rep(1:0, c(400, 5000)), id = c(id1, id0))
  fit <- stats::lm(x ~ g, d); V <- sandwich::vcovCL(fit, cluster = d$id, type = "HC1")
  expect_equal(r$difference, unname(coef(fit)[["g"]]), tolerance = 1e-10)
  expect_equal(r$se, sqrt(V[2L, 2L]), tolerance = 0.05)                                           # same sandwich, up to the small-sample constant
  expect_true(r$low < r$difference && r$difference < r$high)
  expect_true(is.na(cluster_diff_ci(1, 1, x0, id0)$p))
})

test_that("trajectory means equal an independent calculation (previous CALENDAR day), by gap-length group, with controls and Holm p", {
  s <- small_study(-0.5); cfg <- s$cfg; b <- s$b
  tr <- onset_pre_gap_trajectory(b, cfg, days_before = 0:3)
  expect_setequal(as.character(unique(tr$group[tr$group == "all gaps"])), "all gaps"); expect_true(all(tr$days_before %in% 0:3))
  rs <- copy(b$rs)[order(person_id, date)]
  rs[, `:=`(prev_date = shift(date), prev_z = shift(w_act_steps)), by = person_id]
  ev <- rs[y_ge1 == 1L & !is.na(prev_date) & as.integer(date - prev_date) == 1L & !is.na(prev_z)]
  ct <- rs[y_ge1 == 0L & !is.na(prev_date) & as.integer(date - prev_date) == 1L & !is.na(prev_z)]
  r1 <- tr[group == "all gaps" & days_before == 1L]
  expect_equal(r1$mean, mean(ev$prev_z)); expect_equal(r1$n_gaps, nrow(ev)); expect_equal(r1$mean_control, mean(ct$prev_z)); expect_equal(r1$n_control_rows, nrow(ct))
  r0 <- tr[group == "all gaps" & days_before == 0L]
  expect_equal(r0$mean, mean(rs[y_ge1 == 1L & !is.na(w_act_steps), w_act_steps])); expect_equal(r0$difference, r0$mean - r0$mean_control)
  expect_true(all(tr$lo <= tr$mean & tr$mean <= tr$hi, na.rm = TRUE) && all(tr$p >= 0 & tr$p <= 1, na.rm = TRUE))
  # a gap of length L is in exactly one length group; the groups add up to all gaps (day 0 has a z-score for every event)
  g0 <- tr[days_before == 0L & group != "all gaps", sum(n_gaps)]; expect_equal(g0, r0$n_gaps)
  # the number of gaps can only fall with the distance (a day without a record has no z-score)
  expect_true(all(diff(tr[group == "all gaps"][order(days_before), n_gaps]) <= 0))
})

test_that("a partial last day shows as a drop on day 0 only; the plot and the stage unit build", {
  cfg <- cfg_test()
  inp <- suppressMessages(simulate_inputs(cfg, n_per_cohort = c(Insomnia = 100, Diabetes = 100, Healthy = 200, Prostate = 50), beta = 0, seed = 31))
  act <- copy(inp$day_tables$activity); setorder(act, person_id, date)
  act[, next_date := shift(date, type = "lead"), by = person_id]
  last_day <- !is.na(act$next_date) & as.integer(act$next_date - act$date) > 1L
  set.seed(5); act[last_day, steps := steps * stats::runif(.N, 0.15, 0.6)]; act[, next_date := NULL]
  tabs <- inp$day_tables; tabs$activity <- act
  b <- suppressMessages(prepare_stream_bundle("activity", tabs, inp$persons, cfg))
  tr <- cbind(stream = "activity", onset_pre_gap_trajectory(b, cfg))
  tr[, p_holm := stats::p.adjust(p, "holm"), by = .(stream, group)]
  a <- tr[group == "all gaps"]
  expect_lt(a[days_before == 0L, difference], -1); expect_lt(a[days_before == 0L, p_holm], 1e-6)
  expect_lt(max(abs(a[days_before >= 1L, difference])), 0.25)                                       # the earlier days are as usual
  p <- plot_pre_gap_observed(tr, "activity", cfg); expect_silent(suppressWarnings(ggplot2::ggplot_build(p)))
  expect_match(p$labels$caption, "mean steps z-score"); expect_false("series" %in% names(p$data)); expect_equal(length(unique(p$layers[[2]]$aes_params$colour)), 1L)    # only the days followed by a gap are drawn
  expect_null(plot_pre_gap_observed(tr, "hr", cfg)); expect_null(plot_pre_gap_observed(tr[0L], "activity", cfg))
  expect_match(levels(p$data$panel)[length(levels(p$data$panel))], "all gaps: ")
  expect_equal(figure_analysis("pre_gap_activity_observed_hr")$step, 4L); expect_match(caption_for("pre_gap_activity_observed_hr"), "BEFORE a gap")
  # the stage computes it as one unit, resumes it, and collects it with Holm p
  dir <- file.path(tempdir(), "pg_ck"); unlink(dir, recursive = TRUE)
  dk <- data_key(as_day_store(tabs, inp$persons), "activity"); grid <- make_onset_grid(cfg)[stream == "activity" & role == "primary"]
  st <- suppressMessages(run_onset_stage("activity", lazy(function() b), grid, cfg, dir, dk))
  expect_true("extra:pre_gap_trajectory" %in% names(st))
  got <- collect_onset_results(dir, cfg, list(activity = dk))
  expect_true(nz(got$onset_pre_gap_trajectory) && "p_holm" %in% names(got$onset_pre_gap_trajectory))
  again <- suppressMessages(run_onset_stage("activity", lazy(function() stop("must not be built")), grid, cfg, dir, dk))
  expect_equal(unname(again["extra:pre_gap_trajectory"]), "resumed")
})

test_that("the usual of the activity predictors is taken over the days observed in the modelled stream (no offset from days without a sleep record)", {
  cfg <- cfg_test()
  inp <- suppressMessages(simulate_inputs(cfg, n_per_cohort = c(Insomnia = 60, Diabetes = 60, Healthy = 120, Prostate = 30), seed = 17))
  tabs <- inp$day_tables
  # days WITHOUT a sleep record get far fewer steps (partial wear): with all activity days as the reference the sleep days would sit well above the usual
  sl <- tabs$sleep[, .(person_id, date, has_sleep = TRUE)]
  act <- merge(tabs$activity, sl, by = c("person_id", "date"), all.x = TRUE); act[is.na(has_sleep), steps := steps * 0.25]; act[, has_sleep := NULL]
  setkey(act, person_id, date); tabs$activity <- act
  b <- suppressMessages(prepare_stream_bundle("sleep", tabs, inp$persons, cfg))
  expect_lt(abs(mean(b$rs$w_act_steps, na.rm = TRUE)), 0.1)                                         # no offset: the usual is built from the same kind of day
  # independent check of one person-day: usual = mean of log1p(steps) over the person's SLEEP-observed days strictly before the day
  r <- b$rs[!is.na(w_act_steps)][1000L]
  pr <- merge(tabs$activity, tabs$sleep[, .(person_id, date)], by = c("person_id", "date"))[person_id == r$person_id & date < r$date]
  x <- log1p(tabs$activity[person_id == r$person_id & date == r$date, steps])
  sd_w <- b$scales[predictor == "act_steps", sd_within]
  expect_equal(r$w_act_steps, (x - mean(log1p(pr$steps))) / sd_w, tolerance = 1e-9)
  # the activity stream itself is unchanged: its usual is built from its own days
  ba <- suppressMessages(prepare_stream_bundle("activity", inp$day_tables, inp$persons, cfg))
  expect_lt(abs(mean(ba$rs$w_act_steps, na.rm = TRUE)), 0.3)
  # the gap table: the day BEFORE the gap is standardized on the stream's days (as in the onset rows); the days AFTER the gap are any activity days
  # with the usual of all activity days, so no day is dropped because the stream has no record on it (no selection on what happens after the gap)
  g <- build_gap_level_data(tabs$sleep, tabs$activity, inp$persons, cfg)
  gr <- get_gap_runs(tabs$sleep); a_all <- predictor_table(tabs$activity, "act_steps", cfg, trailing = FALSE, keep_x = TRUE, keep_baseline = TRUE); setkey(a_all, person_id, date)
  n_expected <- sum(!is.na(a_all[.(gr$person_id, gr$next_obs_date + 1L), x_act_steps]) & !is.na(a_all[.(gr$person_id, gr$prev_obs_date - 1L), m_act_steps]))        # the usual through the last FULL day d - 1
  expect_equal(sum(!is.na(g$back2_dev)), n_expected)
  rsw <- b$rs[, .(person_id, date, w = w_act_steps)]; setkey(rsw, person_id, date)
  expect_equal(sort(g$w_pre_day_steps), sort(rsw[.(gr$person_id, gr$prev_obs_date - 1L), w]))     # the day before the gap = the last full day d - 1, on the stream's days
  expect_equal(sort(g$alt_w_pre_day_steps), sort(rsw[.(gr$person_id, gr$prev_obs_date), w]))       # sensitivity: the last recorded day d
  st <- as_day_store(inp$day_tables, inp$persons)
  expect_true("activity_reference" %in% names(data_key(st, "sleep"))); expect_false("activity_reference" %in% names(data_key(st, "activity")))
})
