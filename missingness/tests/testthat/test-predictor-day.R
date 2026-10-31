# The predictor day of the onset model (onset_model.R, sensitivity.R): the PRIMARY model reads the predictors from the last FULL day d - 1 (the calendar day
# before the last recorded day d); the day right before the gap (d) is a sensitivity analysis.

test_that("predictor-day data: the primary reads day d - 1 (needs a record there), day_before_gap reads d, same_rows/both use the rows with both days", {
  s <- small_study(-0.5); cfg <- s$cfg; b <- s$b; spec <- primary_spec(cfg)
  expect_equal(cfg$onset$predictor_day, "last_full_day")
  dp <- build_onset_data(b, spec, cfg)                                                                 # primary: d - 1
  dd <- build_onset_data(b, spec, cfg, list(day = "day_before_gap"))                                   # sensitivity: d, every row with complete predictors on d
  ds <- build_onset_data(b, spec, cfg, list(day = "day_before_gap", same_rows = TRUE)); dl <- build_onset_data(b, spec, cfg, list(day = "last_full_day", same_rows = TRUE))
  db <- build_onset_data(b, spec, cfg, list(day = "both"))
  expect_equal(c(nrow(ds), nrow(dl)), rep(nrow(db), 2L)); expect_lt(nrow(db), nrow(dd)); expect_lte(nrow(db), nrow(dp)); expect_lt(nrow(dp), nrow(dd))
  # independent expectation: the previous ROW of the risk set is the previous observed day; it counts only if it is exactly one calendar day earlier
  rs <- copy(b$rs)[order(person_id, date)]
  all_preds <- unique(c(cfg$own_predictors[[b$stream]], cfg$activity_predictors))
  cols <- c(paste0("w_", all_preds), paste0("b_", all_preds))
  rs[, `:=`(prev_date = shift(date), prev_w = shift(w_act_steps)), by = person_id]
  for (cl in cols) rs[, (paste0("prev_", cl)) := shift(get(cl)), by = person_id]
  has_prev <- !is.na(rs$prev_date) & as.integer(rs$date - rs$prev_date) == 1L
  adj <- c("age10", "sex2", "cohort")
  okp <- rs[has_prev & complete.cases(rs[, c(paste0("prev_", cols), adj), with = FALSE])]
  okb <- rs[has_prev & complete.cases(rs[, c(cols, paste0("prev_", cols), adj), with = FALSE])]
  expect_equal(nrow(dp), nrow(okp)); expect_equal(dp$w_act_steps, okp$prev_w)                        # primary: the value of day d - 1
  expect_equal(nrow(db), nrow(okb))
  expect_equal(dl$w_act_steps, okb$prev_w); expect_equal(ds$w_act_steps, okb$w_act_steps)               # same rows: d - 1 / d
  expect_equal(db$lag1_w_act_steps, okb$prev_w); expect_equal(db$w_act_steps, okb$w_act_steps)         # joint mode: both
  expect_equal(nrow(dd), sum(complete.cases(rs[, c(cols, adj), with = FALSE])))                         # day_before_gap: every row with complete predictors on d
  expect_true(all(ds$y_ge1 == dl$y_ge1))                                                                # the outcome is the same: a gap starts on d + 1
  expect_error(build_onset_data(b, spec, cfg, list(day = "nope")), "unknown predictor day")
  # the config switch: the old definition is available as cfg$onset$predictor_day = "day_before_gap"
  cfg2 <- cfg; cfg2$onset$predictor_day <- "day_before_gap"
  expect_equal(build_onset_data(b, spec, cfg2)$w_act_steps, dd$w_act_steps)
})

test_that("partial-wear exclusion and the run-in rows follow the predictor day", {
  s <- small_study(-0.5); cfg <- s$cfg; b <- s$b; spec <- primary_spec(cfg)
  b2 <- b; b2$rs <- copy(b$rs); set.seed(3); b2$rs[, pw_act := sample(c(TRUE, FALSE), .N, TRUE, prob = c(0.1, 0.9))]; b2$rs[, pw_own := pw_act]
  all_preds <- unique(c(cfg$own_predictors[[b$stream]], cfg$activity_predictors)); cols <- c(paste0("w_", all_preds), paste0("b_", all_preds))
  rs <- copy(b2$rs)[order(person_id, date)]
  rs[, `:=`(prev_date = shift(date), prev_pw = shift(pw_act), prev_study_day = shift(study_day)), by = person_id]
  for (cl in cols) rs[, (paste0("prev_", cl)) := shift(get(cl)), by = person_id]
  okp <- rs[!is.na(prev_date) & as.integer(date - prev_date) == 1L & complete.cases(rs[, c(paste0("prev_", cols), "age10", "sex2", "cohort"), with = FALSE])]
  okd <- rs[complete.cases(rs[, c(cols, "age10", "sex2", "cohort"), with = FALSE])]
  a  <- build_onset_data(b2, spec, cfg); a2 <- build_onset_data(b2, spec, cfg, list(drop_partial = TRUE))
  d  <- build_onset_data(b2, spec, cfg, list(day = "day_before_gap")); d2 <- build_onset_data(b2, spec, cfg, list(day = "day_before_gap", drop_partial = TRUE))
  expect_equal(nrow(a) - nrow(a2), sum(okp$prev_pw))                      # the primary model drops the rows whose last FULL day d - 1 is a partial day ...
  expect_equal(nrow(d) - nrow(d2), sum(okd$pw_act))                       # ... the day-before-the-gap analysis the rows whose day d is
  r <- build_onset_data(b, spec, cfg, list(runin_days = cfg$runin$end)); expect_lt(nrow(r), nrow(a)); expect_true(all(!is.na(r$runin_missing)))
  expect_equal(nrow(r), sum(!is.na(okp$runin_missing) & okp$prev_study_day > cfg$runin$end))   # the days the predictors are read from lie after the run-in
})

test_that("the contrast helper agrees with a cluster-robust regression (CR1 sandwich)", {
  set.seed(4); n <- 3000; id <- sample(1:400, n, TRUE); y <- rbinom(n, 1, 0.1); x <- rnorm(n) + 0.3 * y + rnorm(400)[id]
  r <- day_contrast(x, y, id)
  fit <- stats::lm(x ~ y); V <- sandwich::vcovCL(fit, cluster = id, type = "HC1", cadjust = TRUE)
  expect_equal(r$difference, unname(coef(fit)[["y"]]), tolerance = 1e-10); expect_equal(r$se, sqrt(V[2L, 2L]), tolerance = 1e-8)
  expect_equal(r$mean_gap, mean(x[y == 1L])); expect_null(day_contrast(x[y == 0L], y[y == 0L], id[y == 0L]))
})

test_that("a partial last day is detected: with NO true effect, the association sits on day d (sensitivity) and not on the last full day d - 1 (primary)", {
  cfg <- cfg_test()
  inp <- suppressMessages(simulate_inputs(cfg, n_per_cohort = c(Insomnia = 100, Diabetes = 100, Healthy = 200, Prostate = 50), beta = 0, seed = 31))
  act <- copy(inp$day_tables$activity); setorder(act, person_id, date)
  act[, next_date := shift(date, type = "lead"), by = person_id]
  last_day <- !is.na(act$next_date) & as.integer(act$next_date - act$date) > 1L                    # a gap follows this day
  set.seed(5); act[last_day, steps := steps * stats::runif(.N, 0.15, 0.6)]                       # the device was taken off during the day: fewer steps, for a technical reason
  act[, next_date := NULL]
  tabs <- inp$day_tables; tabs$activity <- act
  b <- suppressMessages(prepare_stream_bundle("activity", tabs, inp$persons, cfg))
  spec <- primary_spec(cfg)
  fit <- function(v) suppressWarnings(fit_onset_spec(b, spec, cfg, variant = v))
  main <- fit(list()); cur <- fit(list(day = "day_before_gap", same_rows = TRUE)); lag <- fit(list(day = "last_full_day", same_rows = TRUE)); both <- fit(list(day = "both"))
  dall <- fit(list(day = "day_before_gap"))
  g <- function(f, nm) { x <- f$coef; x[x$term == nm, ] }
  expect_lt(g(dall, "w_act_steps")$or, 0.85); expect_lt(g(cur, "w_act_steps")$p, 1e-6)           # the artifact creates a strong "effect" on day d (the sensitivity analysis) ...
  expect_gt(g(main, "w_act_steps")$or, 0.9); expect_gt(g(main, "w_act_steps")$p, 0.01)           # ... that is absent in the primary model (last full day)
  expect_gt(g(lag, "w_act_steps")$or, 0.9); expect_gt(g(lag, "w_act_steps")$p, 0.01)
  expect_equal(cur$info$n_rows, lag$info$n_rows); expect_equal(cur$info$n_rows, both$info$n_rows)
  i <- both$info; dlt <- i$b_cur - i$b_lag; se <- sqrt(i$v_cur + i$v_lag - 2 * i$cov_cur_lag)
  expect_lt(dlt, 0); expect_lt(2 * stats::pnorm(-abs(dlt / se)), 1e-4)                             # the difference of the two estimates of the joint model is significant
  expect_equal(i$b_cur, g(both, "w_act_steps")$estimate); expect_equal(i$b_lag, g(both, "lag1_w_act_steps")$estimate)
  # the observed pattern: the drop from d - 1 to d is larger before a gap
  ds <- day_describe(b, spec, cfg)
  ch <- ds$descriptive[quantity == "Change from d - 1 to d"]
  expect_lt(ch$difference, -0.5); expect_lt(ch$p, 1e-6); expect_equal(nrow(ds$pattern), 4L)
  expect_true(ds$descriptive[quantity == "Last full day (d - 1)", abs(difference)] < 0.3)
})

test_that("with a real effect and no partial days the primary model recovers it; day d is attenuated; check table, figures and the stage build", {
  s <- small_study(-0.7, seed = 12); cfg <- s$cfg; b <- s$b; spec <- primary_spec(cfg)
  dir <- file.path(tempdir(), "fd_ck"); unlink(dir, recursive = TRUE)
  dk <- data_key(as_day_store(s$inp$day_tables, s$inp$persons), "activity")
  st <- suppressMessages(run_sensitivity_stage("activity", lazy(function() b), NULL, spec, cfg, dir, dk))
  expect_true(all(c("day_before_gap", "day_before_gap_same_rows", "last_full_day_same_rows", "both_days", "day_description") %in% names(st)))
  got <- collect_sensitivity_results(dir, cfg, list(activity = dk))
  chk <- got$day_check
  expect_setequal(as.character(chk$model), DAY_MODEL_LEVELS)
  df <- chk[grepl("^Difference", as.character(model))]
  expect_equal(df$or, exp(df$estimate)); expect_true(df$p >= 0 && df$p <= 1); expect_equal(df$p_holm, df$p)                  # one stream
  expect_true(is.finite(df$rho_d_dminus1) && is.finite(df$lagged_over_current))
  est <- function(m) chk[as.character(model) == m, estimate]
  expect_lt(est(DAY_MODEL_LEVELS[1L]), 0)                                                                                     # the planted effect is recovered on d - 1
  expect_lt(abs(est(DAY_MODEL_LEVELS[3L])), abs(est(DAY_MODEL_LEVELS[2L])))                                                   # and is weaker on d (one day further from the planted day)
  expect_equal(nrow(got$day_descriptive), 3L); expect_equal(nrow(got$day_pattern), 4L)
  tested <- got$sensitivity[analysis %in% c("day_before_gap", "day_before_gap_same_rows", "last_full_day_same_rows", "both_days")]
  expect_equal(uniqueN(tested$analysis), 4L); expect_false(any(duplicated(tested[, .(stream, analysis)])))                   # one tested row per analysis
  p <- plot_day_check(chk, "activity", cfg, got$day_descriptive); expect_silent(suppressWarnings(ggplot2::ggplot_build(p)))
  expect_match(p$labels$caption, "b_d - b_\\(d-1\\)"); expect_match(p$labels$caption, "Mean steps z-score")
  p2 <- plot_day_pattern(got$day_pattern, got$day_descriptive, "activity", cfg); expect_silent(suppressWarnings(ggplot2::ggplot_build(p2)))
  expect_match(p2$labels$caption, "Change from d - 1 to d")
  expect_null(plot_day_check(chk, "sleep", cfg)); expect_null(plot_day_pattern(got$day_pattern, NULL, "sleep", cfg))
  expect_equal(term_label("lag1_w_act_steps"), "Within: Steps, last full day (d - 1)")
  expect_equal(figure_analysis("day_check_activity")$step, 6L); expect_equal(table_analysis("day_check")$step, 6L)
  expect_match(caption_for("day_check_activity"), "partial day"); expect_match(caption_for("day_pattern_hr"), "partial last day")
  # the stage resumes
  again <- suppressMessages(run_sensitivity_stage("activity", lazy(function() stop("must not be built")), NULL, spec, cfg, dir, dk))
  expect_true(all(again[c("day_before_gap", "day_before_gap_same_rows", "last_full_day_same_rows", "both_days", "day_description")] == "resumed"))
})
