test_that("every plot function builds a ggplot that renders", {
  cfg <- cfg_test()
  inp <- simulate_inputs(cfg, n_per_cohort = c(Insomnia = 40, Diabetes = 40, Healthy = 80, Prostate = 30), seed = 12)
  dts <- inp$day_tables; per <- inp$persons
  gaps <- lapply(dts, get_gap_runs)
  pg <- rbindlist(lapply(names(dts), function(s) with_cohort(summarise_person_gaps(dts[[s]], gaps[[s]]), per, s)))
  gaps_all <- rbindlist(lapply(names(dts), function(s) with_cohort(gaps[[s]], per, s)))
  # data warnings (e.g. a violin with <2 points in a tiny stratum) are expected on small data; errors are not
  ok <- function(p) { expect_s3_class(p, "ggplot"); expect_error(suppressWarnings(ggplot2::ggplot_build(p)), NA) }

  pms <- person_metric_stats(dts$activity, descriptive_metrics()$activity)
  long <- person_means_long(pms, per, "activity")
  suppressWarnings({
    ok(plot_metric_violin(long, cfg, "t"))
    for (st in c("sex", "age_group", "race_plot", "cci_group")) ok(plot_metric_by_stratum(long, st, cfg, "t"))
    ok(plot_single_metric_violin(long[metric == "steps"], "value", cfg, "t", "y", ref_line = 10000))
    ok(plot_gap_ecdf(pg, cfg)); ok(plot_gap_length_distribution(gaps_all, cfg))
    counts <- rbindlist(lapply(unique(pg$stream), function(s) transition_person_counts(pg[stream == s])))
    ok(plot_transition(transition_table(counts, cfg), cfg))
    wk <- cbind(stream = "activity", study_week_counts(pg[stream == "activity"], gaps$activity))
    ok(plot_missing_over_time(rate_table(wk, c("stream", "cohort", "week"), cfg), "week", cfg, "Week", "t"))
    dw <- cbind(stream = "activity", dow_counts(pg[stream == "activity"], gaps$activity))
    ok(plot_missing_by_weekday(rate_table(dw, c("stream", "cohort", "dow"), cfg, 1L), cfg))
    ok(plot_missing_heatmap(dts$activity, cfg, n_people = 30L))
    ok(plot_forest(data.table(label = c("a", "b"), or = c(1.5, 1.2), or_low = c(1.2, 0.9), or_high = c(1.9, 1.6)), "t"))
    ok(plot_risk_curve(data.table(stream = "activity", within_z = c(-1, 0, 1), risk = c(0.03, 0.02, 0.015)), cfg))
  })
  f <- tempfile(fileext = ".png"); save_plot(plot_gap_ecdf(pg, cfg), f, 6, 4, cfg); expect_true(file.exists(f))
})

test_that("model-result figures draw from saved result tables and are labelled honestly", {
  cfg <- make_config()
  ok <- function(p) { expect_s3_class(p, "ggplot"); expect_silent(suppressWarnings(ggplot2::ggplot_build(p))) }
  coefs <- data.table(stream = "activity", role = "primary", term = c("w_act_steps", "w_act_mvpa", "b_act_steps"),
                      or = c(1.6, 1.0, 0.9), or_low = c(1.4, 0.9, 0.8), or_high = c(1.8, 1.1, 1.0),
                      n_events = 500L, n_persons = 300L)
  ok(plot_onset_forest(coefs, cfg))
  spec <- data.table(stream = "activity", cohort_stratum = c("Insomnia", "Healthy"), or = c(1.5, 1.7), or_low = c(1.2, 1.4), or_high = c(1.9, 2.1))
  p <- plot_cohort_forest(spec, data.table(stream = "activity", p = 0.4), cfg); ok(p)
  expect_match(p$labels$subtitle, "Not estimable.*Diabetes.*Prostate")               # missing cohorts are named, not dropped silently
  risk <- data.table(stream = "activity", within_z = c(-1, 0, 1), risk = c(0.03, 0.02, 0.015),
                     risk_low = c(0.025, 0.017, 0.012), risk_high = c(0.036, 0.024, 0.019))
  ok(plot_risk_curve(risk, cfg))
  curve <- CJ(stream = "activity", model_id = c("B_back2_unadj", "B_back2_adj"), gap_length = 1:14)[, `:=`(est = -0.01 * gap_length, lo = -0.01 * gap_length - 0.05, hi = -0.01 * gap_length + 0.05)]
  obs <- data.table(stream = "activity", model_id = "B_back2_unadj", x = c(1, 2, 5, 10), mean = c(0, -0.02, -0.04, -0.1),
                    lo = c(-0.05, -0.08, -0.12, -0.3), hi = c(0.05, 0.04, 0.04, 0.1), n_gaps = c(300L, 150L, 60L, 10L))
  ok(plot_gap_dose_response(curve, obs, cfg))
  pred <- CJ(stream = "activity", adjusted = TRUE, within_z = c(-1, 0, 1), length_group = factor(c("1", "2", "3+"), ordered = TRUE))[
    , `:=`(prob = 1 / 3, prob_low = 0.25, prob_high = 0.4)]
  ok(plot_ordinal_probs(pred, cfg))
})

test_that("every figure the pipeline can make has a caption", {
  nm <- c("missing_pct_by_cohort", "longest_gap_by_cohort", "gap_ecdf", "gap_length_distribution", "heatmap_sleep", "transition",
          "missing_by_study_week", "missing_by_year", "missing_by_weekday", "metrics_hr_by_sex", "onset_primary_forest",
          "forest_gap_definition", "risk_curve", "cohort_forest", "forest_sensitivity", "gap_dose_response", "gap_length_probs")
  expect_false(anyNA(vapply(nm, caption_for, character(1))))
  expect_true(is.na(caption_for("not_a_figure")))
  d <- tempfile(); write_figure_captions(c("gap_ecdf", "risk_curve"), d)
  expect_match(paste(readLines(file.path(d, "figures", "figure_captions.md")), collapse = " "), "Population-average")
})
