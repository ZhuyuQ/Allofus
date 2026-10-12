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
  p <- plot_cohort_forest(spec, data.table(stream = "activity", chisq = 2.9, df = 3L, p = 0.4), cfg); ok(p)
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

test_that("every 'missing by ...' line figure uses the same configuration: a line and a shaded band per cohort, one panel per stream", {
  cfg <- make_config()
  rt <- function(x) CJ(stream = c("activity", "sleep"), cohort = cfg$cohorts, x = x)[, `:=`(estimate = 0.05, lo = 0.03, hi = 0.07, n_persons = 50L)]
  mk <- function(name, x) setnames(rt(x), "x", name)[]
  plots <- list(week = plot_missing_over_time(mk("week", 1:10), "week", cfg, "Week", "t"),
                year = plot_missing_over_time(mk("year", 2020:2023), "year", cfg, "Year", "t"),
                dow = plot_missing_by_weekday(mk("dow", DOW_LABELS), cfg),
                age = plot_missing_by_age(mk("age", 40:70), cfg))
  for (p in plots) {
    geoms <- vapply(p$layers, function(l) class(l$geom)[1L], "")
    expect_true("GeomRibbon" %in% geoms); expect_true("GeomLine" %in% geoms)
    expect_s3_class(p$facet, "FacetWrap"); expect_silent(suppressWarnings(ggplot2::ggplot_build(p)))
  }
  expect_equal(plots$age$labels$x, "Age (years)")
  for (p in c(plots, list(ecdf = plot_gap_ecdf(data.table(cohort = factor(rep(cfg$cohorts, each = 5), levels = cfg$cohorts), stream = "activity", gap_max = rep(0:4, 4L)), cfg),
                          risk = plot_risk_curve(data.table(stream = c("activity", "sleep"), within_z = 0, risk = 0.02), cfg))))
    expect_equal(p$facet$params$ncol, 1L)                                     # panels are listed in rows, not columns
  empty <- plot_missing_by_age(mk("age", 40:41)[0L], cfg)                      # nothing reaches 30 people: a message, not an error
  expect_s3_class(empty, "ggplot"); expect_silent(suppressWarnings(ggplot2::ggplot_build(empty)))
})

test_that("stratified % missing figure draws the dodged violins with the within-cohort test labels", {
  cfg <- make_config(); set.seed(2)
  long <- CJ(person_id = 1:240, stream = c("activity", "sleep"))[, `:=`(cohort = factor(rep(cfg$cohorts, each = 60)[person_id %% 240 + 1L], levels = cfg$cohorts),
        sex = factor(sample(c("Female", "Male"), .N, TRUE), levels = c("Female", "Male", "Other")), value = rexp(.N, 0.2), metric = stream)]
  tests <- missing_pct_stratum_tests(long, "sex")
  p <- plot_missing_pct_by_stratum(long, "sex", tests, cfg, "t")
  expect_s3_class(p, "ggplot"); expect_silent(suppressWarnings(ggplot2::ggplot_build(p)))
  expect_equal(p$labels$fill, "Sex"); expect_equal(stratum_label("race_plot"), "Race")
  expect_error(plot_missing_pct_by_stratum(long, "nope", tests, cfg, "t"), "unknown stratum")
})

test_that("estimates with intervals are drawn as box plots (line = estimate, box = 50%, whiskers = 95%); dot styles remain available", {
  z <- stats::qnorm(0.975); q <- stats::qnorm(0.75)
  d <- data.table(est = 1, lo = 1 - z * 0.2, hi = 1 + z * 0.2)                       # a normal estimate: SE = 0.2
  add_interval_columns(d, "est", "lo", "hi"); expect_equal(c(d$lo_, d$q25_, d$est_, d$q75_, d$hi_), 1 + c(-z, -q, 0, q, z) * 0.2)
  dl <- data.table(or = 2, lo = exp(log(2) - z * 0.1), hi = exp(log(2) + z * 0.1))   # odds ratio: on the log scale
  add_interval_columns(dl, "or", "lo", "hi", log = TRUE); expect_equal(c(dl$q25_, dl$q75_), exp(log(2) + c(-q, q) * 0.1))
  dp <- data.table(p = 0.3, lo = 0.2, hi = 0.5); add_interval_columns(dp, "p", "lo", "hi", scale = 100)   # asymmetric interval stays asymmetric and nested
  expect_true(dp$q25_ > dp$lo_ && dp$q25_ < dp$est_ && dp$q75_ > dp$est_ && dp$q75_ < dp$hi_); expect_gt(dp$hi_ - dp$est_, dp$est_ - dp$lo_)
  ds <- data.table(est = 1, lo = 1 - z * 0.2, hi = 1 + z * 0.2); add_interval_columns(ds, "est", "lo", "hi", style = "se"); expect_equal(c(ds$lo_, ds$hi_), c(0.8, 1.2))
  d95 <- data.table(est = 1, lo = 0.6, hi = 1.4); add_interval_columns(d95, "est", "lo", "hi", style = "ci95"); expect_equal(c(d95$lo_, d95$hi_), c(0.6, 1.4))

  cfg <- make_config(); ft <- data.table(label = c("a", "b"), or = c(1.5, 1.2), or_low = c(1.2, 0.9), or_high = c(1.9, 1.6))
  geoms <- function(p) vapply(p$layers, function(l) class(l$geom)[1L], "")
  expect_true("GeomBoxplot" %in% geoms(plot_forest(ft, "t"))); expect_false("GeomPointrange" %in% geoms(plot_forest(ft, "t")))
  for (st in c("se", "ci95")) expect_true("GeomPointrange" %in% geoms(plot_forest(ft, "t", style = st)))
  expect_match(plot_forest(ft, "t")$labels$caption, "box: 50% interval"); expect_match(plot_forest(ft, "t", style = "se")$labels$caption, "1 SD")

  coefs <- data.table(stream = "activity", role = "primary", term = "w_act_steps", or = 1.6, or_low = 1.4, or_high = 1.8, n_events = 5L, n_persons = 3L)
  cs <- data.table(stream = "activity", cohort_stratum = "Healthy", or = 1.5, or_low = 1.2, or_high = 1.9)
  pred <- CJ(stream = "activity", adjusted = TRUE, within_z = c(-1, 0, 1), length_group = factor(c("1", "2"), ordered = TRUE))[, `:=`(prob = 0.5, prob_low = 0.4, prob_high = 0.6)]
  tr <- data.table(stream = "activity", cohort = cfg$cohorts, p_onset = 0.02, p_onset_low = 0.01, p_onset_high = 0.03,
                   p_continue = 0.5, p_continue_low = 0.4, p_continue_high = 0.6)
  inter <- data.table(stream = "activity", term = "w_act_steps", chisq = 3.2, df = 3L, p = 0.36, status = "ok")
  for (p in list(plot_onset_forest(coefs, cfg), plot_cohort_forest(cs, inter, cfg), plot_ordinal_probs(pred, cfg), plot_transition(tr, cfg))) {
    expect_true("GeomBoxplot" %in% geoms(p)); expect_silent(suppressWarnings(ggplot2::ggplot_build(p)))
  }
  cfg_se <- make_config(list(plot = list(interval_style = "se")))
  expect_true("GeomPointrange" %in% geoms(plot_ordinal_probs(pred, cfg_se)))
  txt <- vapply(plot_cohort_forest(cs, inter, cfg)$layers, function(l) if (is.null(l$data) || inherits(l$data, "waiver")) "" else paste(unlist(l$data), collapse = " "), "")
  expect_true(any(grepl("LRT chisq\\(3\\)=3.2, p=0.36", txt)))                                    # the cohort-difference test is labelled with its statistic and df
  expect_match(caption_for("cohort_forest"), "central 50% interval"); expect_match(caption_for("cohort_forest"), "likelihood-ratio test")
  expect_false(grepl("central 50%", caption_for("gap_ecdf")))
})

test_that("violin figures with more than 12 violins in total are placed in 2 or 3 rows with a size that follows the layout", {
  expect_null(violin_layout(4, 3, 3L)$nrow); expect_null(violin_layout(4, 3, 3L)$size)                  # 12 violins: unchanged
  expect_null(violin_layout(30, 1, 3L)$size)                                                              # a single panel: unchanged
  l2 <- violin_layout(4, 5, 2L); expect_equal(c(l2$nrow, l2$ncol), c(2L, 3L))                            # 20 violins: 2 rows
  l3 <- violin_layout(12, 3, 3L); expect_equal(c(l3$nrow, l3$ncol), c(3L, 1L)); expect_gte(l3$size[1L], 9)  # 36 violins, 3 panels: 3 rows of one wide panel
  l7 <- violin_layout(4, 7, 2L); expect_equal(c(l7$nrow, l7$ncol), c(3L, 3L))                            # 28 violins, 7 panels: 3 x 3
  expect_gt(l3$size[2L], l2$size[2L]); expect_lte(max(violin_layout(12, 8, 2L)$size[1L], 0), 24)
  expect_equal(violin_layout(4, 2, 2L)$nrow, NULL)                                                       # 8 violins

  cfg <- make_config(); set.seed(1); n <- 600
  long <- data.table(person_id = 1:n, stream = rep(c("activity", "sleep", "hr"), each = 200), cohort = factor(rep(cfg$cohorts, length.out = n), levels = cfg$cohorts),
                     age_group = factor(sample(c("18-44", "45-64", ">=65"), n, TRUE), levels = c("18-44", "45-64", ">=65")), value = rexp(n, 0.2), metric = rep(c("activity", "sleep", "hr"), each = 200))
  tests <- missing_pct_stratum_tests(long, "age_group")
  p <- plot_missing_pct_by_stratum(long, "age_group", tests, cfg, "t")          # 12 violins x 3 panels = 36 > 24
  expect_equal(p$facet$params$nrow, 3L); expect_equal(p$facet$params$ncol, 1L); expect_equal(length(attr(p, "size")), 2L)
  long$metric <- long$stream; long[, value := value]
  ml <- data.table(cohort = factor(rep(cfg$cohorts, 60), levels = cfg$cohorts), metric = rep(c("a", "b", "c", "d", "e"), each = 48), value = rexp(240))
  expect_equal(plot_metric_violin(ml, cfg, "t")$facet$params$nrow, 2L)             # 4 violins x 5 panels = 20: 2 rows
  small <- ml[metric %in% c("a", "b", "c")]; expect_null(attr(plot_metric_violin(small, cfg, "t"), "size"))   # 12 in total: as before
})
