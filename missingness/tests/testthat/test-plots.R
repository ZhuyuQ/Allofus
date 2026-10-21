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
  gb <- vapply(plot_gap_dose_response(curve, obs, cfg)$layers, function(l) class(l$geom)[1L], "")
  expect_true("GeomPointrange" %in% gb); expect_false(any(c("GeomCol", "GeomLine", "GeomRibbon") %in% gb))       # forest rows, not bars or lines
  dd <- ggplot2::ggplot_build(plot_gap_dose_response(curve, obs, cfg))$plot$data
  expect_equal(sort(as.integer(as.character(unique(dd$day)))), c(1L, 2L, 3L, 5L, 13L))                      # only these gap lengths
  expect_equal(nrow(dd), 5L * 2L); expect_setequal(as.character(dd$model), c("Unadjusted", "Adjusted"))   # each with an unadjusted and an adjusted estimate
  expect_true("GeomBoxplot" %in% vapply(plot_gap_dose_response(curve, obs, make_config(list(plot = list(interval_style = "box"))))$layers, function(l) class(l$geom)[1L], ""))
  pred <- CJ(stream = "activity", adjusted = TRUE, within_z = c(-1, 0, 1), length_group = factor(c("1", "2", "3+"), ordered = TRUE))[
    , `:=`(prob = 1 / 3, prob_low = 0.25, prob_high = 0.4)]
  ok(plot_ordinal_probs(pred, cfg))
})

test_that("every figure the pipeline can make has a caption", {
  nm <- c("missing_pct_by_cohort", "longest_gap_by_cohort", "gap_ecdf", "gap_length_distribution", "heatmap_sleep", "transition",
          "missing_by_study_week", "missing_by_year", "missing_by_weekday", "metrics_hr_by_sex", "onset_primary_forest",
          "forest_gap_definition", "risk_curve", "cohort_forest", "forest_sensitivity", "gap_dose_response", "gap_length_probs",
          "onset_forest_other_activity", "onset_forest_other_hr", "sensitivity_forest_sleep", "gapB_forest_activity", "gapC_forest_hr", "onset_deciles_sleep")
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

test_that("model figures are dot +/- 1 SD by default; box and 95%-interval styles remain available; transition keeps its own setting", {
  z <- stats::qnorm(0.975); q <- stats::qnorm(0.75)
  d <- data.table(est = 1, lo = 1 - z * 0.2, hi = 1 + z * 0.2)                       # a normal estimate: SE = 0.2
  add_interval_columns(d, "est", "lo", "hi"); expect_equal(c(d$lo_, d$q25_, d$est_, d$q75_, d$hi_), 1 + c(-z, -q, 0, q, z) * 0.2)
  dl <- data.table(or = 2, lo = exp(log(2) - z * 0.1), hi = exp(log(2) + z * 0.1))   # odds ratio: on the log scale
  add_interval_columns(dl, "or", "lo", "hi", log = TRUE); expect_equal(c(dl$q25_, dl$q75_), exp(log(2) + c(-q, q) * 0.1))
  dp <- data.table(p = 0.3, lo = 0.2, hi = 0.5); add_interval_columns(dp, "p", "lo", "hi", scale = 100)   # asymmetric interval stays asymmetric and nested
  expect_true(dp$q25_ > dp$lo_ && dp$q25_ < dp$est_ && dp$q75_ > dp$est_ && dp$q75_ < dp$hi_); expect_gt(dp$hi_ - dp$est_, dp$est_ - dp$lo_)
  ds <- data.table(est = 1, lo = 1 - z * 0.2, hi = 1 + z * 0.2); add_interval_columns(ds, "est", "lo", "hi", style = "se"); expect_equal(c(ds$lo_, ds$hi_), c(0.8, 1.2))
  dls <- data.table(or = 2, lo = exp(log(2) - z * 0.1), hi = exp(log(2) + z * 0.1)); add_interval_columns(dls, "or", "lo", "hi", log = TRUE, style = "se")
  expect_equal(c(dls$lo_, dls$hi_), exp(log(2) + c(-0.1, 0.1)))                       # +/- 1 SE on the log scale
  d95 <- data.table(est = 1, lo = 0.6, hi = 1.4); add_interval_columns(d95, "est", "lo", "hi", style = "ci95"); expect_equal(c(d95$lo_, d95$hi_), c(0.6, 1.4))

  cfg <- make_config(); expect_equal(cfg$plot$interval_style, "se"); expect_equal(cfg$plot$transition_style, "box")
  cfg_box <- make_config(list(plot = list(interval_style = "box")))
  ft <- data.table(label = c("a", "b"), or = c(1.5, 1.2), or_low = c(1.2, 0.9), or_high = c(1.9, 1.6))
  geoms <- function(p) vapply(p$layers, function(l) class(l$geom)[1L], "")
  expect_true("GeomPointrange" %in% geoms(plot_forest(ft, "t"))); expect_false("GeomBoxplot" %in% geoms(plot_forest(ft, "t")))
  expect_true("GeomBoxplot" %in% geoms(plot_forest(ft, "t", style = "box")))
  expect_true("GeomPointrange" %in% geoms(plot_forest(ft, "t", style = "ci95")))
  expect_match(plot_forest(ft, "t")$labels$caption, "1 SD"); expect_match(plot_forest(ft, "t", style = "box")$labels$caption, "box: 50% interval")

  coefs <- data.table(stream = "activity", role = "primary", term = "w_act_steps", or = 1.6, or_low = 1.4, or_high = 1.8, n_events = 5L, n_persons = 3L)
  cs <- data.table(stream = "activity", cohort_stratum = "Healthy", or = 1.5, or_low = 1.2, or_high = 1.9)
  pred <- CJ(stream = "activity", adjusted = TRUE, within_z = c(-1, 0, 1), length_group = factor(c("1", "2"), ordered = TRUE))[, `:=`(prob = 0.5, prob_low = 0.4, prob_high = 0.6)]
  tr <- data.table(stream = "activity", cohort = cfg$cohorts, p_onset = 0.02, p_onset_low = 0.01, p_onset_high = 0.03,
                   p_continue = 0.5, p_continue_low = 0.4, p_continue_high = 0.6)
  inter <- data.table(stream = "activity", term = "w_act_steps", chisq = 3.2, df = 3L, p = 0.36, status = "ok")
  for (p in list(plot_onset_forest(coefs, cfg), plot_cohort_forest(cs, inter, cfg), plot_ordinal_probs(pred, cfg))) {
    expect_true("GeomPointrange" %in% geoms(p)); expect_silent(suppressWarnings(ggplot2::ggplot_build(p)))
  }
  for (p in list(plot_onset_forest(coefs, cfg_box), plot_cohort_forest(cs, inter, cfg_box), plot_ordinal_probs(pred, cfg_box), plot_transition(tr, cfg))) {
    expect_true("GeomBoxplot" %in% geoms(p)); expect_silent(suppressWarnings(ggplot2::ggplot_build(p)))      # transition stays a box by default
  }
  expect_true("GeomPointrange" %in% geoms(plot_transition(tr, make_config(list(plot = list(transition_style = "se"))))))
  txt <- vapply(plot_cohort_forest(cs, inter, cfg)$layers, function(l) if (is.null(l$data) || inherits(l$data, "waiver")) "" else paste(unlist(l$data), collapse = " "), "")
  expect_true(any(grepl("LRT chisq\\(3\\)=3.2, p=0.36", txt)))                                    # the cohort-difference test is labelled with its statistic and df
  expect_match(caption_for("cohort_forest", cfg), "1 SD of the estimate"); expect_match(caption_for("cohort_forest", cfg), "likelihood-ratio test")
  expect_match(caption_for("cohort_forest", cfg_box), "central 50% interval"); expect_match(caption_for("transition", cfg), "central 50% interval")
  expect_false(grepl("1 SD of the estimate", caption_for("gap_ecdf", cfg)))
})

test_that("one forest plot with all predictors per stream for the onset, sensitivity, B and C models, and the extra decile figure", {
  cfg <- make_config(); geoms <- function(p) vapply(p$layers, function(l) class(l$geom)[1L], "")
  built <- function(p, sized = TRUE) { expect_s3_class(p, "ggplot"); expect_silent(suppressWarnings(ggplot2::ggplot_build(p))); if (sized) expect_equal(length(attr(p, "size")), 2L) }
  terms <- c("(Intercept)", "w_act_steps", "w_act_sedentary", "b_act_steps", "b_act_sedentary", "cohortDiabetes", "age10", "sex2Male",
             "dow_tTue", "splines::ns(study_day, 3)1", "deviceVersa", "cci", "prev_gap_len_l", "any_prior_gap", "days_since_return_l")
  co <- data.table(stream = "activity", role = "primary", term = terms, estimate = 0.1, or = 1.2, or_low = 1.1, or_high = 1.3, n_events = 50L, n_persons = 30L)
  co <- rbind(co, copy(co)[, stream := "sleep"][term == "w_act_steps", term := "w_sleep_eff"])
  p <- plot_onset_forest_all(co, "activity", cfg); built(p)
  d <- ggplot2::ggplot_build(p)$plot$data
  expect_false(any(grepl("Intercept|splines|year_t", d$term)))                           # intercept and spline basis are not drawn; there are no calendar years
  expect_true(all(c("w_act_steps", "b_act_steps", "age10", "cohortDiabetes", "dow_tTue", "any_prior_gap", "deviceVersa", "cci") %in% d$term))
  expect_equal(levels(d$label)[1L], "Days since return")                 # model order, top to bottom
  expect_equal(rev(levels(d$label))[1L:2L], c("Within: Steps", "Within: Sedentary minutes"))
  expect_s3_class(p$facet, "FacetNull")                                                  # no panel strips on the right
  expect_true(all(c("Device: Versa", "CCI score") %in% levels(d$label)))
  expect_false(any(grepl("[()]", levels(d$label))))                                      # no bracketed text on the y axis
  expect_equal(nrow(d), nrow(co[stream == "activity" & !term %in% c("(Intercept)", "splines::ns(study_day, 3)1")]))   # nothing lost
  expect_true("GeomPointrange" %in% geoms(p)); expect_true("GeomBoxplot" %in% geoms(plot_onset_forest_all(co, "activity", make_config(list(plot = list(interval_style = "box"))))))
  expect_null(plot_onset_forest_all(co, "hr", cfg))
  expect_equal(term_label("cohortDiabetes"), "Cohort: Diabetes"); expect_equal(term_label("w_act_mvpa"), "Within: MVPA minutes"); expect_equal(term_label("b_sleep_eff"), "Between: Sleep efficiency")
  expect_false(any(grepl("[()]", term_label(c("log2_gap", "w_pre_day_steps", "age10", "sex2Male", "dow_tMon", "prev_gap_len_l", "days_since_return_l", "tw_hr_elev")))))

  an <- c("main", "trailing_baseline", "exclude_partial_wear", "first_gap_only", "conditional_logit", "cluster_robust_logit", "spline_vs_linear")
  at <- CJ(stream = "activity", analysis = an, term = c("w_act_steps", "w_act_sedentary"))[, `:=`(estimate = 0.1, or = 1.1, or_low = 1.0, or_high = 1.2)]
  at[analysis == "trailing_baseline", term := sub("^w_", "tw_", term)]; at[analysis == "spline_vs_linear", `:=`(estimate = NA_real_, or = NA_real_)]
  p <- plot_sensitivity_forest_stream(at, "activity", cfg); built(p)
  d <- ggplot2::ggplot_build(p)$plot$data
  expect_false("spline_vs_linear" %in% d$analysis); expect_match(levels(d$panel)[1L], "pre-specified test")
  expect_equal(uniqueN(d$panel), 2L); expect_null(plot_sensitivity_forest_stream(at, "sleep", cfg))

  ba <- CJ(stream = "activity", model_id = names(GAPB_LABELS), term = c("log2_gap", "age10", "cohortDiabetes"))[, `:=`(estimate = -0.02, se = 0.01, conf_low = -0.04, conf_high = 0)]
  ba <- rbind(ba, data.table(stream = "activity", model_id = "B_back2_adj_pre", term = "w_pre_day_steps", estimate = 0.1, se = 0.01, conf_low = 0.08, conf_high = 0.12))
  p <- plot_gapB_forest(ba, "activity", cfg); built(p); expect_equal(uniqueN(ggplot2::ggplot_build(p)$plot$data$panel), length(GAPB_LABELS))
  ca <- CJ(stream = "activity", model_id = names(GAPC_LABELS), term = c("w_pre_day_steps", "age10"))[, `:=`(estimate = 0.1, se = 0.02, or = 0.9, or_low = 0.85, or_high = 0.95)]
  p <- plot_gapC_forest(ca, "activity", cfg); built(p); expect_null(plot_gapC_forest(ca[0L], "activity", cfg))

  dec <- CJ(stream = "activity", cohort = cfg$cohorts, decile = 1:10)[, `:=`(rate = 0.05 - 0.002 * decile, rate_low = 0.04 - 0.002 * decile, rate_high = 0.06 - 0.002 * decile)]
  p <- plot_onset_deciles(dec, "activity", cfg); built(p, sized = FALSE)
  expect_true("GeomPointrange" %in% geoms(p)); expect_true("GeomBoxplot" %in% geoms(plot_onset_deciles(dec, "activity", make_config(list(plot = list(interval_style = "box"))))))
  expect_null(plot_onset_deciles(dec, "hr", cfg))
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

test_that("one direction everywhere: odds ratios are per +1 SD higher, no label or caption speaks of 'per 1 SD below'", {
  root <- Filter(function(d) file.exists(file.path(d, "R", "plots.R")), c(".", "../.."))[1L]
  skip_if(is.na(root), "package sources not found")
  files <- file.path(root, c(file.path("R", list.files(file.path(root, "R"))), "README.md", "STATISTICS.md", "MODEL_EQUATIONS.md"))
  src <- unlist(lapply(files[file.exists(files)], readLines))
  expect_false(any(grepl("1 SD below|1 SD BELOW|SD lower than|lower-than-usual|exp\\(-", src)))
  cfg <- make_config()
  expect_match(ckpt_cfg_part(cfg, "onset")$or_direction, "+1 SD higher", fixed = TRUE)            # results saved with another convention are never reused
  expect_match(ckpt_cfg_part(cfg, "gapC")$or_direction, "+1 SD higher", fixed = TRUE)
  expect_match(plot_forest(data.table(label = "a", or = 0.8, or_low = 0.7, or_high = 0.9), "t")$labels$x, "per +1 SD", fixed = TRUE)
  expect_equal(vapply(c("onset_primary_forest", "onset_forest_other_activity", "sensitivity_forest_hr", "gapC_forest_sleep"),
                      function(n) grepl("+1 SD", caption_for(n, cfg), fixed = TRUE), TRUE), setNames(rep(TRUE, 4), c("onset_primary_forest", "onset_forest_other_activity", "sensitivity_forest_hr", "gapC_forest_sleep")))
})

test_that("all-terms comparison forests: primary first (dark), every term but the intercept and the spline basis, calendar years apart", {
  cfg <- make_config(); geoms <- function(p) vapply(p$layers, function(l) class(l$geom)[1L], "")
  terms <- c("(Intercept)", "w_act_steps", "b_act_steps", "cohortDiabetes", "age10", "dow_tMon", "splines::ns(study_day, 3)1", "splines::ns(study_day, 3)2",
             "deviceVersa", "cci", "prev_gap_len_l")
  mk <- function(a) data.table(stream = "activity", analysis = a, term = terms, or = 1.2, or_low = 1.0, or_high = 1.4, estimate = 0.2, conf_low = 0.1, conf_high = 0.3)
  tab <- rbind(mk("main"), mk("trailing_baseline")[term == "w_act_steps", term := "tw_act_steps"], mk("conditional_logit")[!grepl("^b_|^cohort|^age|^device|^cci", term)])
  gl <- c(main = "Main", trailing_baseline = "Trailing baseline", conditional_logit = "Conditional logit")
  p <- plot_effects_by_group(tab, "activity", cfg, "analysis", gl, "other", "t", subtitle_or(cfg))
  d <- ggplot2::ggplot_build(p)$plot$data
  expect_false(any(grepl("Intercept|splines|year_t", d$term_key)))
  expect_equal(levels(d$grp), unname(gl)); expect_equal(unname(unlist(lapply(ggplot2::ggplot_build(p)$data[1], function(x) 0))), 0)
  expect_true(all(d[grp == "Trailing baseline", term_key] %in% terms) && "w_act_steps" %in% d[grp == "Trailing baseline", term_key])   # tw_ lines up with w_
  expect_equal(nrow(d[grp == "Conditional logit"]), 3L)                                         # only the terms that model has (w_act_steps, weekday, previous gap length)
  expect_equal(unname(p$guides$colour$params$reverse %||% TRUE), TRUE); expect_equal(length(attr(p, "size")), 2L)
  expect_true("GeomPointrange" %in% geoms(p))
  expect_equal(rev(levels(d$label))[1L], "Within: Steps")                                       # the model's own order, top to bottom
  expect_true(all(c("deviceVersa", "cci") %in% d$term_key))
  expect_null(plot_effects_by_group(tab, "sleep", cfg, "analysis", gl, "other", "t"))
  pe <- plot_effects_by_group(tab[!grepl("^b_", term)], "activity", cfg, "analysis", gl, "other", "t", measure = "estimate")   # linear model: estimates, reference 0
  expect_equal(pe$layers[[1]]$data$xintercept %||% pe$layers[[1]]$aes_params$xintercept, 0)
  expect_match(pe$labels$x, "SD units"); expect_match(p$labels$x, "Odds ratio")
  expect_silent(suppressWarnings(ggplot2::ggplot_build(p))); expect_silent(suppressWarnings(ggplot2::ggplot_build(pe)))
})

test_that("collinearity, observed post-gap, run-in and exploratory-grid figures build, with explicit sizes and labels", {
  cfg <- make_config(); geoms <- function(p) vapply(p$layers, function(l) class(l$geom)[1L], "")
  gv <- data.table(stream = "activity", analysis = rep(c("Primary onset model (gap >= 1 day)", "Model C: C_adj"), each = 3L), term = rep(c("w_act_steps", "cohort", "age10"), 2),
                   df = rep(c(1L, 3L, 1L), 2), vif_equiv = c(1.2, 6, 12, 1, 2, 1.1), flag = c("ok", "high", "severe", "ok", "ok", "ok"))
  p <- plot_collinearity(gv, "activity", cfg); expect_silent(suppressWarnings(ggplot2::ggplot_build(p)))
  expect_equal(levels(p$data$flag), c("ok", "high", "severe")); expect_equal(as.character(p$data[term == "cohort", label][1L]), "Cohort [3 df]")
  expect_null(plot_collinearity(gv, "hr", cfg))

  tr <- CJ(stream = "activity", gap_group = factor(c("1", "2", "3", "4-7", "8-14", "15+", "all gaps"), levels = c("1", "2", "3", "4-7", "8-14", "15+", "all gaps")), days_since_return = 0:7)
  tr[, `:=`(mean = 0.05 * days_since_return - 0.1, se = 0.05, n_gaps = 100L, n_persons = 60L, median_gap_length = c(1, 2, 3, 5, 10, 20, 2)[as.integer(gap_group)])]
  tr[, `:=`(lo = mean - 0.1, hi = mean + 0.1)]
  cv <- CJ(stream = "activity", model_id = "B_back2_unadj", gap_length = 1:14)[, `:=`(est = -0.01 * gap_length, lo = -0.01 * gap_length - 0.05, hi = -0.01 * gap_length + 0.05)]
  p <- plot_post_gap_observed(tr, cv, "activity", cfg); expect_silent(suppressWarnings(ggplot2::ggplot_build(p)))
  expect_true(all(c("GeomLine", "GeomPointrange") %in% geoms(p))); expect_equal(sum(geoms(p) == "GeomPointrange"), 2L)       # observed, and the model-based diamonds apart
  expect_equal(unique(p$layers[[4]]$aes_params$shape), 5)                                                                     # model estimate: open diamond
  expect_equal(nrow(p$layers[[4]]$data), 6L)                                                                                  # groups 1, 2, 3, 4-7, 8-14 and all gaps (the 15+ group is beyond the model curve)
  expect_match(levels(p$data$panel)[1L], "1 day: 100 gaps, 60 people")
  expect_equal(range(p$data$days_since_return), c(0, 7)); expect_null(plot_post_gap_observed(tr, cv, "sleep", cfg))

  # run-in figures (runin.R): hand-made result tables
  cats <- runin_category_labels(cfg)
  dist <- CJ(stream = c("activity", "sleep"), cohort = c(cfg$cohorts, "All cohorts"), category = factor(cats, levels = cats))[, `:=`(n_people = 50L, pct = 25, total_people = 200L, median_missing_days_in_category = c(0, 1, 2, 6)[as.integer(category)], mean_missing_days = 2, pct_no_missing_day = 25)]
  p <- plot_runin_distribution(dist, cfg); expect_silent(suppressWarnings(ggplot2::ggplot_build(p))); expect_null(plot_runin_distribution(dist[0L], cfg))
  rate <- CJ(stream = c("activity", "sleep"), cohort = c(cfg$cohorts, "All cohorts"), category = cats)[, `:=`(n_people = 50L, n_rows = 1000L, n_events = 30L, rate = 2 + 0.5 * match(category, cats), rate_low = 1.5 + 0.5 * match(category, cats), rate_high = 2.5 + 0.5 * match(category, cats))]
  assoc <- data.table(stream = "activity", cohort = c(cfg$cohorts, "All cohorts"), n_people = 100L, rho = 0.2, rho_low = 0.1, rho_high = 0.3, p = 0.01, p_holm = 0.03)
  p <- plot_runin_rate(rate, assoc, cfg); expect_silent(suppressWarnings(ggplot2::ggplot_build(p))); expect_match(p$labels$caption, "Spearman"); expect_null(plot_runin_rate(rate[0L], assoc, cfg))
  sg <- CJ(stream = "activity", cohort = cfg$cohorts, category = cats, bin = 1:5)[, `:=`(mean_z = (bin - 3) / 2, n_rows = 500L, n_events = 40L, n_people = 60L, rate = 3 - 0.3 * bin, rate_low = 2.5 - 0.3 * bin, rate_high = 3.5 - 0.3 * bin)]
  p <- plot_runin_steps_by_group(sg, "activity", cfg); expect_silent(suppressWarnings(ggplot2::ggplot_build(p))); expect_null(plot_runin_steps_by_group(sg, "hr", cfg))
  expect_null(plot_runin_steps_by_group(sg[, n_events := 5L], "activity", cfg))                                                   # bins with < 20 gaps are not drawn
  pe <- CJ(model = c("onset", "model C"), stream = c("activity", "sleep"), variant = c("reference", "main", "interaction"))[, `:=`(term = "w_act_steps", effect = 0.7, low = 0.6, high = 0.8, p = 0.01, status = "ok", measure = "odds ratio")]
  tt <- data.table(model = c("onset", "model C"), stream = "activity", analysis = "primary", hypothesis = "H2: steps effect depends on run-in", form = "log2 coding, with interaction",
                   p_lrt = c(0.02, NA), p_wald = c(0.02, 0.4), p = c(0.02, 0.4), p_holm = c(0.04, 0.4))
  p <- plot_runin_primary_effect(pe, tt, cfg); expect_silent(suppressWarnings(ggplot2::ggplot_build(p))); expect_match(p$labels$caption, "H2, steps x run-in interaction")
  expect_setequal(as.character(unique(p$data$variant)), RUNIN_PRIMARY_VARIANTS)                                                    # the primary figure: reference and run-in + interaction only
  ps <- plot_runin_primary_effect(pe, tt, cfg, variants = RUNIN_ONSET_VARIANTS, sensitivity = TRUE); expect_silent(suppressWarnings(ggplot2::ggplot_build(ps)))
  expect_match(ps$labels$title, "Sensitivity"); expect_setequal(as.character(unique(ps$data$variant)), c("reference", "main", "interaction"))     # every form that exists is drawn there
  eff <- CJ(model = c("onset", "model B"), stream = "activity", missing_days = c(0L, 1L, 3L, 7L))[, `:=`(estimate = 0, se = 0.1, effect = 0.7 + 0.02 * missing_days, low = 0.6, high = 0.8, measure = "x")]
  p <- plot_runin_steps_by_missing(eff, tt, cfg); expect_silent(suppressWarnings(ggplot2::ggplot_build(p))); expect_equal(unique(as.character(p$data$model)), "onset")           # odds-ratio models only
  cf <- rbind(CJ(stream = "activity", runin_variant = "categorical", term = paste0("cov0_cat", cats[-1L]))[, `:=`(or = 1.2, or_low = 1.0, or_high = 1.4)],
              data.table(stream = "activity", runin_variant = "main", term = "cov0", estimate = 0.1, se = 0.03), fill = TRUE)
  p <- plot_runin_shape(cf, NULL, dist, data.table(stream = "activity", variant = "categorical", delta_aic_vs_main = -1.5), cfg); expect_silent(suppressWarnings(ggplot2::ggplot_build(p)))
  expect_match(p$labels$caption, "AIC")
  expect_equal(term_label("w_act_steps:cov0"), "Within: Steps x run-in"); expect_equal(term_label("cov0"), "Run-in: missing days (log2)"); expect_equal(term_label("cov0_cat2-3"), "Run-in: 2-3 missing days")

  co <- CJ(stream = "activity", set_name = c("uni_act_steps", "joint"), adjusted = c(FALSE, TRUE), term = c("w_act_steps", "b_act_steps"))[
    , `:=`(k = 1L, role = ifelse(set_name == "joint" & adjusted, "primary", "exploratory"), or = 0.8, or_low = 0.7, or_high = 0.9)]
  p <- plot_exploratory_grid(co, "activity", cfg); expect_silent(suppressWarnings(ggplot2::ggplot_build(p)))
  expect_setequal(levels(p$data$spec), c("All predictors together", "Steps alone")); expect_null(plot_exploratory_grid(co[role == "primary"], "activity", cfg))

  # intervals that cannot be estimated (0 / Inf after separation) are left out: no "infinite values" warning on the log axis
  co2 <- copy(co); co2[1L, `:=`(or_low = 0, or_high = Inf)]
  p2 <- plot_exploratory_grid(co2, "activity", cfg); expect_silent(ggplot2::ggplot_build(p2)); expect_equal(nrow(p2$data), nrow(co) - 1L)
  pe2 <- copy(pe); pe2[1L, `:=`(low = 0, high = Inf)]
  p3 <- plot_runin_primary_effect(pe2, tt, cfg); expect_silent(ggplot2::ggplot_build(p3)); expect_equal(nrow(p3$data), sum(pe$variant %in% RUNIN_PRIMARY_VARIANTS) - 1L)      # an interval that cannot be estimated is not drawn
})

test_that("every figure belongs to one analysis step of the workflow and has a caption; figure_index.tsv lists them", {
  nm <- c("heatmap_sleep", "missing_pct_by_sex", "transition", "runin_distribution", "runin_rate_by_category", "onset_forest_other_activity",
          "onset_deciles_hr", "risk_curve", "gap_dose_response", "gapC_forest_hr", "post_gap_activity_observed_sleep", "runin_primary_effect",
          "runin_onset_all_terms_other_activity", "runin_gapB_all_terms_hr", "runin_steps_by_group_sleep", "runin_steps_by_missing_days", "runin_sensitivity_shape", "sensitivity_forest_hr", "sensitivity_all_terms_other_hr", "forest_sensitivity",
          "gap_definition_all_terms_other_sleep", "exploratory_grid_effects_activity", "collinearity_gvif_hr")
  steps <- vapply(nm, function(n) figure_analysis(n)$step, 1L)
  expect_false(anyNA(steps)); expect_false(anyNA(vapply(nm, caption_for, character(1))))
  expect_equal(unname(steps[c("heatmap_sleep", "transition", "risk_curve", "runin_primary_effect", "sensitivity_forest_hr", "collinearity_gvif_hr")]), c(1L, 2L, 4L, 5L, 6L, 7L))
  expect_true(is.na(figure_analysis("not_a_figure")$step))
  d <- tempfile(); dir.create(file.path(d, "figures"), recursive = TRUE); for (n in nm[1:3]) file.create(file.path(d, "figures", paste0(n, ".png")))
  write_figure_captions(nm[1:3], d, make_config())
  ix <- fread(file.path(d, "figures", "figure_index.tsv"))
  expect_equal(ix$figure, c("heatmap_sleep.png", "missing_pct_by_sex.png", "transition.png")); expect_equal(ix$step, c(1L, 1L, 2L))
  expect_true(all(c("figure", "step", "analysis", "model", "caption") %in% names(ix)))
  expect_match(caption_for("collinearity_gvif_model_B_hr"), "ONE model"); expect_match(caption_for("runin_gvif_onset_activity"), "run-in and its interaction"); expect_equal(figure_analysis("runin_gvif_onset_activity")$step, 5L); expect_equal(figure_analysis("collinearity_gvif_onset_primary_hr")$step, 7L); expect_match(caption_for("runin_sensitivity_shape"), "SENSITIVITY"); expect_match(caption_for("runin_sensitivity_gapB_all_terms_hr"), "SENSITIVITY"); expect_match(caption_for("runin_primary_effect"), "H3"); expect_match(caption_for("runin_rate_by_category"), "H1")
})

test_that("every result table is assigned to an analysis step and table_index.tsv lists them", {
  nm <- c("cohort_missingness", "missing_pct_by_stratum", "metric_summary", "transition", "weekday_tests", "device_tests", "onset_coefficients", "primary_test",
          "model_B_post_gap_activity", "model_C_gap_length", "runin_tests", "sensitivity_all_terms", "model_collinearity", "onset_term_tests", "unit_failures")
  st <- vapply(nm, function(n) table_analysis(n)$step, 1L)
  expect_equal(unname(st), c(1L, 1L, 1L, 2L, 2L, 3L, 4L, 4L, 4L, 4L, 5L, 6L, 7L, 7L, 7L))
  d <- tempfile(); dir.create(file.path(d, "tables"), recursive = TRUE)
  for (n in nm[1:4]) fwrite(data.table(a = 1), file.path(d, "tables", paste0(n, ".tsv")))
  write_table_index(d); ix <- fread(file.path(d, "tables", "table_index.tsv"))
  expect_equal(ix$table, c("cohort_missingness.tsv", "metric_summary.tsv", "missing_pct_by_stratum.tsv", "transition.tsv")); expect_equal(ix$step, c(1L, 1L, 1L, 2L))
})

test_that("the MVPA figure of models B and C compares the specifications without (dark) and with MVPA, for the three terms only", {
  cfg <- make_config()
  mk <- function(ids, terms, est) CJ(stream = "activity", model_id = ids, term = terms)[, `:=`(estimate = est, se = 0.05, conf_low = est - 0.1, conf_high = est + 0.1, or = exp(est), or_low = exp(est - 0.1), or_high = exp(est + 0.1))]
  b <- rbind(mk(c("B_back2_adj", "B_back2_adj_mvpa"), c("log2_gap", "cohortDiabetes"), 0.02), mk("B_back2_adj_mvpa", "w_pre_day_mvpa", 0.1),
             mk(c("B_back2_adj_pre", "B_back2_adj_pre_mvpa"), c("log2_gap", "w_pre_day_steps"), 0.03), mk("B_back2_adj_pre_mvpa", "w_pre_day_mvpa", 0.08), mk("B_back2_unadj", "log2_gap", 0))
  p <- plot_gap_mvpa(b, "activity", cfg, "B"); d <- ggplot2::ggplot_build(p)$plot$data
  expect_silent(suppressWarnings(ggplot2::ggplot_build(p)))
  expect_false(any(grepl("cohort|unadj", paste(d$term_key, d$model_id))))                     # only gap length, steps and MVPA; only the pairs that have an MVPA version
  expect_setequal(as.character(unique(d$model_id)), c("B_back2_adj", "B_back2_adj_mvpa", "B_back2_adj_pre", "B_back2_adj_pre_mvpa"))
  expect_equal(levels(d$grp)[1L], "Adjusted (cohort, age, sex)"); expect_true("MVPA before the gap" %in% levels(d$label))
  expect_match(p$labels$x, "Change in activity"); expect_match(p$labels$title, "MVPA Before the Gap")
  cc <- rbind(mk(c("C_adj", "C_adj_mvpa"), "w_pre_day_steps", 0.0), mk("C_adj_mvpa", "w_pre_day_mvpa", -0.1))
  pc <- plot_gap_mvpa(cc, "activity", cfg, "C"); expect_silent(suppressWarnings(ggplot2::ggplot_build(pc)))
  expect_match(pc$labels$x, "Odds ratio of a longer gap"); expect_equal(length(attr(pc, "size")), 2L)
  expect_null(plot_gap_mvpa(cc[term != "w_pre_day_mvpa"], "activity", cfg, "C")); expect_null(plot_gap_mvpa(cc, "sleep", cfg, "C")); expect_null(plot_gap_mvpa(cc[0L], "activity", cfg, "C"))
  expect_equal(term_label("w_pre_day_mvpa"), "MVPA before the gap")
  expect_false(is.na(caption_for("gapB_mvpa_effect_hr"))); expect_false(is.na(caption_for("gapC_mvpa_effect_sleep")))
  expect_equal(figure_analysis("gapB_mvpa_effect_hr")$step, 4L)
})
