ov_cfg <- function() { cfg <- cfg_test(); cfg$outcome_vars$boot_B <- 60L; cfg$device$mc_B <- 300L; cfg$device$boot_B <- 50L; cfg }

# one row per person with a planted association: Diabetes misses more days, and the missing share rises with x
ov_synth <- function(cfg, n = 1200, seed = 3) {
  set.seed(seed)
  d <- data.table(person_id = seq_len(n), cohort = factor(sample(cfg$cohorts, n, TRUE, prob = c(.2, .2, .5, .1)), levels = cfg$cohorts),
                  sex = factor(sample(c("Female", "Male"), n, TRUE), levels = c("Female", "Male", "Other")), age_first = round(runif(n, 25, 80)),
                  watch = factor(sample(c("Missing", "Charge", "Versa"), n, TRUE), levels = c("Missing", "Charge", "Versa")),
                  act_steps = rnorm(n, 8.8, 0.4))
  d[, act_steps__display := expm1(act_steps)]
  d[, cci_score := stats::rpois(.N, c(Insomnia = 3, Diabetes = 5, Healthy = 1, Prostate = 4)[as.character(cohort)])]
  d[, cci_group := cci_group_factor(as.character(cut(cci_score, c(-Inf, 0, 2, 4, Inf), c("0", "1-2", "3-4", ">=5"))))]
  p <-plogis(-3 + 0.7 * (d$cohort == "Diabetes") + 0.5 * scale(d$act_steps)[, 1])
  d[, `:=`(n_observed = 200L, window_days = 200L)]
  d[, n_missing := stats::rbinom(.N, 200L, p)]; d[, n_observed := 200L - n_missing]; d[, pct_missing := 100 * n_missing / 200]
  d[]
}

test_that("person means of the model predictors: log1p predictors are shown as geometric means, identity ones as means", {
  cfg <- ov_cfg()
  tbl <- data.table(person_id = c(1, 1, 1, 2, 2), date = as.Date("2022-01-01") + c(0:2, 0:1), steps = c(1000, 3000, NA, 500, 500),
                    sedentary_minutes = c(600, 700, 800, 900, 1000), mvpa_minutes = c(10, 20, 30, 0, 0))
  m <- person_predictor_means(list(activity = tbl), c("act_steps", "act_sedentary"), cfg)
  expect_equal(m[person_id == 1, act_steps], mean(log1p(c(1000, 3000))))
  expect_equal(m[person_id == 1, act_steps__display], expm1(mean(log1p(c(1000, 3000)))))
  expect_equal(m[person_id == 2, act_steps__display], 500)
  expect_equal(m[person_id == 1, act_sedentary__display], 700)
})

test_that("Spearman interval and the bootstrap effect-size interval behave", {
  ci <- spearman_ci(0.3, 400); expect_lt(ci[1L], 0.3); expect_gt(ci[2L], 0.3); expect_true(all(abs(ci) < 1))
  expect_true(all(is.na(spearman_ci(NA_real_, 400))))
  set.seed(1); v <- c(rnorm(80, 0), rnorm(80, 1), rnorm(80, 0)); g <- factor(rep(c("a", "b", "c"), each = 80))
  e <- boot_rank_effect(v, g, 80); H <- unname(stats::kruskal.test(v, g)$statistic) / (length(v) - 1)
  expect_lt(e[1L], H); expect_gt(e[2L], H)
  e2 <- boot_rank_effect(v[g != "c"], droplevels(g[g != "c"]), 80); expect_true(all(abs(e2) <= 1))                  # two groups: Cliff's delta
})

test_that("the analysis finds the planted effects, adjusts them, compares each level with the reference and corrects over the variables of the stream", {
  cfg <- ov_cfg(); d <- ov_synth(cfg)
  spec <- list(act_steps = list(label = "Steps per day, geometric mean", display = "act_steps__display", log_x = TRUE))
  a <- outcome_variable_analysis(d, "activity", cfg, spec)
  tt <- a$tests
  expect_setequal(tt$variable, c("cohort", "sex", "cci_group", "watch", "age_first", "cci_score", "act_steps"))
  expect_equal(unique(a$contrasts[variable == "cci_group", reference]), as.character(d[, .N, by = cci_group][which.max(N), cci_group]))   # reference = the most common CCI group
  expect_true(all(tt$p_holm >= tt$p - 1e-12) && all(tt$p_holm <= 1))
  expect_lt(tt[variable == "cohort", p_holm], 0.001); expect_lt(tt[variable == "act_steps", p_holm], 0.001)
  expect_gt(tt[variable == "sex", p_holm], 0.01)                                                              # nothing planted
  expect_gt(tt[variable == "act_steps", effect], 0.1)                                                         # positive Spearman rho
  expect_true(all(tt$effect_low <= tt$effect & tt$effect <= tt$effect_high, na.rm = TRUE))
  ct <- a$contrasts[variable == "cohort"]
  expect_equal(unique(ct$reference), cfg$reference_cohort); expect_gt(ct[level == "Diabetes", hl_difference_pct_points], 0.5)
  expect_true(all(ct$p_holm >= ct$p - 1e-12))
  m <- a$models
  expect_true(all(c("unadjusted", "adjusted for sex, age10") %in% m[variable == "cohort", adjustment]))
  expect_gt(m[variable == "cohort" & term == "Diabetes" & adjustment == "unadjusted", odds_ratio], 1.3)
  s1 <- m[variable == "act_steps" & adjustment != "unadjusted"]; expect_equal(s1$term, "per +1 SD"); expect_gt(s1$or_low, 1)
  ce <- a$cells[variable == "act_steps"]
  expect_lte(uniqueN(ce$key), cfg$outcome_vars$bins); expect_equal(sum(ce$n), nrow(d))                  # every person is in exactly one cohort x bin cell
  expect_true(all(ce$low <= ce$mean_pct & ce$mean_pct <= ce$high))
  expect_true(all(ce$w_low <= ce$q25 & ce$q25 <= ce$median & ce$median <= ce$q75 & ce$q75 <= ce$w_high))
  expect_equal(a$summary[variable == "cohort", sum(n)], nrow(d))
  # small levels are left out of the test and said so
  d2 <- copy(d); d2[1:3, sex := "Other"]
  b <- outcome_variable_analysis(d2, "activity", cfg, spec); expect_match(b$tests[variable == "sex", note], "Other")
})

test_that("per-cohort correlations and the cohort-difference test are produced; a planted cohort-specific slope is detected", {
  cfg <- ov_cfg(); set.seed(5); n <- 1600
  d <- data.table(person_id = seq_len(n), cohort = factor(sample(cfg$cohorts, n, TRUE, prob = c(.25, .25, .4, .1)), levels = cfg$cohorts),
                  sex = factor(sample(c("Female", "Male"), n, TRUE), levels = c("Female", "Male", "Other")), age_first = round(runif(n, 25, 80)), x = rnorm(n))
  d[, x__display := x]
  slope <- c(Insomnia = 0.6, Diabetes = -0.6, Healthy = 0, Prostate = 0)[as.character(d$cohort)]          # the association differs by cohort
  d[, n_missing := stats::rbinom(.N, 200L, plogis(-3 + slope * x))]; d[, `:=`(n_observed = 200L - n_missing, pct_missing = 100 * n_missing / 200)]
  a <- outcome_variable_analysis(d, "activity", cfg, list(x = list(label = "x", display = "x__display", log_x = FALSE)))
  tt <- a$tests[variable == "x"]
  expect_lt(tt$p_cohort_difference, 1e-6); expect_gte(tt$p_cohort_difference_holm, tt$p_cohort_difference)
  bc <- a$by_cohort[variable == "x"]; expect_equal(nrow(bc), 4L)
  expect_gt(bc[cohort == "Insomnia", rho], 0.2); expect_lt(bc[cohort == "Diabetes", rho], -0.2)
  expect_true(!is.na(a$tests[variable == "sex", p_cohort_difference])); expect_true(is.na(a$tests[variable == "cohort", p_cohort_difference]))   # no such test for the cohort variable itself
  expect_equal(a$cells[variable == "sex", sum(n)], nrow(d))
})

test_that("categorical variables are violins and continuous variables are lines with a band; the statistics are written in the caption below the figure", {
  cfg <- ov_cfg(); d <- ov_synth(cfg)
  a <- outcome_variable_analysis(d, "activity", cfg, list(act_steps = list(label = "Steps per day, geometric mean", display = "act_steps__display", log_x = TRUE)))
  ovs <- list(activity = a)
  pc <- plot_outcome_variable(ovs, "sex", cfg); pn <- plot_outcome_variable(ovs, "act_steps", cfg)
  has <- function(p, g) any(vapply(p$layers, function(l) inherits(l$geom, g), logical(1)))
  expect_true(has(pc, "GeomViolin")); expect_false(has(pc, "GeomBoxplot")); expect_false(has(pc, "GeomLine"))
  expect_true(has(pn, "GeomRibbon") && has(pn, "GeomLine")); expect_false(has(pn, "GeomViolin") || has(pn, "GeomBoxplot"))
  for (p in list(pc, pn)) {
    invisible(ggplot2::ggplot_build(p))
    expect_match(p$labels$caption, "Statistics"); expect_match(p$labels$caption, "p\\(Holm\\)")                  # the tests are in the caption ...
    expect_null(p$labels$subtitle)                                                                              # ... not in the panels or subtitle
    expect_false(any(vapply(p$layers, function(l) inherits(l$geom, "GeomText"), logical(1))))
  }
  expect_match(pn$labels$caption, "rho within cohort"); expect_match(pc$labels$caption, "Kruskal-Wallis|Wilcoxon")
})

test_that("ONE figure per variable builds and the stage writes them, the tables and the index entries (with and without device information)", {
  cfg <- ov_cfg()
  s <- small_study(-0.5); persons <- s$inp$persons
  tabs <- build_device_tables(sim_device_rows(persons, cfg), persons, cfg)
  out <- file.path(tempdir(), "ov_stage"); unlink(out, recursive = TRUE)
  r <- suppressMessages(run_missingness_pipeline(as_day_store(s$inp$day_tables, persons, device = tabs), persons, cfg, out, stages = "variables"))
  vars <- unique(r$outcome_variable_tests$variable)
  expect_true(all(c("cohort", "sex", "cci_group", "watch", "age_first", "cci_score", "act_steps", "sleep_inbed", "hr_elev") %in% vars))
  expect_true(all(file.exists(file.path(out, "figures", paste0("outcome_by_variable_", vars, ".png")))))              # ONE figure per variable
  expect_false(file.exists(file.path(out, "figures", "outcome_by_variable_categorical.png")))
  expect_true(all(file.exists(file.path(out, "tables", paste0(c("outcome_variable_tests", "outcome_variable_contrasts", "outcome_variable_models", "outcome_variable_summary", "outcome_variable_cells", "outcome_variable_by_cohort"), ".tsv")))))
  expect_equal(r$outcome_variable_tests[variable == "act_steps", sort(stream)], c("activity", "hr", "sleep"))        # shared variables appear in every stream's column
  expect_equal(r$outcome_variable_tests[variable == "sleep_inbed", stream], "sleep")
  fi <- fread(file.path(out, "figures", "figure_index.tsv")); expect_true(all(fi[grepl("^outcome_by_variable", figure), step] == 1L))
  expect_match(paste(readLines(file.path(out, "figures", "figure_captions.md")), collapse = " "), "ONE model variable")
  # no device information: the watch version is simply left out
  out2 <- file.path(tempdir(), "ov_stage2"); unlink(out2, recursive = TRUE)
  r2 <- suppressMessages(run_missingness_pipeline(as_day_store(s$inp$day_tables, persons), persons, cfg, out2, stages = "variables"))
  expect_false("watch" %in% r2$outcome_variable_tests$variable); expect_true("cohort" %in% r2$outcome_variable_tests$variable)
})
