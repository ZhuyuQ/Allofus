# Comparison of two finished runs on different people (pipeline/compare_bands.R).
if (!exists("cmp_run")) source(file.path(if (file.exists("pipeline/compare_bands.R")) "." else "../..", "pipeline/compare_bands.R"))

test_that("the ratio of odds ratios is the second band relative to the first, with the SE of two independent estimates", {
  prim <- data.table(band = factor(rep(c("A", "B"), each = 2), levels = c("A", "B")), stream = rep(c("activity", "sleep"), 2),
                     estimate = c(-0.50, -0.40, -0.30, -0.45), se = c(0.05, 0.06, 0.12, 0.10), n_events = c(900L, 800L, 150L, 140L))
  prim[, `:=`(or = exp(estimate), or_low = exp(estimate - 1.96 * se), or_high = exp(estimate + 1.96 * se))]
  sds <- data.table(band = factor(rep(c("A", "B"), each = 3), levels = c("A", "B")), stream = rep(names(CMP_STREAM_LABELS), 2), sd_within = c(.30, .25, .20, .36, .30, .24))
  r <- cmp_ror(prim, c("A", "B"), sds)
  expect_equal(r$stream, c("activity", "sleep"))
  expect_equal(r$log_ratio, c(-0.30 - -0.50, -0.45 - -0.40)); expect_equal(r$ratio, exp(r$log_ratio))
  expect_equal(r$se_log_ratio, c(sqrt(0.05^2 + 0.12^2), sqrt(0.06^2 + 0.10^2)))
  expect_equal(r$ratio_low, exp(r$log_ratio - 1.96 * r$se_log_ratio)); expect_equal(r$ratio_high, exp(r$log_ratio + 1.96 * r$se_log_ratio))
  expect_equal(c(r$within_sd_ref, r$within_sd_cmp), c(.30, .25, .36, .30))
  expect_equal(c(r$reference[1], r$compared[1]), c("A", "B")); expect_false(any(grepl("^p", names(r))))     # no p-value by design
  expect_s3_class(plot_cmp_primary(prim, cfg_test(), sds), "ggplot")
})

test_that("the within-person SD is the one the models standardize with, per stream", {
  s <- small_study(-0.5, seed = 3); cfg <- s$cfg
  st <- as_day_store(s$inp$day_tables, s$inp$persons)
  sd <- compare_within_sd(st, cfg)
  expect_equal(sd$stream, c("activity", "sleep", "hr"))
  b <- suppressMessages(prepare_stream_bundle("sleep", s$inp$day_tables, s$inp$persons, cfg))
  expect_equal(sd[stream == "sleep", sd_within], b$scales[predictor == "act_steps", sd_within])         # the same constant as in the bundle of that stream
  a <- suppressMessages(prepare_stream_bundle("activity", s$inp$day_tables, s$inp$persons, cfg))
  expect_equal(sd[stream == "activity", sd_within], a$scales[predictor == "act_steps", sd_within])
})

test_that("two finished runs on different people are compared: tables, ratio and every figure", {
  cfg <- cfg_test()
  run <- function(seed, n) {
    inp <- suppressMessages(simulate_inputs(cfg, n_per_cohort = n, seed = seed)); out <- tempfile()
    invisible(suppressWarnings(suppressMessages(capture.output({
      run_missingness_pipeline(inp$day_tables, inp$persons, cfg, out, stages = c("missingness", "temporal"))
      run_missingness_pipeline(inp$day_tables, inp$persons, cfg, out, stages = "onset", onset_roles = "primary", model_streams = "activity")
      run_missingness_pipeline(inp$day_tables, inp$persons, cfg, out, stages = "sensitivity", model_streams = "activity")}))))
    list(inp = inp, out = out)
  }
  a <- run(5, c(Insomnia = 50, Diabetes = 50, Healthy = 100, Prostate = 30)); b <- run(9, c(Insomnia = 40, Diabetes = 40, Healthy = 80, Prostate = 30))
  out_dirs <- list(`>=70%` = a$out, `50-60%` = b$out); cd <- tempfile()
  stores <- list(`>=70%` = as_day_store(a$inp$day_tables, a$inp$persons), `50-60%` = as_day_store(b$inp$day_tables, b$inp$persons))
  r <- suppressWarnings(suppressMessages(cmp_run(out_dirs, cfg, cd, stores)))
  tabs <- c("compare_within_sd", "compare_primary", "compare_primary_ratio", "compare_sensitivity", "compare_missingness")
  figs <- c("compare_primary_effect", "compare_sensitivity", "compare_missing_pct", "compare_longest_gap", "compare_gaps_per_person", "compare_gap_length_distribution", "compare_transition")
  expect_true(all(file.exists(file.path(cd, "tables", paste0(tabs, ".tsv"))))); expect_true(all(file.exists(file.path(cd, "figures", paste0(figs, ".png")))))
  # the primary table holds the primary term of each run, as saved by that run
  for (bd in names(out_dirs)) {
    saved <- fread(file.path(out_dirs[[bd]], "tables", "primary_test.tsv"))[stream == "activity"]
    expect_equal(r$primary[band == bd & stream == "activity", estimate], saved$estimate)
    expect_equal(r$primary[band == bd & stream == "activity", n_events], saved$n_events)
  }
  expect_equal(r$ratio[stream == "activity", ratio], exp(r$primary[band == "50-60%", estimate] - r$primary[band == ">=70%", estimate]))
  expect_true(all(c("wear_ge16h", "wear_ge20h", "wear_ge22h", "day_before_gap") %in% r$sensitivity[band == "50-60%", analysis]))      # every analysis of both runs is in the comparison
  expect_equal(nlevels(r$missingness$band), 2L); expect_equal(nrow(r$missingness), 3L * 4L * 2L)                # 3 streams x 4 cohorts x 2 bands
  expect_error(cmp_read(list(x = file.path(tempdir(), "nope")), "sensitivity"), "not found")
})
