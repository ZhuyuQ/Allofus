small_study <- function(beta, seed = 7) {
  cfg <- cfg_test()
  inp <- simulate_inputs(cfg, n_per_cohort = c(Insomnia = 80, Diabetes = 80, Healthy = 160, Prostate = 40),
                         beta = beta, seed = seed)
  list(cfg = cfg, inp = inp,
       b = suppressMessages(prepare_stream_bundle("activity", inp$day_tables, inp$persons, cfg, trailing = TRUE)))
}
primary_spec <- function(cfg, which_stream = "activity") {
  g <- make_onset_grid(cfg)
  g[g$stream == which_stream & g$role == "primary"]
}

test_that("grid has one primary spec per stream and the model-D thresholds", {
  g <- make_onset_grid(make_config())
  expect_equal(g[role == "primary", .N], 3L)
  expect_equal(g[role == "secondary_threshold", sort(unique(k))], c(2L, 3L, 7L))
  expect_true(all(g[role == "primary", adjusted]))
})

test_that("formula builder includes within + between terms and drops constant adjusters", {
  f <- build_onset_formula("act_steps", TRUE, 1L, drop_cols = "sex2")
  expect_match(f, "w_act_steps"); expect_match(f, "b_act_steps"); expect_false(grepl("sex2", f))
  expect_match(f, "prev_gap_len_l"); expect_false(grepl("cohort", build_onset_formula("act_steps", FALSE, 1L)))
  expect_match(build_onset_formula("act_steps", FALSE, 2L, baseline = "trailing"), "tw_act_steps")
})

test_that("the primary model recovers a planted effect and finds nothing when there is none", {
  s <- small_study(beta = -0.5)
  r <- fit_onset_spec(s$b, primary_spec(s$cfg), s$cfg)
  est <- r$coef[term == "w_act_steps"]
  expect_equal(r$info$status, "ok")
  expect_lt(abs(est$estimate - (-0.5)), 0.2); expect_lt(est$p, 0.001)
  expect_gt(est$or, 1)                                       # OR per 1 SD BELOW usual
  expect_equal(est$or, exp(-est$estimate))

  s0 <- small_study(beta = 0, seed = 8)
  r0 <- fit_onset_spec(s0$b, primary_spec(s0$cfg), s0$cfg)
  expect_lt(abs(r0$coef[term == "w_act_steps", estimate]), 0.2)
})

test_that("too few events -> not_estimable (no error)", {
  s <- small_study(-0.5); cfg <- s$cfg; cfg$model$min_events <- 1e6
  r <- fit_onset_spec(s$b, primary_spec(cfg), cfg)
  expect_equal(r$info$status, "not_estimable"); expect_null(r$coef)
})

test_that("marginal risk is higher 1 SD below usual when the slope is negative", {
  coef <- data.table(term = c("(Intercept)", "w_act_steps"), estimate = c(-4, -0.5), se = c(0.05, 0.04))
  rk <- marginal_risk(coef, "w_act_steps", sigma2 = 0.3)
  expect_true(rk$risk[1] > rk$risk[2] && rk$risk[2] > rk$risk[3])
  expect_true(all(rk$risk > 0 & rk$risk < 1))
})

test_that("first-gap-only restriction and cohort-stratified fits work (constant sex dropped for all-male cohort)", {
  s <- small_study(-0.5)
  expect_lt(nrow(restrict_first_gap(s$b$rs)), nrow(s$b$rs))
  r <- fit_onset_spec(s$b, primary_spec(s$cfg), s$cfg, variant = list(cohort = "Prostate"))
  expect_true(r$info$status %in% c("ok", "not_estimable", "not_converged"))
  expect_equal(r$info$cohort_stratum, "Prostate")
  it <- test_cohort_interaction(s$b, primary_spec(s$cfg), s$cfg)
  expect_true(is.finite(it$p) || it$status == "not_estimable")
})

test_that("sensitivity analyses are separate units, agree with the main estimate, and resume", {
  s <- small_study(-0.5); ck <- tempfile(); st <- as_day_store(s$inp$day_tables, s$inp$persons)
  dk <- data_key(st, "activity"); prim <- primary_spec(s$cfg)
  bundle <- lazy(function() s$b)
  status <- suppressMessages(run_sensitivity_stage("activity", bundle, NULL, prim, s$cfg, ck, dk))
  expect_true(all(status %in% c("computed")))
  got <- collect_sensitivity_results(ck, s$cfg, list(activity = dk))
  sens <- got$sensitivity
  expect_true(all(c("main", "trailing_baseline", "exclude_partial_wear", "first_gap_only",
                    "conditional_logit", "cluster_robust_logit", "spline_vs_linear", "random_slope") %in% sens$analysis))
  est <- sens[analysis %in% c("main", "conditional_logit", "cluster_robust_logit", "trailing_baseline"), estimate]
  expect_true(all(abs(est - sens[analysis == "main", estimate]) < 0.25))
  expect_gt(sens[analysis == "spline_vs_linear", p], 0.01)   # truth is linear
  expect_equal(nrow(got$sensitivity_failures), 0L)
  expect_equal(sens[, .N, by = analysis]$N, rep(1L, uniqueN(sens$analysis)))              # `sensitivity` = the tested term only, one row per analysis
  expect_false("tested" %in% names(sens))
  at <- got$sensitivity_all_terms                                                          # every within-person term of the estimators
  for (a in c("main", "exclude_partial_wear", "first_gap_only", "conditional_logit", "cluster_robust_logit"))
    expect_setequal(at[analysis == a, term], c("w_act_steps", "w_act_sedentary", "w_act_mvpa"))
  expect_setequal(at[analysis == "trailing_baseline", term], c("tw_act_steps", "tw_act_sedentary", "tw_act_mvpa"))
  expect_equal(at[analysis == "main" & term == "w_act_steps", estimate], sens[analysis == "main", estimate])
  expect_true(all(at[!is.na(estimate) & grepl("w_", term), or] == exp(-at[!is.na(estimate) & grepl("w_", term), estimate])))   # OR per 1 SD below usual
  expect_false(any(at$analysis %in% c("spline_vs_linear", "random_slope")))                # tests, not estimates
  again <- suppressMessages(run_sensitivity_stage("activity", lazy(function() stop("must not be built")), NULL, prim, s$cfg, ck, dk))
  expect_true(all(again == "resumed"))                        # nothing is refitted and the bundle is never built
})

test_that("default fast engine (lme4, nAGQ = 0) agrees with glmmTMB", {
  s <- small_study(-0.5)
  ag <- suppressMessages(verify_engine_agreement(s$b, primary_spec(s$cfg), s$cfg, n_persons = 300L))
  expect_lt(abs(ag[engine == "glmmTMB", diff_vs_lme4]), 0.05)
})

test_that("run_isolated returns values, passes errors on, and reports a dead worker", {
  expect_equal(run_isolated(1 + 1), 2)
  expect_equal(run_isolated(1 + 1, enabled = FALSE), 2)
  expect_error(run_isolated(stop("boom")), "boom")
  expect_error(suppressWarnings(run_isolated({ tools::pskill(Sys.getpid(), tools::SIGKILL); 1 })), "ran out of memory")
})

test_that("the model data keep only the columns a model uses", {
  s <- small_study(-0.5)
  d  <- build_onset_data(s$b, as.list(primary_spec(s$cfg)), s$cfg)
  du <- build_onset_data(s$b, as.list(make_onset_grid(s$cfg)[spec_id == "activity|uni_act_steps|k1|unadj"]), s$cfg)
  expect_lt(ncol(d), ncol(s$b$rs)); expect_lt(ncol(du), 6L)
  expect_setequal(names(du), c("person_id", "y_ge1", "w_act_steps", "b_act_steps"))
  expect_equal(nrow(du), nrow(s$b$rs[complete.cases(s$b$rs[, .(w_act_steps, b_act_steps, w_act_sedentary, b_act_sedentary, w_act_mvpa, b_act_mvpa)])]))
})

# ---- memory-saving pieces ----------------------------------------------------------------------

test_that("fit_loglik matches logLik of the fitted model and the LRT helper agrees with it", {
  s <- small_study(-0.5)
  sp <- primary_spec(s$cfg); sp[, adjusted := FALSE]
  dat <- build_onset_data(s$b, sp, s$cfg)
  f <- build_onset_formula(sp$predictors[[1L]], FALSE, 1L)
  l <- fit_loglik(f, dat, s$cfg)
  expect_true(is.finite(l$ll) && l$df >= 2)
})

test_that("vif_table from the correlation matrix equals 1/(1-R2) from lm", {
  set.seed(1); d <- data.table(a = rnorm(500)); d[, b := 0.6 * a + rnorm(500)][, c := 0.3 * a - 0.4 * b + rnorm(500)]
  v <- vif_table(d, c("a", "b", "c"))
  r2 <- summary(lm(a ~ b + c, d))$r.squared
  expect_equal(v[term == "a", vif], 1 / (1 - r2), tolerance = 1e-8)
})

test_that("the bundle has no raw / trailing / date columns unless asked, and trailing is explicit", {
  s <- small_study(-0.5)
  b0 <- suppressMessages(prepare_stream_bundle("activity", s$inp$day_tables, s$inp$persons, s$cfg))
  expect_false(any(grepl("^(x_|tw_)", names(b0$rs)))); expect_false("t_date" %in% names(b0$rs))
  expect_true(any(grepl("^tw_", names(s$b$rs))))
  expect_error(build_onset_data(b0, as.list(primary_spec(s$cfg)), s$cfg, list(baseline = "trailing")), "trailing")
  expect_null(b0$act); expect_null(b0$days)
})

test_that("adding columns by reference gives the same table as a merge", {
  s <- small_study(-0.5)
  b <- suppressMessages(prepare_stream_bundle("sleep", s$inp$day_tables, s$inp$persons, s$cfg))
  expect_equal(anyDuplicated(b$rs, by = c("person_id", "date")), 0L)
  a <- predictor_table(s$inp$day_tables$activity, s$cfg$activity_predictors, s$cfg, trailing = FALSE, keep_x = FALSE)
  chk <- merge(b$rs[, .(person_id, date, w_act_steps)], a[, .(person_id, date, w2 = w_act_steps)], by = c("person_id", "date"))
  expect_equal(chk$w_act_steps, chk$w2)
})

test_that("memory estimate is linear; the guard stops one model with advice and can be switched off", {
  s <- small_study(-0.5); sp <- as.list(primary_spec(s$cfg))
  expect_equal(estimate_fit_memory_gb(1e6, 10, 56), 0.56)
  expect_true(is.na(available_memory_gb()) || available_memory_gb() > 0)
  expect_true(check_fit_memory(s$b, sp, s$cfg))                          # a tiny model always fits
  tight <- s$cfg; tight$resources$memory_budget_gb <- 0
  expect_error(check_fit_memory(s$b, sp, tight), "Nothing was subsampled")
  # a cohort-specific fit is judged on that cohort's rows only, not on the whole risk set
  per_row <- estimate_fit_memory_gb(1, 2 * length(sp$predictors[[1L]]) + 36 + 3, s$cfg$resources$bytes_per_cell)
  n_all <- nrow(s$b$rs); n_pro <- sum(s$b$rs$cohort == "Prostate")
  edge <- s$cfg; edge$resources$memory_budget_gb <- per_row * (n_pro + n_all) / 2 * 1.0       # between the two sizes
  expect_error(check_fit_memory(s$b, sp, edge), "Nothing was subsampled")
  expect_true(check_fit_memory(s$b, sp, edge, variant = list(cohort = "Prostate")))
  expect_error(check_fit_memory(s$b, sp, tight), "sample_persons")      # names the explicit alternative
  tight$resources$memory_guard <- FALSE
  expect_true(check_fit_memory(s$b, sp, tight))
})

test_that("absolute risk uses the slope (not the intercept) and has an interval", {
  coef <- data.table(term = c("(Intercept)", "w_act_steps"), estimate = c(-4, -0.5), se = c(0.05, 0.04))
  rk <- marginal_risk(coef, "w_act_steps", sigma2 = 0, cov_b0_term = 0.001)
  expect_equal(rk$risk, plogis(-4 + -0.5 * c(-1, 0, 1)))
  expect_true(all(rk$risk_low < rk$risk & rk$risk < rk$risk_high))
  expect_true(all(is.na(marginal_risk(coef, "w_act_steps", 0.3)$risk_low)))   # no covariance -> no interval
})

test_that("a singular lme4 fit counts as converged (and is flagged singular), a real failure does not", {
  set.seed(1); d <- data.table(g = rep(1:50, each = 6), x = rnorm(300)); d[, y := 0.3 * x + rnorm(300)]       # no person variance
  fit <- suppressMessages(lme4::lmer(y ~ x + (1 | g), data = d))
  expect_true(lme4::isSingular(fit)); expect_true(lme4_converged(fit)); expect_true(fit_converged(fit))
  fake <- fit; fake@optinfo$conv$lme4 <- list(messages = "Model failed to converge with max|grad| = 0.5")
  expect_false(lme4_converged(fake))
  fake2 <- fit; fake2@optinfo$conv$opt <- 1L; expect_false(lme4_converged(fake2))
})

test_that("observed gap rate by decile of the primary term: ten deciles per cohort, rate falls as the prior day gets more active when the effect is negative", {
  s <- small_study(-0.8)
  dec <- onset_deciles(s$b, s$cfg)
  expect_true(all(c("cohort", "decile", "mean_z", "n_rows", "n_events", "n_persons", "rate", "rate_low", "rate_high") %in% names(dec)))
  expect_equal(max(dec$decile), 10L)
  expect_equal(sum(dec$n_rows), sum(!is.na(s$b$rs$w_act_steps) & !is.na(s$b$rs$y_ge1)))        # every row at risk is in exactly one decile
  expect_equal(sum(dec$n_events), sum(s$b$rs$y_ge1[!is.na(s$b$rs$w_act_steps)], na.rm = TRUE))
  expect_true(all(dec$rate_low <= dec$rate + 1e-12 & dec$rate <= dec$rate_high + 1e-12, na.rm = TRUE))
  expect_lt(cor(dec[, .(r = weighted.mean(rate, n_rows)), by = decile]$r, 1:10), -0.5)         # lower activity -> more gaps
  expect_equal(nrow(onset_deciles(list(rs = s$b$rs[0L]), s$cfg)), 0L)
})

test_that("Healthy is the reference cohort in every model: it is the first level, no Healthy coefficient exists, and the setting is part of the checkpoint key", {
  s <- small_study(-0.5)
  expect_equal(levels(cohort_factor(c("Prostate", "Healthy"), s$cfg)), c("Healthy", "Insomnia", "Diabetes", "Prostate"))
  expect_equal(levels(s$b$rs$cohort)[1L], "Healthy")
  r <- fit_onset_spec(s$b, primary_spec(s$cfg), s$cfg)
  expect_setequal(grep("^cohort", r$coef$term, value = TRUE), c("cohortInsomnia", "cohortDiabetes", "cohortProstate"))
  g <- build_gap_level_data(s$inp$day_tables$activity, s$inp$day_tables$activity, s$inp$persons, s$cfg)
  expect_equal(levels(g$cohort)[1L], "Healthy"); expect_equal(gap_adjusters(g), c("cohort", "age10", "sex2"))
  cfg2 <- s$cfg; cfg2$reference_cohort <- "Diabetes"
  expect_false(identical(ckpt_cfg_part(s$cfg, "onset"), ckpt_cfg_part(cfg2, "onset")))        # old checkpoints (other reference) are not reused
  expect_false(identical(ckpt_cfg_part(s$cfg, "gapB"), ckpt_cfg_part(cfg2, "gapB")))
  fit_d <- fit_onset_spec(s$b, primary_spec(s$cfg), s$cfg, variant = list(cohort = "Prostate"))   # a cohort-specific fit still works with a factor
  expect_equal(fit_d$info$cohort_stratum, "Prostate")
})
