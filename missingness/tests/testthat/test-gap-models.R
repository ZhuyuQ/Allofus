test_that("sparse length groups merge into the next shorter group", {
  f <- factor(c(rep("1", 100), rep("2", 50), rep("3", 30), rep("4-7", 10), rep("8-14", 5)),
              levels = c("1", "2", "3", "4-7", "8-14"), ordered = TRUE)
  g <- collapse_sparse_bins(f, 20)
  expect_equal(levels(g), c("1", "2", "3+4-7+8-14")); expect_equal(as.integer(table(g)), c(100, 50, 45))
})

gap_fixture <- function(seed = 4, n = c(Insomnia = 40, Diabetes = 40, Healthy = 80, Prostate = 20)) {
  cfg <- cfg_test(); inp <- suppressMessages(simulate_inputs(cfg, n_per_cohort = n, seed = seed))
  list(cfg = cfg, inp = inp, store = as_day_store(inp$day_tables, inp$persons))
}
gap_data_for <- function(fx, stream) with(fx, build_gap_level_data(inp$day_tables[[stream]], inp$day_tables$activity, inp$persons, cfg))

test_that("gap-level data lines up with the activity z-scores before and after each gap, and keeps only needed columns", {
  fx <- gap_fixture(); g <- gap_data_for(fx, "sleep")
  expect_equal(nrow(g), nrow(get_gap_runs(fx$inp$day_tables$sleep)))
  expect_false(any(c("prev_obs_date", "next_obs_date", "gap_start", "gap_end", "gap_index") %in% names(g)))
  a <- predictor_table(fx$inp$day_tables$activity, "act_steps", fx$cfg, trailing = FALSE, keep_x = TRUE, keep_baseline = TRUE)
  sdw <- attr(a, "scales")$sd_within
  gr <- get_gap_runs(fx$inp$day_tables$sleep)
  i <- which(!is.na(g$w_pre_day_steps) & !is.na(g$back2_dev))[1]
  key <- g[i, .(person_id, gap_length)]
  row <- gr[person_id == key$person_id & gap_length == key$gap_length][1]
  pre <- a[person_id == row$person_id & date == row$prev_obs_date]
  expect_equal(g$w_pre_day_steps[i], pre$w_act_steps)                                     # the day before the gap: against the usual BEFORE that day
  # the days after the gap: against the usual THROUGH the last day before the gap (data observed before the gap only)
  expect_equal(g$back2_dev[i], (a[person_id == row$person_id & date == row$next_obs_date + 1L, x_act_steps] - pre$m_act_steps) / sdw)
  expect_equal(g$back1_dev[i], g$post_k0[i]); expect_equal(g$back2_dev[i], g$post_k1[i])
  expect_equal(g$log2_gap[i], log2(g$gap_length[i]))
  expect_true(all(paste0("post_k", 0:7) %in% names(g)))
})

test_that("gap data from the activity-only predictor equal those from the full bundle's activity predictors", {
  fx <- gap_fixture(); b <- suppressMessages(prepare_stream_bundle("sleep", fx$inp$day_tables, fx$inp$persons, fx$cfg))
  g <- gap_data_for(fx, "sleep"); gr <- get_gap_runs(fx$inp$day_tables$sleep)
  a <- b$rs[, .(person_id, date, w = w_act_steps)]; setkey(a, person_id, date)
  expect_equal(sort(g$w_pre_day_steps), sort(a[.(gr$person_id, gr$prev_obs_date), w]))
  # the full-period variant: the days after a gap are the plain whole-period z-scores
  gf <- build_gap_level_data(fx$inp$day_tables$sleep, fx$inp$day_tables$activity, fx$inp$persons, fx$cfg, standardization = "full_period")
  full <- predictor_table(fx$inp$day_tables$activity, "act_steps", fx$cfg, trailing = FALSE, keep_x = FALSE, standardization = "full_period")
  setkey(full, person_id, date)
  expect_equal(sort(gf$back2_dev), sort(full[.(gr$person_id, gr$next_obs_date + 1L), w_act_steps]))
  expect_false("post_k0" %in% names(gf))
  expect_false(isTRUE(all.equal(sort(g$back2_dev), sort(gf$back2_dev))))                  # the two standardizations are different quantities
})

test_that("NO LOOK-AHEAD in the gap table: days after the post-gap window never change a gap's z-scores, days before the gap do", {
  cfg <- cfg_test(); set.seed(2)
  days <- data.table(person_id = rep(1:20, each = 80), day = rep(1:80, 20))[!(day %in% 40:42)]       # every person: one gap of 3 days (days 40-42)
  days[, `:=`(date = as.Date("2024-01-01") + day - 1L, steps = round(runif(.N, 3000, 9000)), sedentary_minutes = 800, mvpa_minutes = 10)]
  setkey(days, person_id, date)
  persons <- data.table(person_id = 1:20, cohort = rep(c("Healthy", "Diabetes"), 10), sex = rep(c("Female", "Male"), each = 10), age_first = 50, cci_score = rep(0:3, 5))
  gaps_of <- function(d) build_gap_level_data(d, d, persons, cfg)[order(person_id)]
  cols <- c("w_pre_day_steps", paste0("post_k", 0:7))
  g1 <- gaps_of(days)
  late <- copy(days); late[person_id == 7L & day > 50L, steps := steps * 5]                          # return day is day 43; the window ends on day 50
  early <- copy(days); early[person_id == 7L & day < 40L, steps := steps * 5]                        # before the gap
  ratio <- function(g2) { m <- as.matrix(g1[person_id == 7L, ..cols]) / as.matrix(g2[person_id == 7L, ..cols]); m[is.finite(m)] }
  r_late <- ratio(gaps_of(late)); r_early <- ratio(gaps_of(early))
  expect_lt(diff(range(r_late)), 1e-9)          # only the global scale constant moves (one common factor); the person's own deviations are unchanged
  expect_gt(diff(range(r_early)), 1e-3)         # a change before the gap does move the z-scores
})

test_that("the observed post-gap trajectory is a clustered mean by gap-length group and day since return, with sample sizes", {
  fx <- gap_fixture(); g <- gap_data_for(fx, "sleep")
  tr <- observed_post_gap_trajectory(g, fx$cfg)
  expect_setequal(levels(tr$gap_group), c("1", "2", "3", "4-7", "8-14", "15+", "all gaps"))
  expect_setequal(unique(tr$days_since_return), 0:7)
  a <- tr[gap_group == "all gaps" & days_since_return == 1L]
  expect_equal(a$mean, mean(g$post_k1, na.rm = TRUE)); expect_equal(a$n_gaps, sum(!is.na(g$post_k1)))
  expect_equal(a$n_persons, uniqueN(g$person_id[!is.na(g$post_k1)]))
  expect_true(all(tr$lo <= tr$mean & tr$mean <= tr$hi, na.rm = TRUE))
  expect_error(observed_post_gap_trajectory(g[, !grepl("^post_k", names(g)), with = FALSE], fx$cfg), "post-gap columns")
})

test_that("each model B / C is its own checkpointed unit; resuming skips finished ones; one failure leaves the others", {
  fx <- gap_fixture(); ck <- tempfile(); dk <- data_key(fx$store, "sleep")
  built <- 0; gd <- lazy(function() { built <<- built + 1; gap_data_for(fx, "sleep") })
  st1 <- suppressWarnings(suppressMessages(run_gap_stage("sleep", gd, fx$cfg, ck, dk)))
  no_full <- setdiff(names(gap_unit_specs()), c("B_back2_adj_fullperiod", "C_adj_fullperiod"))     # run only with gap_data_full
  expect_setequal(names(st1), no_full); expect_true(all(st1 == "computed")); expect_equal(built, 1)
  files <- list.files(ck); expect_length(files, length(no_full))

  f <- file.path(ck, files[grepl("^gapC__sleep__C_adj__", files)]); unlink(f)           # one unit lost
  t0 <- file.mtime(file.path(ck, files[grepl("^gapB__sleep__B_back2_adj__", files)]))
  st2 <- suppressWarnings(suppressMessages(run_gap_stage("sleep", gd, fx$cfg, ck, dk)))
  expect_equal(unname(st2[c("C_adj")]), "computed"); expect_true(all(st2[names(st2) != "C_adj"] == "resumed"))
  expect_equal(file.mtime(file.path(ck, files[grepl("^gapB__sleep__B_back2_adj__", files)])), t0)

  gd2 <- lazy(function() stop("must not be built when everything is finished"))
  expect_true(all(suppressMessages(run_gap_stage("sleep", gd2, fx$cfg, ck, dk)) == "resumed"))

})

test_that("a failure in model C leaves model B untouched, is recorded, and is retried on the next call", {
  fx <- gap_fixture(); ck <- tempfile(); dk <- data_key(fx$store, "sleep")
  broken <- lazy(function() { g <- gap_data_for(fx, "sleep"); g[, length_bin := NULL]; g })          # only C needs length_bin
  st <- suppressWarnings(suppressMessages(run_gap_stage("sleep", broken, fx$cfg, ck, dk)))
  expect_true(all(st[grepl("^B_", names(st))] == "computed")); expect_true(all(st[grepl("^C_", names(st))] == "failed"))
  r <- collect_gap_results(ck, fx$cfg, list(sleep = dk))
  expect_equal(nrow(r$model_B_post_gap_activity), 6L); expect_equal(nrow(r$model_C_gap_length), 0L)     # the six B specifications (incl. the two with MVPA)
  expect_setequal(r$gap_failures$id, c("C_unadj", "C_adj", "C_adj_mvpa")); expect_true(all(nzchar(r$gap_failures$error)))

  st2 <- suppressWarnings(suppressMessages(run_gap_stage("sleep", lazy(function() gap_data_for(fx, "sleep")), fx$cfg, ck, dk)))
  expect_true(all(st2[grepl("^B_", names(st2))] == "resumed")); expect_true(all(st2[grepl("^C_", names(st2))] == "computed"))
  r2 <- collect_gap_results(ck, fx$cfg, list(sleep = dk))
  expect_equal(nrow(r2$gap_failures), 0L); expect_equal(nrow(r2$model_C_gap_length), 3L)
})

test_that("results collect into tables, curves and predicted probabilities that are coherent", {
  fx <- gap_fixture(); ck <- tempfile(); dk <- data_key(fx$store, "sleep")
  suppressWarnings(suppressMessages(run_gap_stage("sleep", lazy(function() gap_data_for(fx, "sleep")), fx$cfg, ck, dk)))
  r <- collect_gap_results(ck, fx$cfg, list(sleep = dk))
  expect_setequal(r$model_B_post_gap_activity$model_id, c("B_back2_unadj", "B_back2_adj", "B_back2_adj_pre", "B_back1_adj", "B_back2_adj_mvpa", "B_back2_adj_pre_mvpa"))
  expect_true(all(c("converged", "singular", "status", "note") %in% names(r$model_B_post_gap_activity)))
  cv <- r$model_B_curve[model_id == "B_back2_unadj"]
  expect_equal(nrow(cv), fx$cfg$predict$gap_curve_max); expect_true(all(cv$lo < cv$est & cv$est < cv$hi))
  b <- r$model_B_post_gap_activity[model_id == "B_back2_unadj"]                      # the curve is the fitted line
  expect_equal(cv$est[cv$gap_length == 2] - cv$est[cv$gap_length == 1], b$estimate, tolerance = 1e-6)
  expect_true(all(r$model_B_observed$n_gaps > 0))
  ba <- r$model_B_all_terms                                                          # every fixed effect (no intercept) of each specification
  expect_false("(Intercept)" %in% ba$term)
  expect_setequal(ba[model_id == "B_back2_unadj", term], "log2_gap")
  expect_true(all(c("log2_gap", "w_pre_day_steps") %in% ba[model_id == "B_back2_adj_pre", term]))
  expect_equal(ba[model_id == "B_back2_unadj" & term == "log2_gap", estimate], b$estimate)    # same estimate as the one-row table
  expect_true(all(ba$conf_low < ba$estimate & ba$estimate < ba$conf_high))
  ca <- r$model_C_all_terms
  if (nrow(ca) > 0L) {
    expect_false(any(grepl("\\|", ca$term)))                                         # thresholds are not predictors
    cc1 <- r$model_C_gap_length[model_id == "C_unadj" & status %in% c("ok", "ok_bootstrap_se")]
    if (nrow(cc1)) expect_equal(ca[model_id == "C_unadj" & term == "w_pre_day_steps", or], cc1$or)
  }
  pr <- r$model_C_predicted_probs
  if (nrow(pr) > 0L) {
    sums <- pr[, .(s = sum(prob)), by = .(model_id, within_z)]
    expect_equal(sums$s, rep(1, nrow(sums)), tolerance = 1e-8)                       # probabilities are coherent
    expect_true(all(pr$prob_low <= pr$prob & pr$prob <= pr$prob_high))
  }
})

test_that("models B and C are null when gap length is unrelated to activity (negative control)", {
  fx <- gap_fixture(seed = 4, n = c(Insomnia = 120, Diabetes = 120, Healthy = 240, Prostate = 60)); g <- gap_data_for(fx, "sleep")
  b <- suppressWarnings(fit_gap_dose_response(g, fx$cfg, "back2_dev", adjusted = FALSE, id = "x"))
  expect_true(b$coef$status %in% c("ok", "not_converged")); expect_lt(abs(b$coef$estimate), 0.15)
  c1 <- suppressWarnings(suppressMessages(fit_gap_length_ordinal(g, fx$cfg, FALSE, "y")))
  expect_true(c1$coef$status %in% c("ok", "ok_bootstrap_se", "po_rejected"))
  if (c1$coef$status %in% c("ok", "ok_bootstrap_se")) expect_lt(abs(log(c1$coef$or)), 0.2)
})

test_that("planted association between prior activity and gap length is recovered by model C's sign", {
  set.seed(3); n <- 600
  g <- data.table(person_id = rep(1:150, each = 4), w_pre_day_steps = rnorm(n))
  g[, lat := -0.8 * w_pre_day_steps + rnorm(.N)]                               # lower activity -> longer gap
  g[, gap_length := as.integer(cut(lat, c(-Inf, -0.5, 0.3, 1, Inf), labels = c(1, 2, 4, 9)))]
  g[, gap_length := c(1L, 2L, 4L, 9L)[gap_length]]
  g[, `:=`(length_bin = collapse_sparse_bins(gap_length_bin(gap_length, cfg_test()), 20), log2_gap = log2(gap_length),
           back1_dev = 0, back2_dev = 0, cohort = "A", age10 = 0, sex2 = factor("Female", levels = c("Female", "Male")))]
  r <- suppressWarnings(suppressMessages(fit_gap_length_ordinal(g, cfg_test(), FALSE, "z")))
  expect_true(r$coef$status %in% c("ok", "ok_bootstrap_se")); expect_lt(r$coef$or, 1 / 1.3)   # OR of a longer gap per +1 SD HIGHER: lower activity -> longer gaps, so below 1
  expect_equal(r$coef$or, exp(r$coef$estimate)); expect_equal(r$coef_all[term == "w_pre_day_steps", or], r$coef$or)
  expect_equal(nrow(r$pred), 3L * nlevels(droplevels(g$length_bin)))
})

test_that("too few gaps -> not_estimable; a failed fit is never reported as a result", {
  cfg <- cfg_test(); cfg$gap_models$min_gaps <- 1e6
  g <- data.table(person_id = 1:3, log2_gap = 1, back2_dev = 0, back1_dev = 0, w_pre_day_steps = 0, gap_length = 2L,
                  length_bin = factor("1", ordered = TRUE), cohort = "A", age10 = 0, sex2 = NA)
  expect_equal(fit_gap_dose_response(g, cfg)$coef$status, "not_estimable")
  expect_equal(fit_gap_length_ordinal(g, cfg)$coef$status, "not_estimable")
  expect_true(value_failed(list(info = data.table(status = "failed: boom"))))
  expect_false(value_failed(list(info = data.table(status = "not_estimable"))))
})

test_that("cluster-clustered mean interval is wider than the naive one when gaps cluster within people", {
  set.seed(2); id <- rep(1:40, each = 10); y <- rnorm(40)[id] + rnorm(400, sd = 0.3)
  ci <- cluster_mean_ci(y, id); naive <- 1.96 * sd(y) / sqrt(length(y))
  expect_gt(ci$hi - ci$mean, naive); expect_equal(ci$n_persons, 40L)
})

test_that("sparse random-effect matrices give the same ordinal mixed model as the dense default", {
  fx <- gap_fixture(seed = 6, n = c(Insomnia = 60, Diabetes = 60, Healthy = 120, Prostate = 30))
  g <- gap_data_for(fx, "sleep")
  d <- g[!is.na(w_pre_day_steps) & !is.na(length_bin), .(person_id, length_bin, w_pre_day_steps)]; d[, length_bin := droplevels(length_bin)]
  f <- length_bin ~ w_pre_day_steps + (1 | person_id)
  dense  <- ordinal::clmm(f, data = d, control = ordinal::clmm.control(useMatrix = FALSE))
  sparse <- ordinal::clmm(f, data = d, control = ordinal::clmm.control(useMatrix = TRUE))
  expect_equal(unname(sparse$beta), unname(dense$beta), tolerance = 1e-3)
  expect_equal(as.numeric(sparse$logLik), as.numeric(dense$logLik), tolerance = 1e-5)
  expect_equal(unname(sparse$alpha), unname(dense$alpha), tolerance = 1e-3)
})

test_that("the full-period sensitivity units of models B and C are fitted from their own gap table, only when it is provided, and are reused", {
  fx <- gap_fixture(); ck <- tempfile(); dk <- data_key(fx$store, "sleep")
  built <- c(primary = 0, full = 0)
  gd  <- lazy(function() { built[["primary"]] <<- built[["primary"]] + 1; gap_data_for(fx, "sleep") })
  gdf <- lazy(function() { built[["full"]] <<- built[["full"]] + 1; build_gap_level_data(fx$inp$day_tables$sleep, fx$inp$day_tables$activity, fx$inp$persons, fx$cfg, standardization = "full_period") })
  st <- suppressWarnings(suppressMessages(run_gap_stage("sleep", gd, fx$cfg, ck, dk, gdf)))
  expect_setequal(names(st), names(gap_unit_specs())); expect_equal(built, c(primary = 1, full = 1))
  r <- collect_gap_results(ck, fx$cfg, list(sleep = dk))
  expect_true("B_back2_adj_fullperiod" %in% r$model_B_post_gap_activity$model_id)
  b1 <- r$model_B_post_gap_activity[model_id == "B_back2_adj"]; b2 <- r$model_B_post_gap_activity[model_id == "B_back2_adj_fullperiod"]
  expect_false(isTRUE(all.equal(b1$estimate, b2$estimate)))                      # a different standardization gives a different estimate
  expect_true(nrow(r$model_B_observed_trajectory) > 0L)
  expect_true(all(c("model", "model_id", "term", "df", "gvif", "vif_equiv", "flag") %in% names(r$gap_collinearity)))
  expect_true(all(c("B", "C") %in% r$gap_collinearity$model) || nrow(r$model_C_gap_length[status %in% c("ok", "ok_bootstrap_se")]) == 0L)
  expect_true(nrow(r$gap_term_tests) > 0L)
  gd2 <- lazy(function() stop("not needed")); gdf2 <- lazy(function() stop("not needed"))
  expect_true(all(suppressMessages(run_gap_stage("sleep", gd2, fx$cfg, ck, dk, gdf2)) == "resumed"))
})

# ---- MVPA before the gap as a covariate of models B and C ----------------------------------------------------------------

test_that("MVPA is aligned with the day before the gap: the MVPA z-score of the last observed day, against the usual BEFORE that day (never a day after the gap)", {
  fx <- gap_fixture(); g <- gap_data_for(fx, "activity"); act <- fx$inp$day_tables$activity; gr <- get_gap_runs(act)
  pt <- predictor_table(act, c("act_steps", "act_mvpa"), fx$cfg, trailing = FALSE, keep_x = TRUE)
  setkey(pt, person_id, date)
  i <- which(!is.na(g$w_pre_day_mvpa))[c(1, 50, 120)]
  for (j in i) {
    row <- gr[person_id == g$person_id[j] & gap_length == g$gap_length[j]][1]
    pre  <- pt[.(row$person_id, row$prev_obs_date)]; post <- pt[.(row$person_id, row$next_obs_date)]
    expect_equal(g$w_pre_day_mvpa[j], pre$w_act_mvpa)                                   # MVPA of the day before the gap
    expect_false(isTRUE(all.equal(g$w_pre_day_mvpa[j], post$w_act_mvpa)))               # not the day the person comes back
    # by hand: the day's log1p(MVPA) minus the mean of the person's EARLIER days, in units of the pooled within-person SD
    x <- log1p(act[person_id == row$person_id & date <= row$prev_obs_date, mvpa_minutes]); n <- length(x)
    sdw <- attr(pt, "scales")[predictor == "act_mvpa", sd_within]
    expect_equal(g$w_pre_day_mvpa[j], (x[n] - mean(x[-n])) / sdw)
  }
  expect_equal(g$w_pre_day_mvpa[is.na(g$w_pre_day_mvpa)] * 0, g$w_pre_day_mvpa[is.na(g$w_pre_day_mvpa)])        # NA, never imputed
  expect_true(all(is.na(g$w_pre_day_mvpa[is.na(g$w_pre_day_steps)])))                   # no baseline for steps -> none for MVPA either (same days)
  expect_gt(cor(g$w_pre_day_steps, g$w_pre_day_mvpa, use = "complete.obs"), 0.3)        # steps and MVPA of the same day are correlated
})

test_that("NO LOOK-AHEAD for MVPA: changing MVPA on or after the gap leaves the pre-gap MVPA deviation unchanged (only the common scale moves)", {
  cfg <- cfg_test(); set.seed(5)
  days <- data.table(person_id = rep(1:20, each = 80), day = rep(1:80, 20))[!(day %in% 40:42)]
  days[, `:=`(date = as.Date("2024-01-01") + day - 1L, steps = round(runif(.N, 3000, 9000)), sedentary_minutes = 800, mvpa_minutes = round(runif(.N, 0, 60)))]
  setkey(days, person_id, date)
  persons <- data.table(person_id = 1:20, cohort = rep(c("Healthy", "Diabetes"), 10), sex = rep(c("Female", "Male"), each = 10), age_first = 50, cci_score = rep(0:3, 5))
  gaps_of <- function(d) build_gap_level_data(d, d, persons, cfg)[order(person_id)]
  g1 <- gaps_of(days)
  after <- copy(days); after[person_id == 7L & day >= 43L, mvpa_minutes := mvpa_minutes * 4]        # from the return day on
  before <- copy(days); before[person_id == 7L & day %in% 10:30, mvpa_minutes := 0]                 # earlier days: a long quiet spell before the gap
  r <- function(g2) g1[person_id == 7L, w_pre_day_mvpa] / g2[person_id == 7L, w_pre_day_mvpa]
  expect_lt(abs(r(gaps_of(after)) - 1), 0.05)                       # the scale constant moves slightly (one person's later days), the deviation itself does not
  expect_gt(abs(r(gaps_of(before)) - 1), 0.1)                       # earlier MVPA does move it
  # the steps columns of the gap table do not depend on the MVPA column at all
  noM <- copy(days)[, mvpa_minutes := NULL]
  g0 <- gaps_of(noM)
  expect_true(all(is.na(g0$w_pre_day_mvpa)))
  expect_equal(g0[, .(person_id, gap_length, w_pre_day_steps, back2_dev, post_k3)], g1[, .(person_id, gap_length, w_pre_day_steps, back2_dev, post_k3)])
})

test_that("adding MVPA changes nothing in the existing specifications: same rows, same fit, new units only", {
  fx <- gap_fixture(); g <- gap_data_for(fx, "sleep")
  b0 <- fit_gap_dose_response(g, fx$cfg, "back2_dev", TRUE, FALSE, "B_back2_adj")
  b1 <- fit_gap_dose_response(g[, !"w_pre_day_mvpa"], fx$cfg, "back2_dev", TRUE, FALSE, "B_back2_adj")          # the same fit on a table without the MVPA column
  expect_equal(b0$coef$estimate, b1$coef$estimate); expect_equal(b0$coef$n_gaps, b1$coef$n_gaps); expect_false(b0$coef$adjust_mvpa)
  expect_equal(b0$coef$n_gaps_no_mvpa, 0L)
  expect_equal(names(gap_unit_specs())[!grepl("mvpa", names(gap_unit_specs()))],
               c("B_back2_unadj", "B_back2_adj", "B_back2_adj_pre", "B_back1_adj", "B_back2_adj_fullperiod", "B_observed_trajectory", "C_unadj", "C_adj", "C_adj_fullperiod"))
  expect_null(gap_unit_specs()$B_back2_adj$adjust_mvpa)                 # the old specifications are untouched (so their checkpoints stay valid)
  expect_equal(gap_unit_specs()$C_adj_mvpa, list(model = "C", adjusted = TRUE, adjust_mvpa = TRUE))
})

test_that("the MVPA coefficient of model B is recovered, and omitting MVPA biases the steps coefficient when MVPA is the true driver", {
  set.seed(11); n <- 1500
  g <- data.table(person_id = rep(1:500, each = 3), gap_length = sample(1:8, n, TRUE))
  g[, w_pre_day_mvpa := rnorm(.N)]; g[, w_pre_day_steps := 0.7 * w_pre_day_mvpa + sqrt(1 - 0.49) * rnorm(.N)]      # correlated like steps and MVPA of one day
  g[, back2_dev := 0.30 * w_pre_day_mvpa + 0.00 * w_pre_day_steps + rnorm(.N, sd = 0.9)]                          # the truth: MVPA matters, steps do not
  g[, `:=`(back1_dev = back2_dev, log2_gap = log2(gap_length), cohort = "A", age10 = 0, sex2 = factor("Female", levels = c("Female", "Male")))]
  cfg <- cfg_test()
  without <- fit_gap_dose_response(g, cfg, "back2_dev", FALSE, TRUE, "p")
  withm   <- fit_gap_dose_response(g, cfg, "back2_dev", FALSE, TRUE, "pm", adjust_mvpa = TRUE)
  st0 <- without$coef_all[term == "w_pre_day_steps"]; st1 <- withm$coef_all[term == "w_pre_day_steps"]; mv <- withm$coef_all[term == "w_pre_day_mvpa"]
  expect_gt(st0$estimate, 0.12); expect_lt(st0$conf_low, st0$estimate)                           # steps pick up the MVPA effect when MVPA is left out
  expect_lt(abs(st1$estimate), 0.1)                                                                # ... and not when it is in the model
  expect_true(mv$conf_low < 0.30 && 0.30 < mv$conf_high)                                           # the planted MVPA effect is inside its 95% CI
  expect_gt(withm$gvif[term == "w_pre_day_mvpa", vif_equiv], 1.5); expect_lt(withm$gvif[term == "w_pre_day_mvpa", vif_equiv], 3)   # VIF of r = 0.7 is about 2
  expect_true(withm$coef$adjust_mvpa)
  g2 <- copy(g); g2[1:30, w_pre_day_mvpa := NA]
  expect_equal(fit_gap_dose_response(g2, cfg, "back2_dev", FALSE, TRUE, "pm", adjust_mvpa = TRUE)$coef$n_gaps_no_mvpa, 30L)   # gaps without MVPA are dropped and counted
})

test_that("model C with MVPA: the covariate enters the proportional-odds check, the model and the predicted probabilities; its odds ratio is exp(coefficient)", {
  set.seed(3); n <- 1200
  g <- data.table(person_id = rep(1:300, each = 4)); g[, w_pre_day_mvpa := rnorm(.N)]; g[, w_pre_day_steps := 0.5 * w_pre_day_mvpa + sqrt(0.75) * rnorm(.N)]
  g[, lat := -0.7 * w_pre_day_mvpa + rnorm(.N)]                                                  # higher MVPA -> shorter gaps; steps have no effect of their own
  g[, gap_length := c(1L, 2L, 4L, 9L)[as.integer(cut(lat, c(-Inf, -0.6, 0.2, 1, Inf), labels = FALSE))]]
  cfg <- cfg_test()
  g[, `:=`(length_bin = collapse_sparse_bins(gap_length_bin(gap_length, cfg), 20), log2_gap = log2(gap_length), back1_dev = 0, back2_dev = 0,
           cohort = sample(c("A", "B"), .N, TRUE), age10 = rnorm(.N), sex2 = factor(sample(c("Female", "Male"), .N, TRUE), levels = c("Female", "Male")))]
  r <- suppressWarnings(suppressMessages(fit_gap_length_ordinal(g, cfg, TRUE, "C_adj_mvpa", adjust_mvpa = TRUE)))
  expect_true(r$coef$status %in% c("ok", "ok_bootstrap_se")); expect_true(r$coef$adjust_mvpa)
  m <- r$coef_all[term == "w_pre_day_mvpa"]
  expect_equal(m$or, exp(m$estimate)); expect_lt(m$or_high, 1)                                     # OR of a longer gap per +1 SD MVPA: clearly below 1
  expect_equal(m$or_low, exp(m$conf_low)); expect_equal(m$or_high, exp(m$conf_high))
  s0 <- suppressWarnings(suppressMessages(fit_gap_length_ordinal(g, cfg, TRUE, "C_adj")))
  expect_lt(s0$coef_all[term == "w_pre_day_steps", or], r$coef_all[term == "w_pre_day_steps", or] - 0.1)                  # steps look protective only when MVPA is left out
  expect_equal(nrow(r$pred), 3L * nlevels(droplevels(g$length_bin)))
  expect_true("w_pre_day_mvpa" %in% r$gvif$term)
})

test_that("gap_mvpa_effects lists the MVPA term and the main terms with and without MVPA, with intervals, sample sizes and VIF", {
  fx <- gap_fixture(); ck <- tempfile(); dk <- data_key(fx$store, "sleep")
  suppressWarnings(suppressMessages(run_gap_stage("sleep", lazy(function() gap_data_for(fx, "sleep")), fx$cfg, ck, dk)))
  r <- collect_gap_results(ck, fx$cfg, list(sleep = dk)); m <- r$gap_mvpa_effects
  expect_true(all(c("model", "stream", "pair", "model_id", "mvpa_adjusted", "term", "effect_measure", "estimate", "conf_low", "conf_high", "p", "n_gaps", "n_persons", "n_gaps_no_mvpa", "vif_equiv") %in% names(m)))
  expect_true(all(m[term == "w_pre_day_mvpa", mvpa_adjusted]))                                    # MVPA only appears in the MVPA-adjusted units
  expect_setequal(m[model == "B", unique(pair)], c("B_back2_adj_mvpa", "B_back2_adj_pre_mvpa"))
  expect_true(all(m[model == "B" & term == "log2_gap", .N, by = .(pair)]$N == 2L))                  # the gap-length term with and without MVPA
  expect_true(all(m[term == "w_pre_day_mvpa", conf_low < estimate & estimate < conf_high]))
  expect_true(all(m[model == "B" & term == "w_pre_day_mvpa", !is.na(vif_equiv)]))
  expect_true(all(m$n_gaps > 0))
  if (any(m$model == "C")) expect_equal(m[model == "C" & !is.na(or), or], exp(m[model == "C" & !is.na(or), estimate]))
})
