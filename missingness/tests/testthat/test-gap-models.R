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
  a <- predictor_table(fx$inp$day_tables$activity, "act_steps", fx$cfg, trailing = FALSE, keep_x = FALSE)
  gr <- get_gap_runs(fx$inp$day_tables$sleep)
  i <- which(!is.na(g$w_pre_day_steps) & !is.na(g$back2_dev))[1]
  key <- g[i, .(person_id, gap_length)]
  row <- gr[person_id == key$person_id & gap_length == key$gap_length][1]
  expect_equal(g$w_pre_day_steps[i], a[person_id == row$person_id & date == row$prev_obs_date, w_act_steps])
  expect_equal(g$back2_dev[i], a[person_id == row$person_id & date == row$next_obs_date + 1L, w_act_steps])
  expect_equal(g$log2_gap[i], log2(g$gap_length[i]))
})

test_that("gap data from the activity-only predictor equal those from the full bundle's activity predictors", {
  fx <- gap_fixture(); b <- suppressMessages(prepare_stream_bundle("sleep", fx$inp$day_tables, fx$inp$persons, fx$cfg))
  g <- gap_data_for(fx, "sleep"); gr <- get_gap_runs(fx$inp$day_tables$sleep)
  a <- b$rs[, .(person_id, date, w = w_act_steps)]; setkey(a, person_id, date)
  expect_equal(sort(g$w_pre_day_steps), sort(a[.(gr$person_id, gr$prev_obs_date), w]))
  full <- predictor_table(fx$inp$day_tables$activity, "act_steps", fx$cfg, trailing = FALSE, keep_x = FALSE)   # rs lacks each person's last day
  setkey(full, person_id, date)
  expect_equal(sort(g$back2_dev), sort(full[.(gr$person_id, gr$next_obs_date + 1L), w_act_steps]))
})

test_that("each model B / C is its own checkpointed unit; resuming skips finished ones; one failure leaves the others", {
  fx <- gap_fixture(); ck <- tempfile(); dk <- data_key(fx$store, "sleep")
  built <- 0; gd <- lazy(function() { built <<- built + 1; gap_data_for(fx, "sleep") })
  st1 <- suppressWarnings(suppressMessages(run_gap_stage("sleep", gd, fx$cfg, ck, dk)))
  expect_setequal(names(st1), names(gap_unit_specs())); expect_true(all(st1 == "computed")); expect_equal(built, 1)
  files <- list.files(ck); expect_length(files, 6L)

  f <- file.path(ck, files[grepl("^gapC__sleep__C_adj", files)]); unlink(f)           # one unit lost
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
  expect_equal(nrow(r$model_B_post_gap_activity), 4L); expect_equal(nrow(r$model_C_gap_length), 0L)
  expect_setequal(r$gap_failures$id, c("C_unadj", "C_adj")); expect_true(all(nzchar(r$gap_failures$error)))

  st2 <- suppressWarnings(suppressMessages(run_gap_stage("sleep", lazy(function() gap_data_for(fx, "sleep")), fx$cfg, ck, dk)))
  expect_true(all(st2[grepl("^B_", names(st2))] == "resumed")); expect_true(all(st2[grepl("^C_", names(st2))] == "computed"))
  r2 <- collect_gap_results(ck, fx$cfg, list(sleep = dk))
  expect_equal(nrow(r2$gap_failures), 0L); expect_equal(nrow(r2$model_C_gap_length), 2L)
})

test_that("results collect into tables, curves and predicted probabilities that are coherent", {
  fx <- gap_fixture(); ck <- tempfile(); dk <- data_key(fx$store, "sleep")
  suppressWarnings(suppressMessages(run_gap_stage("sleep", lazy(function() gap_data_for(fx, "sleep")), fx$cfg, ck, dk)))
  r <- collect_gap_results(ck, fx$cfg, list(sleep = dk))
  expect_setequal(r$model_B_post_gap_activity$model_id, c("B_back2_unadj", "B_back2_adj", "B_back2_adj_pre", "B_back1_adj"))
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
  expect_true(r$coef$status %in% c("ok", "ok_bootstrap_se")); expect_gt(r$coef$or, 1.3)       # OR per 1 SD BELOW usual
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
