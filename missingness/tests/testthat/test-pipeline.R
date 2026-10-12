test_that("day tables round-trip through RDS storage", {
  cfg <- cfg_test(); inp <- simulate_inputs(cfg, n_per_cohort = c(Insomnia = 10, Diabetes = 10, Healthy = 10, Prostate = 10), seed = 1)
  d <- tempfile(); write_day_table(inp$day_tables$activity, d, "activity", "rds")
  back <- read_day_table(d, "activity", "rds")
  expect_equal(nrow(back), nrow(inp$day_tables$activity)); expect_s3_class(back$date, "Date")
})

# Write a simulated study as local CSV exports in the bucket layout (activity and meta in two
# shards, sleep and hr as single files). `membership` says which cohorts export which people.
write_sim_exports <- function(sim, cfg, tmp, membership) {
  by_cohort <- split(membership$person_id, membership$cohort)
  write_export <- function(dt, cohort, stream, sharded) {
    prefix <- cfg$cohort_prefix[[cohort]]; stem <- if (stream == "meta") "meta" else cfg$streams[[stream]]$stem
    folder <- file.path(tmp, "bq_exports", "me", "20260511", paste0(prefix, "_", stem)); dir.create(folder, recursive = TRUE)
    if (sharded) {
      half <- ceiling(nrow(dt) / 2)
      fwrite(dt[1:half], file.path(folder, paste0(prefix, "_", stem, "_000000000000.csv")))
      fwrite(dt[(half + 1):nrow(dt)], file.path(folder, paste0(prefix, "_", stem, "_000000000001.csv")))
    } else fwrite(dt, file.path(folder, paste0(prefix, "_", stem, ".csv")))
  }
  for (co in cfg$cohorts) {
    ids <- by_cohort[[co]]
    write_export(sim$raw$activity[person_id %in% ids], co, "activity", TRUE)
    write_export(sim$raw$sleep[person_id %in% ids], co, "sleep", FALSE)
    write_export(sim$raw$hr[person_id %in% ids], co, "hr", FALSE)
    write_export(sim$meta[person_id %in% ids], co, "meta", TRUE)
  }
  fwrite(membership, file.path(tmp, "final_pid.txt"), sep = "\t")
  fwrite(sim$cci, file.path(tmp, "cci.txt"), sep = "\t")
}

etl_cfg <- function(tmp, ...) make_config(c(list(boot = list(B = 30L), export = list(bucket = tmp, owner = "me", date = "20260511")), list(...)))
run_sim_etl <- function(cfg, tmp) {
  res <- NULL
  invisible(capture.output(suppressMessages(
    res <- run_etl(cfg, file.path(tmp, "L1"), lister = Sys.glob, reader = function(f, ...) fread(f, ...),
                   final_pid_path = file.path(tmp, "final_pid.txt"), cci_path = file.path(tmp, "cci.txt")))))
  res
}

test_that("run_etl reads sharded and single-file exports and matches building the tables directly", {
  tmp <- tempfile(); dir.create(tmp); cfg <- etl_cfg(tmp)
  sim <- simulate_study(n_per_cohort = c(Insomnia = 15, Diabetes = 15, Healthy = 30, Prostate = 10), seed = 21)
  write_sim_exports(sim, cfg, tmp, sim$final_pid)
  res <- run_sim_etl(cfg, tmp)
  loaded <- load_inputs(cfg, file.path(tmp, "L1"))                      # a day store: nothing big is read yet
  expect_s3_class(loaded$day_tables, "day_store")
  direct <- build_day_table(sim$raw$activity, "activity", cfg)
  act <- loaded$day_tables$get("activity")
  expect_equal(nrow(act), nrow(direct)); expect_equal(sum(act$steps), sum(direct$steps))
  expect_equal(nrow(loaded$day_tables$get("sleep")), nrow(build_day_table(sim$raw$sleep, "sleep", cfg)))
  expect_equal(nrow(res$persons), nrow(sim$final_pid)); expect_equal(nrow(loaded$persons), nrow(res$persons))
  expect_setequal(res$summary$stream, c("activity", "sleep", "hr"))     # the ETL returns summaries, not the tables
  expect_true(all(file.exists(file.path(tmp, "L1", c("manifest.tsv", "dropped_columns.tsv", "cohort_overlap.tsv")))))
})

test_that("people who belong to two cohorts are counted once: excluded or kept in one cohort", {
  tmp <- tempfile(); dir.create(tmp)
  sim <- simulate_study(n_per_cohort = c(Insomnia = 15, Diabetes = 15, Healthy = 30, Prostate = 10), seed = 21)
  ins <- sim$final_pid[cohort == "Insomnia", person_id][1]
  membership <- rbind(sim$final_pid, data.table(person_id = ins, cohort = "Diabetes"))   # one person in two cohorts
  cfg <- etl_cfg(tmp); write_sim_exports(sim, cfg, tmp, membership)

  r1 <- run_sim_etl(cfg, tmp)                                                  # default policy: exclude
  d1 <- load_inputs(cfg, file.path(tmp, "L1"))$day_tables$get("activity")
  expect_false(ins %in% r1$persons$person_id); expect_false(ins %in% d1$person_id)
  expect_equal(nrow(r1$persons), nrow(sim$final_pid) - 1L)
  expect_equal(anyDuplicated(d1, by = c("person_id", "date")), 0L)
  ov <- fread(file.path(tmp, "L1", "cohort_overlap.tsv")); expect_equal(ov$n_people, 1L)

  cfg2 <- etl_cfg(tmp, overlap = list(policy = "priority")); cfg2$overlap$priority <- c("Prostate", "Diabetes", "Insomnia", "Healthy")
  r2 <- run_sim_etl(cfg2, tmp)
  expect_equal(as.character(r2$persons[person_id == ins, cohort]), "Diabetes")   # Diabetes outranks Insomnia
  d2 <- load_inputs(cfg2, file.path(tmp, "L1"))$day_tables$get("activity")
  expect_true(ins %in% d2$person_id)
  expect_equal(anyDuplicated(d2, by = c("person_id", "date")), 0L)
  expect_equal(nrow(r2$persons), nrow(sim$final_pid))
})

test_that("pipeline writes tables, figures and metadata; primary effect recovered end to end", {
  cfg <- cfg_test()
  inp <- suppressMessages(simulate_inputs(cfg, n_per_cohort = c(Insomnia = 60, Diabetes = 60, Healthy = 120, Prostate = 30), seed = 5))
  out <- tempfile()
  res <- suppressWarnings(suppressMessages(run_missingness_pipeline(
    inp$day_tables, inp$persons, cfg, out, stages = c("missingness", "temporal", "descriptive", "onset"),
    onset_roles = "primary")))
  expect_true(all(file.exists(file.path(out, "tables", c("cohort_missingness.tsv", "transition.tsv", "weekday_tests.tsv",
                                                          "metric_summary.tsv", "onset_coefficients.tsv", "primary_test.tsv")))))
  expect_true(file.exists(file.path(out, "figures", "missing_pct_by_cohort.png")))
  expect_true(all(file.exists(file.path(out, c("run_config.json", "session_info.txt", "input_manifest.tsv")))))
  expect_equal(nrow(res$primary_test), 3L)
  expect_true(all(res$primary_test$or > 1 & res$primary_test$p_holm < 0.05))
})

test_that("the onset stage resumes from checkpoints, runs one stream at a time and collects everything finished", {
  cfg <- cfg_test()
  inp <- suppressMessages(simulate_inputs(cfg, n_per_cohort = c(Insomnia = 60, Diabetes = 60, Healthy = 120, Prostate = 30), seed = 5))
  out <- tempfile(); ck <- file.path(out, "checkpoints")
  quiet <- function(expr) suppressWarnings(suppressMessages(capture.output(r <- expr)))
  r1 <- NULL; invisible(quiet(r1 <- run_missingness_pipeline(inp$day_tables, inp$persons, cfg, out, stages = "onset",
                                                             onset_roles = "primary", model_streams = "activity")))
  expect_equal(r1$primary_test$stream, "activity")                         # only the requested stream so far
  expect_true(all(r1$unit_run_status$status == "computed"))
  f <- file.path(ck, list.files(ck, pattern = "^onset__activity__activity_joint"))[1]; t0 <- file.mtime(f)
  r2 <- NULL; invisible(quiet(r2 <- run_missingness_pipeline(inp$day_tables, inp$persons, cfg, out, stages = "onset",
                                                             onset_roles = "primary", model_streams = c("activity", "sleep"))))
  expect_equal(sort(r2$primary_test$stream), c("activity", "sleep"))      # sleep added, activity not refitted
  expect_true(all(r2$unit_run_status[grepl("activity", stream), status] == "resumed"))
  expect_true(all(r2$unit_run_status[grepl("sleep", stream), status] == "computed"))
  expect_equal(file.mtime(f), t0)
  expect_equal(r2$primary_test[stream == "activity", estimate], r1$primary_test[stream == "activity", estimate])
  expect_true(all(c("vif", "cohort_interaction", "cohort_specific", "risk_unadjusted") %in% names(r2)))
  expect_true(all(c("checkpoint_status", "unit_run_status") %in% names(r2))); expect_true(all(r2$checkpoint_status$state == "valid"))
  expect_true(file.exists(file.path(out, "results_summary.md"))); expect_true(file.exists(file.path(out, "figures", "figure_captions.md")))

  # settings that change the result -> the old checkpoints are not used, and are reported as stale
  cfg_other <- cfg; cfg_other$onset$days_since_cap <- 30L
  dks <- list(activity = data_key(as_day_store(inp$day_tables, inp$persons), "activity"))
  expect_equal(nrow(collect_onset_results(ck, cfg_other, dks)$onset_models), 0L)
  expect_true(all(checkpoint_report(ck, cfg_other, dks)$state == "stale"))
})

test_that("a crash is survivable: lost checkpoints are redone, the rest skipped; changed input data invalidates", {
  cfg <- cfg_test()
  inp <- suppressMessages(simulate_inputs(cfg, n_per_cohort = c(Insomnia = 60, Diabetes = 60, Healthy = 120, Prostate = 30), seed = 5))
  out <- tempfile(); ck <- file.path(out, "checkpoints")
  run <- function(tables = inp$day_tables, c_ = cfg) { r <- NULL; invisible(suppressWarnings(suppressMessages(capture.output(
    r <- run_missingness_pipeline(tables, inp$persons, c_, out, stages = "onset", onset_roles = "primary", model_streams = "activity"))))); r }
  r1 <- run(); n_units <- nrow(r1$unit_run_status)
  unlink(file.path(ck, list.files(ck, pattern = "cohort_Healthy")))                       # simulate a crash that lost one unit
  writeBin(as.raw(1:9), file.path(ck, list.files(ck, pattern = "cohort_interaction")))   # ... and damaged another
  r2 <- run(); st <- setNames(r2$unit_run_status$status, r2$unit_run_status$unit)
  expect_equal(sort(names(st)[st == "computed"]), sort(c("extra:cohort_Healthy", "extra:cohort_interaction")))
  expect_equal(sum(st == "resumed"), n_units - 2L)
  expect_equal(r2$primary_test$estimate, r1$primary_test$estimate)

  changed <- copy(inp$day_tables); changed$activity <- copy(changed$activity); changed$activity[1:50, steps := steps + 1]
  r3 <- run(changed); expect_true(all(r3$unit_run_status$status == "computed"))          # new data -> everything re-run

  tight <- cfg; tight$resources$memory_budget_gb <- 0                                    # every model refused by the memory guard
  inp2 <- changed
  r4 <- run(inp2, tight)
  expect_true(all(r4$unit_run_status$status == "resumed"))                               # results exist already: guard is not even consulted
})

test_that("a model that cannot run is reported as failed in the tables and never as a result", {
  cfg <- cfg_test(); cfg$resources$memory_budget_gb <- 0
  inp <- suppressMessages(simulate_inputs(cfg, n_per_cohort = c(Insomnia = 40, Diabetes = 40, Healthy = 80, Prostate = 20), seed = 5))
  out <- tempfile()
  res <- NULL; invisible(suppressWarnings(suppressMessages(capture.output(
    res <- run_missingness_pipeline(inp$day_tables, inp$persons, cfg, out, stages = "onset", onset_roles = "primary", model_streams = "activity")))))
  expect_true(any(res$unit_run_status$status == "failed"))
  expect_true(nrow(res$unit_failures) > 0L); expect_match(res$unit_failures$error[1], "Nothing was subsampled")
  expect_true(is.null(res$primary_test) || nrow(res$primary_test) == 0L)
  expect_true(file.exists(file.path(out, "tables", "unit_failures.tsv")))
  expect_match(paste(readLines(file.path(out, "results_summary.md")), collapse = "\n"), "FAILED units")
})

test_that("results accumulate across calls: gap-length models, sensitivity files and a different sample size", {
  cfg <- cfg_test()
  inp <- suppressMessages(simulate_inputs(cfg, n_per_cohort = c(Insomnia = 60, Diabetes = 60, Healthy = 120, Prostate = 30), seed = 5))
  out <- tempfile(); ck <- file.path(out, "checkpoints")
  quiet <- function(expr) invisible(suppressWarnings(suppressMessages(capture.output(r <- expr))))

  # gap-length models: one stream per call, the tables on disk must keep every stream
  for (s in c("activity", "sleep")) {
    r <- NULL; quiet(r <- run_missingness_pipeline(inp$day_tables, inp$persons, cfg, out, stages = "gap_models", model_streams = s))
  }
  expect_setequal(r$model_B_post_gap_activity$stream, c("activity", "sleep"))
  expect_setequal(fread(file.path(out, "tables", "model_B_post_gap_activity.tsv"))$stream, c("activity", "sleep"))
  expect_setequal(fread(file.path(out, "tables", "model_C_gap_length.tsv"))$stream, c("activity", "sleep"))

  # a different sample size (the speed knob for exploratory models) must not hide earlier results
  r1 <- NULL; quiet(r1 <- run_missingness_pipeline(inp$day_tables, inp$persons, cfg, out, stages = "onset", onset_roles = "primary", model_streams = "activity"))
  cfg_s <- cfg; cfg_s$model$exploratory_sample_persons <- 120L
  r2 <- NULL; quiet(r2 <- run_missingness_pipeline(inp$day_tables, inp$persons, cfg_s, out, stages = "onset", onset_roles = "exploratory", model_streams = "activity"))
  expect_true(all(c("primary", "exploratory") %in% r2$onset_models$role))
  expect_equal(r2$primary_test$estimate, r1$primary_test$estimate)
})

test_that("descriptive-only runs write no results summary, an empty failure table and no model checkpoints", {
  cfg <- cfg_test()
  inp <- suppressMessages(simulate_inputs(cfg, n_per_cohort = c(Insomnia = 20, Diabetes = 20, Healthy = 40, Prostate = 10), seed = 5))
  out <- tempfile(); res <- NULL
  invisible(suppressWarnings(suppressMessages(capture.output(
    res <- run_missingness_pipeline(inp$day_tables, inp$persons, cfg, out, stages = c("missingness", "temporal"))))))
  expect_false(file.exists(file.path(out, "results_summary.md")))        # no model ran: nothing to summarise
  expect_equal(nrow(res$unit_failures), 0L); expect_equal(nrow(res$checkpoint_status), 0L)
  expect_false(file.exists(file.path(out, "tables", "unit_failures.tsv")))
})

test_that("the 'missing_figures' stage redraws the missingness summaries only: no models, no checkpoints, other results untouched", {
  cfg <- cfg_test()
  inp <- suppressMessages(simulate_inputs(cfg, n_per_cohort = c(Insomnia = 40, Diabetes = 40, Healthy = 80, Prostate = 20), seed = 5))
  out <- tempfile(); dir.create(out); writeLines("sentinel", file.path(out, "results_summary.md"))
  res <- NULL
  invisible(suppressWarnings(suppressMessages(capture.output(
    res <- run_missingness_pipeline(inp$day_tables, inp$persons, cfg, out, stages = "missing_figures")))))
  expect_true(all(file.exists(file.path(out, "figures", paste0(c("missing_by_study_week", "missing_by_year", "missing_by_weekday", "missing_by_age",
         "missing_pct_by_sex", "missing_pct_by_age_group", "missing_pct_by_race_plot", "missing_pct_by_cci_group"), ".png")))))
  expect_true(all(file.exists(file.path(out, "tables", c("missing_pct_by_stratum.tsv", "missing_pct_by_stratum_tests.tsv", "missing_pct_by_stratum_models.tsv", "study_week_rates.tsv",
                                                          "missing_trend_tests.tsv", "missing_trend_heterogeneity.tsv", "missing_weekday_tests.tsv", "missing_weekday_odds_ratios.tsv")))))
  expect_false(dir.exists(file.path(out, "checkpoints")))                                   # no model was touched
  expect_equal(readLines(file.path(out, "results_summary.md")), "sentinel")                  # results summary not overwritten
  expect_false(any(c("transition", "weekday_tests", "metric_summary", "onset_coefficients") %in% names(res)))
  expect_match(paste(readLines(file.path(out, "figures", "figure_captions.md")), collapse = " "), "5-year moving window")
})
