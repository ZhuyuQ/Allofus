mock_run <- function(cfg, fp, existing = character(), fail_on = NULL, ...) {
  calls <- list()
  executor <- function(sql, path) {
    if (!is.null(fail_on) && grepl(fail_on, path)) stop("simulated BigQuery error")
    calls[[length(calls) + 1L]] <<- list(sql = sql, path = path)
    invisible(path)
  }
  lister <- function(p) if (any(grepl(p, existing, fixed = TRUE))) "x.csv" else if (length(calls) && any(vapply(calls, function(z) z$path == p, logical(1)))) "x.csv" else stop("none")
  res <- run_bq_exports(cfg, fp, executor = executor, lister = lister, ...)
  list(res = res, calls = calls)
}

test_that("16 exports are issued, each with its own cohort's ids and the path the ETL reads", {
  cfg <- make_config(list(export = list(bucket = "gs://b", owner = "me@x.org", date = "20261006")))
  fp <- data.table(person_id = 1:8, cohort = rep(c("Diabetes", "Insomnia", "Prostate", "Healthy"), each = 2))
  r <- suppressMessages(capture.output(out <- mock_run(cfg, fp)))
  expect_length(out$calls, 16L); expect_equal(sum(out$res$status == "exported"), 16L)
  paths <- vapply(out$calls, `[[`, "", "path")
  expect_setequal(paths, as.vector(outer(cfg$cohorts, c("activity", "sleep", "hr", "meta"),
                                         Vectorize(function(co, s) export_path(cfg, co, s)))))
  d <- out$calls[[which(paths == export_path(cfg, "Diabetes", "activity"))]]$sql
  expect_match(d, "IN \\(1, 2\\)"); expect_false(grepl("\\b3\\b", d))
  expect_match(d, "`activity_summary` activity_summary")
  expect_match(out$calls[[which(paths == export_path(cfg, "Prostate", "meta"))]]$sql, "IN \\(5, 6\\)")
})

test_that("SQL asks for every raw column the ETL needs", {
  need <- list(activity = c("person_id", "date", "steps", "calories_out", "sedentary_minutes", "lightly_active_minutes",
                            "fairly_active_minutes", "very_active_minutes"),
               sleep = c("person_id", "sleep_date", "is_main_sleep", "minute_in_bed", "minute_asleep"),
               hr = c("person_id", "date", "zone_name", "min_heart_rate", "max_heart_rate", "minute_in_zone", "calorie_count"),
               meta = c("person_id", "date_of_birth", "race", "sex_at_birth", "ethnicity", "self_reported_category"))
  for (s in names(need)) {
    sql <- bq_stream_sql(s, c(10, 20))
    for (cl in need[[s]]) expect_true(grepl(cl, sql, fixed = TRUE), info = paste(s, cl))
  }
  expect_error(bq_stream_sql("nope", 1), "unknown stream")
})

test_that("existing exports are skipped unless overwrite = TRUE; a failure is reported and stops at the end", {
  cfg <- make_config(list(export = list(bucket = "gs://b", owner = "me", date = "20261006")))
  fp <- data.table(person_id = 1:4, cohort = c("Diabetes", "Diabetes", "Healthy", "Healthy"))
  have <- export_path(cfg, "Diabetes", "hr")
  out <- suppressMessages(capture.output(o <- mock_run(cfg, fp, existing = have, cohorts = c("Diabetes", "Healthy"))))
  expect_equal(o$res[cohort == "Diabetes" & stream == "hr", status], "skipped")
  expect_length(o$calls, 7L)
  invisible(suppressMessages(capture.output(o2 <- mock_run(cfg, fp, existing = have, overwrite = TRUE, cohorts = "Diabetes"))))
  expect_length(o2$calls, 4L)                       # overwrite = TRUE redoes all four Diabetes exports
  expect_error(suppressMessages(capture.output(mock_run(cfg, fp, fail_on = "healthy_sleep", cohorts = c("Diabetes", "Healthy")))),
               "1 export\\(s\\) failed.*simulated BigQuery error")
})

test_that("the default executor says what is wrong with the CDR dataset or billing project", {
  cfg <- make_config(list(export = list(cdr = "", billing = "my-project")))
  expect_error(bq_executor("select 1", "gs://x/y.csv", cfg), "must look like .project.dataset.*WORKSPACE_CDR")
  cfg$export$cdr <- "a.b.c"
  expect_error(bq_executor("select 1", "gs://x/y.csv", cfg), "must look like")
  cfg2 <- make_config(list(export = list(cdr = "proj.dataset", billing = "")))
  expect_error(bq_executor("select 1", "gs://x/y.csv", cfg2), "no billing project")
})

test_that("a failed export puts the real reason in the error message", {
  cfg <- make_config(list(export = list(bucket = "gs://b", owner = "me", date = "20261006", cdr = "", billing = "p")))
  fp <- data.table(person_id = 1:2, cohort = "Prostate")
  expect_error(suppressMessages(capture.output(run_bq_exports(cfg, fp, cohorts = "Prostate", lister = function(p) stop("none")))),
               "4 export\\(s\\) failed. Reason.*project.dataset")
})

test_that("the optional device export asks for device_version and src_id, uses its own path, and is not part of the default run", {
  sql <- bq_stream_sql("device", c(10, 20))
  for (cl in c("person_id", "device_date", "device_version", "device_type")) expect_true(grepl(cl, sql, fixed = TRUE), info = cl)
  expect_false(grepl("src_id", sql, fixed = TRUE))
  expect_match(sql, "`device` device"); expect_match(sql, "IN \\(10, 20\\)")
  cfg <- make_config(list(export = list(bucket = "gs://b", owner = "me@x.org", date = "20261006")))
  expect_match(export_path(cfg, "Prostate", "device"), "pca_device/pca_device_\\*\\.csv$")
  fp <- data.table(person_id = 1:8, cohort = rep(c("Diabetes", "Insomnia", "Prostate", "Healthy"), each = 2))
  invisible(suppressMessages(capture.output(out <- mock_run(cfg, fp, streams = "device"))))
  expect_length(out$calls, 4L); expect_true(all(grepl("_device_\\*\\.csv$", vapply(out$calls, `[[`, "", "path"))))
})

test_that("the device probe reports which wanted columns exist and writes nothing", {
  cfg <- make_config()
  ok <- bq_probe_device_table(cfg, sampler = function(sql) { expect_match(sql, "FROM `device` LIMIT 5"); data.frame(person_id = 1, device_date = "2022-01-01", device_version = "Charge 4", device_type = "TRACKER", src_id = "x", battery = "High") })
  expect_true(ok$all_found); expect_true(all(ok$columns_check$in_table))
  bad <- bq_probe_device_table(cfg, sampler = function(sql) data.frame(person_id = 1, device_date = "d", device_type = "TRACKER"))
  expect_false(bad$all_found); expect_equal(bad$columns_check[in_table == FALSE, column], "device_version")
})

test_that("device coverage is asked per cohort with that cohort's ids and reported as percentages", {
  cfg <- make_config()
  fp <- data.table(person_id = 1:6, cohort = rep(c("Diabetes", "Healthy"), each = 3))
  seen <- character()
  cov <- bq_device_coverage(cfg, fp, cohorts = c("Diabetes", "Healthy"), sampler = function(sql) {
    seen <<- c(seen, sql)
    data.frame(n_with_activity = 3, n_with_device_row = 3, n_with_version = 2, n_activity_with_version = 2, n_activity_days = 300, n_activity_days_with_version = 150)
  })
  expect_match(seen[1L], "IN \\(1, 2, 3\\)"); expect_match(seen[2L], "IN \\(4, 5, 6\\)")
  expect_equal(cov$cohort, c("Diabetes", "Healthy")); expect_equal(cov$n_people, c(3, 3))
  expect_equal(cov$pct_activity_with_version, c(66.7, 66.7)); expect_equal(cov$pct_days_with_version, c(50, 50))
  expect_match(seen[1L], "!= 'SCALE'", fixed = TRUE)                            # scales do not count as a wearable
  v <- bq_device_versions(cfg, fp, sampler = function(sql) { expect_match(sql, "GROUP BY device_type, device_version"); expect_match(sql, "IN \\(1, 2, 3, 4, 5, 6\\)")
    data.frame(device_type = c("TRACKER", "SCALE"), device_version = c("Charge 4", "Aria Air"), n_people = c(5, 3), n_rows = c(50, 9)) })
  expect_equal(v$device_version, c("Charge 4", "Aria Air"))
})
