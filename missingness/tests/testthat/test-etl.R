test_that("drop_uninformative_columns drops all-NA / constant columns and records why", {
  dt <- data.table(person_id = 1:10, date = Sys.Date() + 0:9, a = 1:10, allna = NA_real_, const = 5, mostly = c(1, rep(NA, 9)))
  out <- drop_uninformative_columns(dt, na_frac = 0.9)
  expect_setequal(names(out), c("person_id", "date", "a"))
  expect_setequal(attr(out, "dropped_columns")$column, c("allna", "const", "mostly"))
})

test_that("activity: duplicates collapse to the largest step count, mvpa is fairly + very", {
  raw <- data.table(person_id = c(1, 1, 1), date = c("2024-01-01", "2024-01-01", "2024-01-02"),
                    steps = c(100, 900, 500), fairly_active_minutes = c(1, 2, 3), very_active_minutes = c(0, 5, NA),
                    sedentary_minutes = c(10, 20, 30), elevation = NA_real_)
  out <- standardize_activity(raw)
  expect_equal(nrow(out), 2L); expect_equal(out$steps, c(900, 500))
  expect_equal(out$mvpa_minutes, c(7, 3)); expect_false("elevation" %in% names(out))
})

test_that("sleep: main + nap collapse to one day; efficiency uses day totals; zero-asleep gets NA", {
  raw <- data.table(person_id = 1, sleep_date = c("2024-01-01", "2024-01-01", "2024-01-02", "2024-01-03"),
                    is_main_sleep = c("true", "false", "true", "true"),
                    minute_in_bed = c(400, 50, 300, 300), minute_asleep = c(360, 40, 270, 0),
                    minute_awake = c(40, 10, 30, 300))
  out <- standardize_sleep(raw)
  expect_equal(nrow(out), 3L)
  d1 <- out[date == as.Date("2024-01-01")]
  expect_equal(d1$minute_in_bed, 450); expect_equal(d1$minute_asleep, 400)
  expect_equal(d1$sleep_efficiency, 400 / 450 * 100); expect_true(d1$has_nap && d1$has_main_sleep)
  expect_true(is.na(out[date == as.Date("2024-01-03")]$sleep_efficiency))
  expect_true(all(c("has_main_sleep", "has_nap", "n_records") %in% names(out)))   # protected from the drop rule
})

test_that("hr: four zone rows become one wide row; absent zone is 0; totals are sums", {
  raw <- data.table(person_id = 1, date = "2024-01-01",
                    zone_name = c("Out of Range", "Fat Burn", "Cardio"),         # Peak absent that day
                    min_heart_rate = c(30, 100, 130), max_heart_rate = c(100, 130, 160),
                    minute_in_zone = c(1300, 100, 30), calorie_count = c(2000, 300, 200))
  raw2 <- copy(raw)[, date := "2024-01-02"][, minute_in_zone := c(1200, 150, 40)]
  raw2 <- rbind(raw2, data.table(person_id = 1, date = "2024-01-02", zone_name = "Peak", min_heart_rate = 160,
                                 max_heart_rate = 220, minute_in_zone = 10, calorie_count = 90))
  out <- standardize_hr(rbind(raw, raw2))
  expect_equal(nrow(out), 2L)
  d1 <- out[date == as.Date("2024-01-01")]; d2 <- out[date == as.Date("2024-01-02")]
  expect_equal(d1$mins_peak, 0); expect_equal(d1$total_hr_min, 1430); expect_equal(d1$elevated_min, 130)
  expect_equal(d2$total_hr_min, 1400); expect_equal(d2$elevated_min, 200)
  expect_equal(d1$total_cal, 2500)
})

test_that("selection check warns when the upstream rule is violated", {
  cfg <- cfg_test()
  days <- data.table(person_id = 1, date = as.Date("2024-01-01") + c(0:5, 40:60))      # 34-day gap
  expect_warning(res <- check_selection_rule(days, "activity", cfg), "selection rule")
  expect_equal(res$observed_max_gap, 34L)
})

test_that("export paths follow the bucket layout, including the pca prefix", {
  cfg <- make_config(list(export = list(bucket = "gs://b", owner = "me@x.org", date = "20260511")))
  expect_equal(export_path(cfg, "Prostate", "activity"), "gs://b/bq_exports/me@x.org/20260511/pca_activity/pca_activity_*.csv")
  expect_equal(export_path(cfg, "Healthy", "sleep"), "gs://b/bq_exports/me@x.org/20260511/healthy_sleep/healthy_sleep.csv")
  expect_equal(export_path(cfg, "Diabetes", "meta"), "gs://b/bq_exports/me@x.org/20260511/diabetes_meta/diabetes_meta_*.csv")
})
