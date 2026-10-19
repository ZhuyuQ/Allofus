cfg_dev <- function() { cfg <- cfg_test(); cfg$device$mc_B <- 500L; cfg$device$boot_B <- 100L; cfg }

dev_persons <- function() data.table(person_id = 1:8, cohort = factor(rep(c("Healthy", "Diabetes"), each = 4), levels = c("Insomnia", "Diabetes", "Healthy", "Prostate")),
                                     age_first = c(40, 50, 60, 45, 55, 65, 35, 52), sex = factor(c("Female", "Male")[c(1, 2, 1, 2, 1, 2, 1, 2)], levels = c("Female", "Male", "Other")))

test_that("a device version is grouped into a family; scales and empty versions are not wearables", {
  cfg <- cfg_dev()
  expect_equal(device_family(c("Charge 5", "Charge HR", "Versa 4", "Blaze", "Google Pixel Watch", NA), cfg),
               c("Charge", "Charge", "Versa", "Other wrist-worn", "Other wrist-worn", NA))
  expect_equal(device_form(c("Versa 4", "Blaze", "Alta HR", "Flex 2", "Charge 5", "One", "Zip", "MobileTrack", "Aria Air", "Weird Thing", NA), cfg),
               c("smartwatch", "smartwatch", "band", "band", "band", "clip-on", "clip-on", "phone", "unclassified", "unclassified", NA))
  raw <- data.table(person_id = c(1, 2, 3, 3), device_date = c("2022-01-02", "2022-01-03", "2022-01-04", NA), device_version = c("Charge 5", "Aria Air", "  ", "Versa 3"),
                    device_type = c("TRACKER", "SCALE", "TRACKER", "tracker"))
  n <- normalise_device_rows(raw, cfg)
  expect_equal(n[person_id == 2, is_scale], TRUE); expect_true(is.na(n[person_id == 3 & device_version != "Versa 3" | is.na(device_version), device_version][1L]))
  expect_equal(n[person_id == 3 & !is.na(device_version), is_scale], FALSE)       # type compared case-insensitively
})

test_that("person table: main version, status of everybody (also people without any row), timeline of dated wearable rows", {
  cfg <- cfg_dev(); persons <- dev_persons()
  raw <- rbind(
    data.table(person_id = 1, device_date = as.Date("2022-01-01") + c(0, 5, 9), device_version = c("Charge 4", "Charge 5", "Charge 5"), device_type = "TRACKER"),   # most records: Charge 5
    data.table(person_id = 2, device_date = as.Date("2022-02-01"), device_version = "Aria Air", device_type = "SCALE"),                                            # scale only
    data.table(person_id = 3, device_date = as.Date("2022-03-01"), device_version = NA_character_, device_type = "TRACKER"),                                        # row without a version
    data.table(person_id = 5, device_date = as.Date(NA), device_version = "Versa 4", device_type = "TRACKER"),                                                      # known, but no date
    data.table(person_id = 99, device_date = as.Date("2022-01-01"), device_version = "Sense", device_type = "TRACKER"))                                            # not in the person table
  t <- build_device_tables(raw, persons, cfg)
  p <- t$person
  expect_equal(nrow(p), 8L); expect_false(99 %in% p$person_id)
  expect_equal(p[person_id == 1, .(device_version_main, device_family, n_wearable_versions, status)],
               data.table(device_version_main = "Charge 5", device_family = "Charge", n_wearable_versions = 2L, status = DEVICE_STATUS[["known"]]))
  expect_equal(p[person_id == 2, status], DEVICE_STATUS[["scale_only"]]); expect_equal(p[person_id == 3, status], DEVICE_STATUS[["no_version"]])
  expect_equal(p[person_id == 4, status], DEVICE_STATUS[["no_record"]]); expect_equal(p[person_id == 5, status], DEVICE_STATUS[["known"]])
  expect_equal(p[, sum(device_missing)], 6L)                                          # 2, 3, 4, 6, 7, 8
  expect_equal(p[device_missing == TRUE, unique(device_family)], cfg$device$missing_label)
  expect_equal(nrow(t$timeline), 3L); expect_false(5 %in% t$timeline$person_id)          # person 5 has a version but no date
  # tie in the number of records: the most recent record wins
  r2 <- data.table(person_id = 1, device_date = as.Date(c("2022-01-01", "2022-06-01")), device_version = c("Charge 4", "Versa 3"), device_type = "TRACKER")
  expect_equal(build_device_tables(r2, persons, cfg)$person[person_id == 1, device_version_main], "Versa 3")
})

test_that("only WATCH data counts: phone tracking and clip-ons are not a watch version, and the wrist definition is configurable", {
  cfg <- cfg_dev(); persons <- dev_persons()
  raw <- rbind(data.table(person_id = 1, device_date = as.Date("2022-01-01") + 0:3, device_version = "MobileTrack", device_type = "TRACKER"),                   # phone only
               data.table(person_id = 2, device_date = as.Date("2022-01-01") + 0:1, device_version = "Zip", device_type = "TRACKER"),                           # clip-on only
               data.table(person_id = 3, device_date = as.Date("2022-01-01") + 0:5, device_version = "MobileTrack", device_type = "TRACKER"),                   # phone + a watch
               data.table(person_id = 3, device_date = as.Date("2022-02-01"), device_version = "Versa 3", device_type = "TRACKER"),                            # -> the watch is the main device
               data.table(person_id = 4, device_date = as.Date("2022-01-01"), device_version = "Alta HR", device_type = "TRACKER"))                           # band
  p <- build_device_tables(raw, persons, cfg)$person
  expect_equal(p[person_id %in% 1:2, status], rep(DEVICE_STATUS[["non_wrist"]], 2L)); expect_true(all(p[person_id %in% 1:2, device_missing]))
  expect_equal(p[person_id == 3, .(device_version_main, device_form_main, n_wearable_versions)], data.table(device_version_main = "Versa 3", device_form_main = "smartwatch", n_wearable_versions = 1L))
  expect_equal(p[person_id == 4, device_family], "Alta")
  tl <- build_device_tables(raw, persons, cfg)$timeline; expect_false(any(tl$device_version %in% c("MobileTrack", "Zip")))
  cfg2 <- cfg; cfg2$device$wrist_forms <- "smartwatch"                                      # smartwatches only: the band is no longer a watch
  expect_equal(build_device_tables(raw, persons, cfg2)$person[person_id == 4, status], DEVICE_STATUS[["non_wrist"]])
  cfg3 <- cfg; cfg3$device$wrist_forms <- c("smartwatch", "band", "clip-on")
  expect_equal(build_device_tables(raw, persons, cfg3)$person[person_id == 2, device_version_main], "Zip")
})

test_that("the device of a day is the most recent record on or before it; the days before the first record are NA unless back-filled", {
  tl <- data.table(person_id = c(1, 1, 2), device_date = as.Date(c("2022-01-10", "2022-03-01", "2022-01-05")),
                   device_version = c("Charge 4", "Charge 5", "Versa 3"), device_family = c("Charge", "Charge", "Versa"))
  days <- data.table(person_id = c(1, 1, 1, 1, 2, 3), date = as.Date(c("2022-01-09", "2022-01-10", "2022-02-15", "2022-03-01", "2022-06-01", "2022-01-01")))
  a <- device_at_dates(days, tl)
  expect_equal(a$device_version, c(NA, "Charge 4", "Charge 4", "Charge 5", "Versa 3", NA))
  expect_equal(a$person_id, days$person_id); expect_equal(a$date, days$date)
  b <- device_at_dates(days, tl, backfill = TRUE)
  expect_equal(b$device_version, c("Charge 4", "Charge 4", "Charge 4", "Charge 5", "Versa 3", NA))
})

test_that("Wilson, Newcombe and the bias-corrected Cramer's V behave", {
  w <- wilson_ci(0, 250); expect_equal(w$low, 0); expect_gt(w$high, 0.01); expect_lt(w$high, 0.02)
  w2 <- wilson_ci(50, 100); expect_equal(round(c(w2$low, w2$high), 3), c(0.404, 0.596))
  n <- newcombe_diff_ci(30, 100, 20, 100); expect_equal(n$diff, 0.1); expect_lt(n$low, 0.1); expect_gt(n$high, 0.1)
  expect_equal(cramers_v(matrix(c(50, 50, 50, 50), 2)), 0)                               # independence
  expect_gt(cramers_v(matrix(c(90, 10, 10, 90), 2)), 0.7)                                # strong association
  set.seed(1); ci <- cramers_v_ci(matrix(c(90, 10, 10, 90), 2), B = 200); expect_lt(ci[1L], cramers_v(matrix(c(90, 10, 10, 90), 2))); expect_gt(ci[2L], ci[1L])
  x <- chisq_table(matrix(c(30, 20, 10, 40, 25, 25), 3), B_mc = 500, B_boot = 50)
  expect_equal(x$df, 2); expect_true(x$p_monte_carlo > 0 && x$p_monte_carlo <= 1)
})

test_that("device_summary: missing share per cohort with intervals, the tests, the adjusted model (a cohort without events is not estimable) and the composition", {
  cfg <- cfg_dev()
  s <- small_study(-0.5); persons <- s$inp$persons
  raw <- sim_device_rows(persons, cfg)
  t <- build_device_tables(raw, persons, cfg)
  ds <- device_summary(t$person, persons, cfg)
  mb <- ds$device_missing_by_cohort
  expect_equal(mb[cohort == "All cohorts", n_people], nrow(persons)); expect_equal(mb[cohort != "All cohorts", sum(n_people)], nrow(persons))
  expect_true(all(mb$pct_low <= mb$pct_missing & mb$pct_missing <= mb$pct_high))
  expect_equal(mb[cohort == "Prostate", n_missing], 0L)                                          # simulated: nobody missing
  expect_equal(ds$device_missing_vs_reference[cohort == "Prostate", odds_ratio], 0)
  expect_true(all(ds$device_missing_vs_reference$p_holm >= ds$device_missing_vs_reference$p_fisher - 1e-12))
  m <- ds$device_missing_model
  expect_true(any(grepl("not estimable", m[term == "cohortProstate", note])))
  expect_true(all(c("cohortDiabetes", "cohortInsomnia", "age10") %in% m$term))
  comp <- ds$device_composition_by_cohort
  expect_equal(comp[, .(s = sum(pct)), by = cohort]$s, rep(100, uniqueN(comp$cohort)), tolerance = 1e-8)      # each cohort sums to 100%
  expect_equal(levels(comp$category)[1L], cfg$device$missing_label)
  tt <- ds$device_tests; expect_equal(nrow(tt), 3L); expect_true(all(tt$cramers_v >= 0 & tt$v_low <= tt$v_high, na.rm = TRUE))
  expect_equal(ds$device_composition_residuals[, sum(observed)], nrow(persons))
  expect_true(all(c("device_version", "pct_of_known") %in% names(ds$device_version_by_cohort)))
  expect_false(any(is.na(ds$device_switching_by_cohort$n_known)))
})

test_that("the combined figure builds and the stage 'device' writes its tables, the figure and the index entries", {
  cfg <- cfg_dev()
  s <- small_study(-0.5); persons <- s$inp$persons
  tabs <- build_device_tables(sim_device_rows(persons, cfg), persons, cfg)
  ds <- device_summary(tabs$person, persons, cfg)
  p <- plot_device_composition(ds, cfg); expect_s3_class(p, "ggplot"); invisible(ggplot2::ggplot_build(p))
  expect_false(any(vapply(p$layers, function(l) inherits(l$geom, "GeomErrorbar"), logical(1))))                # no interval bars on the barplot
  p2 <- plot_device_cohort_mix(ds, cfg); expect_s3_class(p2, "ggplot"); invisible(ggplot2::ggplot_build(p2))
  expect_equal(levels(ds$device_cohort_mix$category)[1L], "All people")
  expect_equal(ds$device_cohort_mix[category != "All people", sum(pct_of_category), by = category]$V1, rep(100, uniqueN(ds$device_cohort_mix$category) - 1L), tolerance = 1e-8)
  expect_equal(device_missing_label <- cfg$device$missing_label, "Missing")
  store <- as_day_store(s$inp$day_tables, persons, device = tabs)
  out <- file.path(tempdir(), "dev_stage"); unlink(out, recursive = TRUE)
  r <- suppressMessages(run_missingness_pipeline(store, persons, cfg, out, stages = "device"))
  expect_true(file.exists(file.path(out, "figures", "device_version_missing_and_composition.png")))
  expect_true(all(file.exists(file.path(out, "tables", paste0(c("device_missing_by_cohort", "device_tests", "device_composition_by_cohort", "device_missing_model"), ".tsv")))))
  fi <- fread(file.path(out, "figures", "figure_index.tsv")); expect_equal(fi[figure == "device_version_missing_and_composition.png", step], 3L)
  expect_true(file.exists(file.path(out, "figures", "device_cohort_mix_by_version.png")))
  expect_match(readLines(file.path(out, "figures", "figure_captions.md")) |> paste(collapse = " "), "cohort mix of all people")
  ti <- fread(file.path(out, "tables", "table_index.tsv")); expect_true(all(ti[grepl("^device_", table), step] == 3L))
  # without device information the stage says what to run
  expect_error(suppressMessages(run_missingness_pipeline(as_day_store(s$inp$day_tables, persons), persons, cfg, out, stages = "device")), "Cell 8b")
})

test_that("run_etl_device reads the four exports, writes device_person / device_timeline next to the day tables and fails with advice when an export is missing", {
  cfg <- cfg_dev(); cfg$export$bucket <- "gs://b"; cfg$export$owner <- "me"; cfg$export$date <- "20261006"
  s <- small_study(-0.5); persons <- s$inp$persons
  raw <- sim_device_rows(persons, cfg)
  by_co <- split(raw, persons$cohort[match(raw$person_id, persons$person_id)])
  paths <- vapply(cfg$cohorts, function(co) export_path(cfg, co, "device"), "")
  lister <- function(p) { co <- names(paths)[paths == p]; if (!length(co)) stop("no files found for: ", p); paste0(p, ".part") }
  reader <- function(f, ...) { co <- names(paths)[paste0(paths, ".part") == f]; by_co[[co]] }
  dir <- file.path(tempdir(), "dev_etl"); unlink(dir, recursive = TRUE); dir.create(dir); saveRDS(persons, file.path(dir, "persons.rds"))
  r <- suppressMessages(run_etl_device(cfg, dir, lister = lister, reader = reader))
  expect_true(all(file.exists(file.path(dir, c("device_person.rds", "device_timeline.rds", "device_status.tsv")))))
  expect_equal(nrow(readRDS(file.path(dir, "device_person.rds"))), nrow(persons))
  st <- open_day_store(cfg, dir) |> tryCatch(error = function(e) NULL)                          # no day tables here; the accessor is what matters
  expect_equal(nrow(r$person), nrow(persons))
  expect_error(suppressMessages(run_etl_device(cfg, dir, lister = function(p) stop("none"), reader = reader)), "Run Cell 6c")
})
