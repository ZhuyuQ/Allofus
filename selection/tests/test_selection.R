# Run from selection/:   Rscript tests/test_selection.R
suppressMessages({ library(testthat); library(data.table) })
root <- normalizePath("..")
miss <- if (dir.exists(file.path(root, "final_missing"))) "final_missing" else "missingness"   # Workbench / this repo
source(file.path(root, miss, "R/utils.R"))
source(file.path(root, "functions/cci_function.R"))
suppressMessages(source(file.path(root, "functions/ind_descriptive_func.R")))     # filter_pid_func, the earlier rule
source(file.path(root, "selection/selection_sql.R"))
source(file.path(root, "selection/selection_core.R"))

d0 <- as.Date("2020-01-01")
# observed day offsets -> long table
mk <- function(id, offs) data.table(person_id = id, date = d0 + offs)

test_that("wear_from_dates: days, span and longest interior gap", {
  w <- wear_from_dates(rbind(mk(1, c(0, 1, 5, 6, 7, 20)), mk(2, 3), mk(3, c(0, 0, 1))))
  expect_equal(w[person_id == 1, .(n_days, max_gap)], data.table(n_days = 6L, max_gap = 12L))   # days 8-19 are missing
  expect_equal(w[person_id == 2, max_gap], 0L)
  expect_equal(w[person_id == 3, n_days], 2L)                    # duplicate dates count once
})

test_that("the rules reproduce filter_pid_func, including the boundaries", {
  set.seed(7)
  cases <- list(
    `21 days`   = 0:20,                              # exactly 21 days, 100%
    `20 days`   = 0:19,
    `gap 14`    = c(0:10, 25:45),                    # 14 missing days in a row, 32 of 46 days = 69.6%
    `gap 14 ok` = c(0:20, 35:60),                    # 14 missing, 47 of 61 days = 77%
    `gap 15`    = c(0:20, 36:60),
    `70%`       = c(0:20, 22:29)[1:30][c(0:20, 22:29) < 30],     # filled below
    `69%`       = sort(c(0:19, 21:29))[1:29])
  cases[["70%"]]  <- setdiff(0:29, c(1, 4, 7, 10, 13, 16, 19, 22, 25))      # 21 of 30 = 70% exactly
  cases[["69%"]]  <- setdiff(0:29, c(1, 4, 7, 10, 13, 16, 19, 22, 25, 28))   # 20 of 30 -> below 21 days
  cases[["70% b"]] <- setdiff(0:99, sort(sample(1:98, 30)))                  # 70 of 100 = 70% (gaps may exceed 14)
  rnd <- lapply(1:150, function(i) {
    span <- sample(15:250, 1); p <- runif(1, 0.3, 1)
    o <- which(runif(span) < p) - 1L
    if (length(o) > 0 && runif(1) < 0.4) o <- o[!(o %in% sample(0:span, sample(5:30, 1)))]
    o
  })
  rnd <- Filter(length, rnd); names(rnd) <- paste0("r", seq_along(rnd))
  all_cases <- c(cases, rnd)
  long <- rbindlist(Map(function(o, i) mk(i, o), all_cases, seq_along(all_cases)))
  ref <- sort(suppressMessages(filter_pid_func(path.out = paste0(tempdir(), "/"), df = as.data.frame(long),
                                               cohort = "T", type = "x")))
  w <- sel_prepare_wear(cbind(wear_from_dates(long), cohort = "T", stream = "activity"))
  mine <- sel_select(w, data.table(cohort = "T", person_id = unique(long$person_id), cci_score = 0),
                     data.table(cohort = "T", person_id = unique(long$person_id)), streams = "activity")$pids$person_id
  expect_equal(sort(mine), ref)
  expect_true(all(c(1L, 6L) %in% mine))                         # 21 days and 70% exactly pass
})

test_that("a person must pass in every stream; compliance bands are [lo, hi)", {
  w <- rbindlist(list(
    data.table(cohort = "A", person_id = 1:3, stream = "activity", n_days = c(60L, 60L, 60L), max_gap = 2L, followup_days = c(100L, 100L, 100L)),
    data.table(cohort = "A", person_id = 1:3, stream = "sleep",    n_days = c(60L, 60L, 60L), max_gap = 2L, followup_days = c(100L, 100L, 100L)),
    data.table(cohort = "A", person_id = 1:2, stream = "hr",       n_days = c(60L, 80L),      max_gap = 2L, followup_days = c(100L, 100L))))
  pool <- data.table(cohort = "A", person_id = 1:4)
  cci  <- data.table(cohort = "A", person_id = 1:4, cci_score = c(0, 0, 0, 11))
  expect_equal(sel_select(w, cci, pool, comp_min = 60)$pids$person_id, 1:2)                   # 3 has no hr data, 4 has no data
  expect_equal(sel_select(w, cci, pool, comp_min = 50, comp_max = 60)$pids$person_id, integer())   # 60 is not < 60
  expect_equal(sel_select(w, cci, pool, comp_min = 60, comp_max = 61)$pids$person_id, 1L)           # person 2 has 80% in hr
  r <- sel_select(w, cci, pool, comp_min = 60)
  expect_equal(r$attrition[step == "pool", n], 4L)
  expect_true(all(diff(r$attrition$n) <= 0))                    # steps only remove people
  expect_equal(sel_select(w, cci, pool, comp_min = 60, cci_max = -1)$pids$person_id, integer())      # CCI 0 > -1
})

test_that("age is taken at the first record over the three streams", {
  p <- data.table(cohort = "A", person_id = 1:2, dob = as.Date(c("1960-06-15", "1960-06-15")))
  w <- data.table(cohort = "A", person_id = c(1, 1, 2), stream = c("activity", "hr", "sleep"),
                  first_date = as.Date(c("2015-06-14", "2012-01-01", "2020-06-15")))
  pp <- sel_persons(p, w)
  expect_equal(pp$first_date, as.Date(c("2012-01-01", "2020-06-15")))
  expect_equal(pp$age_first, c(51L, 60L))
})

test_that("CCI uses cci_calc_func with the new age", {
  codes <- sel_cci_codes()
  expect_gt(length(codes), 700); expect_true(all(c("57054005", "128462008") %in% codes))
  p <- data.table(cohort = c("A", "A", "A", "B"), person_id = c(1, 2, 3, 1), age_first = c(55L, 72L, 40L, 55L))
  cc <- data.table(cohort = c("A", "A", "A", "B"), person_id = c(1, 1, 3, 1),
                   standard_concept_code = c("57054005", "57054005", "128462008", "57054005"))   # MI twice; metastatic = 6
  r <- sel_cci(p, cc)
  expect_equal(r[cohort == "A" & person_id == 1, cci_score], 1 + 1)       # MI 1 + age 55 -> 1
  expect_equal(r[cohort == "A" & person_id == 2, cci_score], 3)           # no condition, age 72 -> 3
  expect_equal(r[cohort == "A" & person_id == 3, cci_score], 6)           # metastatic 6, age 40 -> 0
  expect_equal(r[cohort == "B" & person_id == 1, cci_score], 2)
})

test_that("the SQL carries the cohort definitions and the wear logic", {
  n_flags <- function(co) lengths(regmatches(sel_sql_pool(co), gregexpr("has_fitbit\\w* = 1", sel_sql_pool(co))))
  expect_equal(c(n_flags("Diabetes"), n_flags("Insomnia"), n_flags("Healthy"), n_flags("Prostate")), c(8L, 8L, 8L, 9L))
  pr <- sel_sql_pool("Prostate")
  expect_true(all(grepl(c("has_fitbit_device_information"), pr), grepl("4163261, 200970", pr), grepl("sex_at_birth_concept_id IN \\(45878463\\)", pr)))
  expect_false(grepl("cr.concept_id IN \\(45878463\\)", pr))            # sex is excluded from `person`, not through cb_criteria
  hy <- sel_sql_pool("Healthy")
  expect_true(grepl("person_id NOT IN", hy)); expect_false(grepl("include", hy))
  expect_equal(length(SEL_COHORTS$Healthy$exclude), 26L); expect_false(anyDuplicated(SEL_COHORTS$Healthy$exclude) > 0)
  expect_true(grepl("436962, 42538607, 37110488, 4102985, 439013", sel_sql_pool("Insomnia")))
  expect_true(grepl("LAG\\(dt\\)", sel_sql_wear("Diabetes", "sleep")) && grepl("sleep_daily_summary", sel_sql_wear("Diabetes", "sleep")))
  expect_true(grepl("heart_rate_summary", sel_sql_wear("Healthy", "hr")))
  expect_true(grepl("'57054005'", sel_sql_conditions("Diabetes", c("57054005", "1"))))
})

test_that("sel_compare and sel_verify", {
  a <- data.table(cohort = "A", person_id = 1:3); b <- data.table(cohort = "A", person_id = 2:4)
  expect_equal(sel_compare(a, b)[, .(n_both, only_new, only_old)], data.table(n_both = 2L, only_new = 1L, only_old = 1L))
  long <- rbind(mk(1, c(0, 1, 9)), mk(2, 0:5))
  w <- sel_prepare_wear(cbind(wear_from_dates(long), cohort = "A", stream = "activity"))
  expect_silent(invisible(capture.output(sel_verify(w, function(s, ids) long[person_id %in% ids]))))
  w[1, max_gap := 99L]
  expect_error(suppressMessages(capture.output(sel_verify(w, function(s, ids) long[person_id %in% ids]))), "differ")
})
cat("selection tests done\n")
