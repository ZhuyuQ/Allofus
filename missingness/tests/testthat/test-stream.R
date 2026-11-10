# Stream states: which streams are present on each day (R/stream_states.R, R/plots_stream_states.R). The new code is sourced here (the loader is unchanged).
for (f in c("R/stream_states.R", "R/plots_stream_states.R")) source(file.path(if (file.exists(f)) "." else "../..", f))

d0 <- as.Date("2024-03-04")            # a Monday

# A hand-built study of five people; every expected number below is counted by hand.
make_study <- function() {
  ids <- 1:5; days <- function(first, last) seq(d0 + first, d0 + last, by = "day")
  grid <- rbindlist(lapply(ids, function(i) data.table(person_id = i, date = days(0, 9))))
  act <- copy(grid)[!(person_id == 5L & date < d0 + 2L)][, `:=`(steps = 5000, sedentary_minutes = 600, lightly_active_minutes = 200, fairly_active_minutes = 10, very_active_minutes = 5, mvpa_minutes = 15)]
  hr <- copy(grid)[, total_hr_min := 1300]; slp <- copy(grid)[, `:=`(has_main_sleep = TRUE, has_nap = FALSE)]
  # person 2: day 6 sedentary the whole day and no heart rate (filtered); day 7 sedentary the whole day with heart rate (flagged); no sleep on day 6
  act[person_id == 2L & date == d0 + 6L, `:=`(steps = 0, sedentary_minutes = 1440, lightly_active_minutes = 0, fairly_active_minutes = 0, very_active_minutes = 0, mvpa_minutes = 0)]
  act[person_id == 2L & date == d0 + 7L, `:=`(steps = 0, sedentary_minutes = 1440, lightly_active_minutes = 0, fairly_active_minutes = 0, very_active_minutes = 0, mvpa_minutes = 0)]
  hr <- hr[!(person_id == 2L & date == d0 + 6L)]; slp <- slp[!(person_id == 2L & date == d0 + 6L)]
  # person 1: day 6 no heart rate and no sleep (state 3), day 7 no sleep (state 4); in the middle, so that every stream stays active to the end of the window
  hr <- hr[!(person_id == 1L & date == d0 + 6L)]; slp <- slp[!(person_id == 1L & date %in% (d0 + 6:7))]
  persons <- data.table(person_id = ids, cohort = "Healthy", age_group = "18-44", device_family = "Charge")
  dev <- data.table(person_id = c(1L, 1L, 2L, 3L, 5L), device_date = as.character(d0 + c(0, 4, 5, 5, 0)), device_version = c("Charge", "Versa 3", "Versa 3", "Zip", "Versa 3"), device_type = "TRACKER")
  list(store = as_day_store(list(activity = act, sleep = slp, hr = hr), persons), persons = persons, dev = dev)
}

test_that("the length of a day: 1,380 minutes where clocks go forward, 1,500 where they go back", {
  m <- st_day_minutes(as.Date(c("2024-03-10", "2024-11-03", "2024-06-01", "2021-03-14", "2021-11-07")))
  expect_equal(m$min, c(1380L, 1440L, 1440L, 1380L, 1440L)); expect_equal(m$max, c(1440L, 1500L, 1440L, 1440L, 1500L))
})

test_that("device versions without a heart-rate sensor", {
  cfg <- st_config()
  v <- c("Charge", "Charge 5", "Alta", "Alta HR", "Inspire", "Inspire 2", "Fitbit Zip", "MobileTrack", "Versa 3", "Surge", NA)
  expect_equal(st_version_has_hr(v, cfg), c(FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, NA))
})

test_that("a hand-built study: windows, the three removal reasons in order, the activity filter and the states", {
  s <- make_study(); cfg <- cfg_test(); st <- st_config()
  trk <- st_trackers(s$dev, s$persons, cfg)
  D <- st_build_days(s$store, s$persons, trk, st, cfg)
  expect_equal(nrow(D), 50L)                                                                    # ten days each: the window is first to last record in any stream
  expect_equal(D[, .N, by = reason][order(reason)]$N, c(24L, 15L, 9L, 2L)); expect_equal(sort(unique(D$reason)), c("analysed", "device_unknown_removed", "device_without_hr", "stream_not_active"))
  expect_equal(D[reason == "stream_not_active", unique(person_id)], 5L)                         # activity starts on day 2
  expect_equal(D[reason == "device_without_hr", .N, by = person_id][order(person_id)]$N, c(4L, 5L))   # person 1 on the Charge, person 3 on the Zip
  expect_equal(D[reason == "device_unknown_removed", .N, by = person_id][order(person_id)]$N, c(5L, 10L))   # person 3 (usual device Zip), person 4 (no device at all)
  expect_equal(sum(D[person_id == 2L, unknown_kept]), 5L)                                       # person 2: unknown before the first record, usual device measures heart rate
  expect_equal(D[reason == "analysed", c(sum(act_filtered), sum(motionless_worn))], c(1L, 1L))
  expect_equal(D[reason == "analysed", .N, by = state][order(state)]$N, c(21L, 1L, 1L, 1L))     # states 1, 2, 3, 4 (no 'other')
  f <- D[person_id == 2L & date == d0 + 6L]; expect_equal(c(f$combo_before, f$combo), c("100", "000"))   # activity set to missing because no heart rate that day
  g <- D[person_id == 2L & date == d0 + 7L]; expect_equal(c(g$combo, g$state, g$motionless_worn), c("111", "1", "TRUE"))
  expect_equal(D[person_id == 1L & date == d0 + 6L, state], 3L); expect_equal(D[person_id == 1L & date == d0 + 7L, state], 4L)
  # the check without unknown-device days removes every such day
  st2 <- modifyList(st, list(drop_unknown_days = TRUE)); D2 <- st_assign_reason(copy(D), st2)
  expect_equal(D2[reason == "analysed", .N], 24L - 5L)
  # parts A and B on it
  a <- st_part_a(D, modifyList(st, list(min_people = 1L)))
  expect_equal(a$flow[step == "analysed days", n_days_at_step], 24L); expect_equal(a$flow[step == "removed: device without heart rate", n_days_at_step], 9L)
  b <- st_part_b(D[reason == "analysed"], modifyList(st, list(min_people = 1L)))
  expect_equal(b$people$n_people_affected, 1L); expect_equal(b$people$n_flagged_motionless, 1L)
})

test_that("a day on the fall-back date with 1,440 sedentary minutes covers the whole day for the shortest zone, a day with movement does not", {
  cfg <- st_config(); D <- data.table(act_row = TRUE, sedentary_minutes = c(1380, 1440, 1380, 1440), lightly_active_minutes = 0, fairly_active_minutes = 0, very_active_minutes = c(0, 0, 0, 30),
                                      hr_row = FALSE, hr_min = NA_real_, slp_row = FALSE, has_main_sleep = NA,
                                      date = as.Date(c("2024-03-10", "2024-03-10", "2024-06-01", "2024-06-01")))
  D <- st_flags(D, cfg)
  expect_equal(D$covers_day, c(TRUE, TRUE, FALSE, FALSE))      # 1,380 covers the day only where clocks go forward; 30 very active minutes are movement
})

test_that("the nesting checks use the data before the filter; the date check finds a shifted sleep table", {
  set.seed(4); n <- 40; days <- 60; d <- CJ(person_id = 1:n, i = 1:days)[, date := d0 + i - 1L]
  keep <- d[runif(.N) < 0.75]                                                                    # activity and heart rate: the same 75% of the days
  act <- copy(keep)[, `:=`(steps = 4000, sedentary_minutes = 600, mvpa_minutes = 10)]; hr <- copy(keep)[, total_hr_min := 1300]
  persons <- data.table(person_id = 1:n, cohort = "Healthy", age_group = "18-44", device_family = "Charge")
  dev <- data.table(person_id = 1:n, device_date = as.character(d0), device_version = "Versa 3", device_type = "TRACKER")
  run <- function(shift) {
    slp <- copy(keep)[, `:=`(date = date + shift, has_main_sleep = TRUE, has_nap = FALSE)]
    cfg <- cfg_test(); st <- st_config(); D <- st_build_days(as_day_store(list(activity = act, sleep = slp, hr = hr), persons), persons, st_trackers(dev, persons, cfg), st, cfg)
    list(D = D, abc = st_run_abc(D, modifyList(st, list(min_people = 1L))), st = st)
  }
  ok <- run(0L); bad <- run(1L)
  expect_true(all(ok$abc$gate[1L, passed])); expect_equal(unique(ok$abc$date_check$best_lag), 0L)
  expect_false(bad$abc$gate[1L, passed]); expect_true(all(bad$abc$date_check[shifted_stream == "sleep", best_lag] == 1L))      # taking the sleep of the NEXT day restores the nesting (the table is one day late)
  expect_error(st_require_gate(bad$abc$gate), "checks failed")
  ck <- ok$abc$nesting_checks[grepl("heart rate present", check)]; expect_equal(ck$share_pct, 0)                                  # heart rate and activity share their days
})

test_that("heart-rate minutes above the length of the day are counted, 1,500 minutes are allowed on the fall-back date", {
  D <- data.table(person_id = 1:4, cohort = factor("Healthy"), hr_row = TRUE, hr_min = c(1441, 1500, 1500, 1000), day_max = c(1440L, 1500L, 1440L, 1440L))
  m <- st_minutes_check(D, modifyList(st_config(), list(min_people = 1L)))
  expect_equal(m$overall$n_events, 2L); expect_equal(m$overall$share_pct, 50)
})

test_that("state 3 against state 1 within the same type of day, people with at least three days of each", {
  mk <- function(id, steps3, wd3 = TRUE) {
    dd <- data.table(person_id = id, date = d0 + 0:27)[, dtype := fifelse(wday(date) %in% 2:6, "weekday", "weekend")]
    s3 <- if (wd3) which(dd$dtype == "weekday")[1:4] else which(dd$dtype == "weekend")[1:4]
    dd[, `:=`(state = 1L, steps = 8000, mvpa_minutes = 30)]; dd[s3, `:=`(state = 3L, steps = steps3, mvpa_minutes = 5)]; dd[]
  }
  Dk <- rbindlist(list(mk(1L, 2000), mk(2L, 3000), mk(3L, 1000, wd3 = FALSE)))[, `:=`(cohort = factor("Healthy"), act_after = TRUE)]
  cmp <- st_state3_compare(Dk, modifyList(st_config(), list(min_people = 1L)))
  wk <- cmp[metric == "steps" & day_type == "weekday"]; we <- cmp[metric == "steps" & day_type == "weekend"]
  expect_equal(wk$n_people, 2L); expect_equal(wk$mean, mean(c(2000, 3000)) - 8000); expect_equal(we$n_people, 1L); expect_equal(we$mean, 1000 - 8000)    # each state-3 day is compared with the same type of day
  expect_equal(cmp[metric == "mvpa" & day_type == "weekday", mean], 5 - 30)
})

test_that("pooled days and the person average weigh differently; habitual night non-wear needs enough worn days", {
  Dk <- rbind(data.table(person_id = 1L, state = 1L, n = 100)[rep(1, 100)], data.table(person_id = 2L, state = 4L, n = 10)[rep(1, 10)])[, `:=`(cohort = factor("Healthy"), age_group = "18-44", device_family = "Charge")]
  cfg <- modifyList(st_config(), list(min_people = 1L))
  comp <- st_comp_one(Dk, "all", cfg)[state == 1L]
  expect_equal(comp[grepl("^all days", weighting), share_pct], 100 * 100 / 110); expect_equal(comp[grepl("^person", weighting), share_pct], 50)
  nh <- st_night_habit(Dk[, state := fifelse(person_id == 2L, 4L, 1L)], cfg)
  expect_equal(nrow(nh$persons), 1L)                                                                        # person 2 has only 10 worn days (< 30)
})

test_that("nights need a worn previous day; 'probably worn' needs high heart-rate wear on both days", {
  D <- data.table(person_id = 1L, cohort = factor("Healthy"), reason = "analysed", state = c(1L, 4L, 4L, 2L, 4L), hr_min = 60 * c(23, 23, 21, 23, 23))
  N <- st_nights(D)
  expect_equal(nrow(N), 2L); expect_equal(N$state, c(4L, 4L)); expect_equal(N$hr_prev, c(23, 23)); expect_equal(N$hr_now, c(23, 21))           # day 5 follows a state-2 day: not a night
  f <- st_part_f(D, modifyList(st_config(), list(min_people = 1L)))$table
  expect_equal(f[cutoff_hours == 20, n_events], 2L); expect_equal(f[cutoff_hours == 22, n_events], 1L)
})

test_that("cells with 1-19 people are suppressed", {
  tab <- data.table(n_people = c(5L, 25L), share_pct = c(1, 2), low = c(0, 1), high = c(2, 3))
  s <- st_suppress(tab, c("share_pct", "low", "high"), st_config())
  expect_equal(s$people_label, c("<20", "25")); expect_true(all(is.na(s[1L, .(share_pct, low, high, n_people)]))); expect_equal(s$share_pct[2L], 2)
})

test_that("a simulated study: the removed days add up and every figure builds", {
  cfg <- cfg_test(); st <- modifyList(st_config(), list(min_people = 5L))
  inp <- suppressMessages(simulate_inputs(cfg, n_per_cohort = c(Insomnia = 40, Diabetes = 40, Healthy = 80, Prostate = 30), seed = 5))
  p <- copy(inp$persons)[, device_family := "Charge"]
  D <- st_build_days(as_day_store(inp$day_tables, inp$persons), p, st_trackers(sim_device_rows(inp$persons, cfg), p, cfg), st, cfg)
  abc <- st_run_abc(D, st); dg <- st_run_dg(D, st)
  fl <- abc$flow[cohort == "Healthy"]
  expect_equal(fl[step == "analysed days", n_days_at_step] + fl[step == "removed: a stream not active", n_days_at_step] + fl[step == "removed: device without heart rate", n_days_at_step] +
               fl[step == "removed: device unknown (usual device without heart rate, or none)", n_days_at_step], fl[step == "all days in the windows", n_days_at_step])
  expect_equal(abc$combos[stage == "after the filter" & cohort == "Healthy", sum(n_events)], abc$flow[cohort == "Healthy" & step == "analysed days", n_days_at_step])
  figs <- st_figures(abc, dg, cfg, st)
  expect_equal(length(figs), 18L); expect_true(all(vapply(figs, inherits, logical(1), what = "ggplot")))
  d <- tempfile(); st_write_tables(c(abc, dg), d); st_save_figures(figs[1:2], d)
  expect_false(any(grepl("person_share|steps_data|day_sample|nights|hr_hours|night_habit_persons", list.files(file.path(d, "tables")))))     # person-level data never written
  expect_true(file.exists(file.path(d, "figures", "stream_a_flow.png")))
})
