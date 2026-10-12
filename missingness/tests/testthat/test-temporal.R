test_that("weekday counts equal a brute-force count", {
  set.seed(2)
  first <- as.Date("2023-03-01") + sample(0:60, 40, TRUE); W <- sample(1:100, 40, TRUE)
  wc <- window_weekday_counts(first, W)
  for (i in seq_along(first)) {
    truth <- table(factor(DOW_LABELS[wday(seq(first[i], by = "day", length.out = W[i]))], levels = DOW_LABELS))
    got <- wc[row == i]
    expect_equal(as.integer(got$den), as.integer(truth))
  }
})

make_pg <- function(days, persons_cohort = "Healthy") {
  g <- get_gap_runs(days); pg <- summarise_person_gaps(days, g); pg[, `:=`(cohort = factor(persons_cohort), stream = "activity")]
  list(pg = pg, g = g)
}

test_that("day-of-week, study-week and calendar-year counts add up to the window and the gaps", {
  days <- random_days(20, seed = 3); x <- make_pg(days)
  dw <- dow_counts(x$pg, x$g); wk <- study_week_counts(x$pg, x$g); yr <- calendar_year_counts(x$pg, x$g)
  tot <- x$pg[, .(win = sum(window_days), mis = sum(n_missing))]
  expect_equal(sum(dw$den), tot$win); expect_equal(sum(dw$num), tot$mis)
  expect_equal(sum(wk$den), tot$win); expect_equal(sum(wk$num), tot$mis)
  expect_equal(sum(yr$den), tot$win); expect_equal(sum(yr$num), tot$mis)
})

test_that("transition counts: continuation days + onsets = missing days", {
  days <- random_days(20, seed = 4); x <- make_pg(days)
  tc <- transition_person_counts(x$pg)
  expect_true(all(tc$continue_days + tc$onsets == tc$miss_days))
  tr <- transition_table(tc, cfg_test())
  expect_true(all(tr$p_onset >= 0 & tr$p_onset <= 1))
  expect_equal(tr$p_onset, sum(tc$onsets) / sum(tc$at_risk))
})
