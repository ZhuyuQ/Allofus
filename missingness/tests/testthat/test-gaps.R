test_that("worked 10-day example gives the documented gaps and risk-set rows", {
  cfg <- cfg_test()
  d0 <- as.Date("2024-01-01")
  days <- data.table(person_id = 1, date = d0 + c(1,2,3,6,7,9,10) - 1)
  g <- get_gap_runs(days)
  expect_equal(g$gap_length, c(2L, 1L))
  expect_equal(g$gap_start, d0 + c(3, 7))
  expect_equal(g$prev_obs_date, d0 + c(2, 6)); expect_equal(g$next_obs_date, d0 + c(5, 8))
  rs <- build_risk_set(days, cfg)
  expect_equal(as.integer(rs$date - d0) + 1L, c(1, 2, 3, 6, 7, 9))          # days 4,5,8 and the last day 10 are not rows
  expect_equal(rs$y_ge1, c(0, 0, 1, 0, 1, 0)); expect_equal(rs$y_ge2, c(0, 0, 1, 0, 0, 0))
  expect_equal(rs$prev_gap_len, c(0, 0, 0, 2, 2, 1))
  expect_equal(rs$days_since_return, c(0, 0, 0, 0, 1, 0))
  expect_equal(rs$any_prior_gap, c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE))
  expect_equal(restrict_first_gap(rs)$date, rs$date[1:3])                     # up to and incl. first onset day
})

test_that("gap runs, expanded gap days and the risk set match a brute-force calendar grid", {
  cfg <- cfg_test(); days <- random_days(25, seed = 11)
  g <- get_gap_runs(days); gd <- expand_gap_days(g); rs <- build_risk_set(days, cfg)
  for (pid in unique(days$person_id)) {
    dts <- days[person_id == pid, date]
    ng <- naive_gaps(dts)
    gp <- g[person_id == pid]
    expect_equal(gp$gap_start, ng$gap_start); expect_equal(gp$gap_length, ng$gap_length)
    # missing days
    full <- seq(min(dts), max(dts), by = "day")
    expect_equal(sort(gd[person_id == pid, date]), full[!(full %in% dts)])
    # risk set: every observed day but the last; y = next day missing
    r <- rs[person_id == pid]
    obs_not_last <- dts[dts < max(dts)]
    expect_equal(r$date, obs_not_last)
    expect_equal(r$y_ge1, as.integer(!((obs_not_last + 1L) %in% dts)))
    # history: last gap ending before d
    for (i in seq_len(nrow(r))) {
      ended <- ng[(gap_start + gap_length - 1L) < r$date[i]]
      exp_len <- if (nrow(ended)) ended$gap_length[nrow(ended)] else 0L
      expect_equal(r$prev_gap_len[i], exp_len)
    }
  }
})

test_that("person gap summary is internally consistent and finds the longest streak", {
  days <- random_days(20, seed = 5); g <- get_gap_runs(days); pg <- summarise_person_gaps(days, g)
  expect_true(all(pg$n_missing == pg$n_missing_gaps))
  expect_true(all(pg$window_days == pg$n_observed + pg$n_missing))
  x <- days[person_id == pg$person_id[1]]
  expect_equal(pg$longest_streak[1], max(diff(c(0L, which(diff(c(as.integer(x$date), Inf)) != 1L)))))
})

test_that("accounting check stops on duplicated dates", {
  days <- data.table(person_id = 1, date = as.Date("2024-01-01") + c(0, 1, 1, 5))
  expect_error(summarise_person_gaps(days, get_gap_runs(days)), "inconsistent gap accounting")
})
