test_that("Cliff's delta has the right sign and bounds", {
  expect_equal(cliffs_delta(c(4, 5, 6), c(1, 2, 3)), 1); expect_equal(cliffs_delta(c(1, 2, 3), c(4, 5, 6)), -1)
  expect_equal(cliffs_delta(1:10, 1:10), 0)
})

make_pg_cohorts <- function() {
  set.seed(3)
  mk <- function(co, p, offset, n = 60) {
    days <- rbindlist(lapply(seq_len(n), function(i) {
      W <- 100; obs <- c(1L, which(runif(W - 1L) < p) + 1L)
      data.table(person_id = offset + i, date = as.Date("2024-01-01") + obs - 1L)
    }))
    g <- get_gap_runs(days); pg <- summarise_person_gaps(days, g)
    pg[, `:=`(cohort = co, stream = "activity")]
    pg
  }
  pg <- rbind(mk("Healthy", 0.95, 1000L), mk("Diabetes", 0.80, 2000L))
  pg[, cohort := factor(cohort, levels = c("Insomnia", "Diabetes", "Healthy", "Prostate"))]
  pg
}

test_that("cohort comparison detects a real difference; table and fraction model are consistent", {
  pg <- make_pg_cohorts()
  cmp <- compare_cohorts(pg, "pct_missing", "activity")
  expect_lt(cmp$overall$p, 0.001); expect_equal(nrow(cmp$pairwise), 1L)
  expect_gt(cmp$pairwise$cliffs_delta, 0)                                  # Diabetes (group1) misses more than Healthy
  tab <- cohort_missingness_table(pg); expect_equal(nrow(tab), 2L)
  persons <- unique(pg[, .(person_id, cohort)])[, `:=`(sex = factor("Female", levels = c("Female", "Male", "Other")), age_first = 50L)]
  fm <- fit_missing_fraction(pg, persons, FALSE, ref = "Healthy")
  expect_gt(fm[term == "cohortDiabetes", or], 1)
})
