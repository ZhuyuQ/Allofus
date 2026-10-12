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

# ---- % missing by race / sex / CCI / age group -----------------------------------------------------------------------

test_that("rank_test_groups uses Kruskal-Wallis for >= 3 groups and Wilcoxon for 2, with the documented effect sizes", {
  set.seed(1); v <- c(rnorm(40, 0), rnorm(40, 0.8), rnorm(40, 0.3)); g <- rep(c("a", "b", "c"), each = 40)
  r3 <- rank_test_groups(v, g); kw <- kruskal.test(v ~ factor(g))
  expect_equal(r3$test, "Kruskal-Wallis"); expect_equal(r3$statistic, unname(kw$statistic)); expect_equal(r3$df, 2)
  expect_equal(r3$p, kw$p.value); expect_equal(r3$effect, unname(kw$statistic) / (120 - 1)); expect_equal(r3$effect_name, "epsilon-squared")
  r2 <- rank_test_groups(v[g != "c"], g[g != "c"]); w <- wilcox.test(v[g == "a"], v[g == "b"], exact = FALSE)
  expect_equal(r2$test, "Wilcoxon rank-sum"); expect_equal(r2$statistic, unname(w$statistic)); expect_equal(r2$p, w$p.value)
  expect_equal(r2$effect, cliffs_delta(v[g == "a"], v[g == "b"])); expect_lt(r2$effect, 0)      # group a is lower than b
  few <- rank_test_groups(c(v[g == "a"], v[g == "b"][1:3]), c(g[g == "a"], g[g == "b"][1:3]), min_n = 5L)
  expect_equal(few$test, "none"); expect_equal(few$n_groups, 1L)                              # groups below min_n are left out
  expect_true(is.na(rank_test_groups(c(1, NA, 3), c("a", "a", NA))$p))
})

test_that("stratified % missing: both test families, Holm within stream x stratum x family, and tidy summaries", {
  cfg <- cfg_test(); inp <- suppressMessages(simulate_inputs(cfg, n_per_cohort = c(Insomnia = 60, Diabetes = 60, Healthy = 120, Prostate = 30), seed = 5))
  gaps <- lapply(inp$day_tables, get_gap_runs)
  pg <- rbindlist(lapply(names(gaps), function(s) with_cohort(summarise_person_gaps(inp$day_tables[[s]], gaps[[s]]), inp$persons, s)))
  ms <- missing_pct_by_stratum(pg, inp$persons)
  expect_equal(nrow(ms$long), nrow(pg)); expect_setequal(ms$summary$stratum_var, c("sex", "age_group", "race_plot", "cci_group"))
  tt <- ms$tests
  expect_setequal(tt$family, c("groups_within_cohort", "cohorts_within_level"))
  one <- tt[stream == "sleep" & stratum_var == "age_group" & family == "groups_within_cohort" & test != "none"]
  expect_equal(one$p_holm, p.adjust(one$p, "holm"))
  # a test reproduces direct computation: Healthy cohort, sleep stream, age groups
  d <- ms$long[stream == "sleep" & cohort == "Healthy" & !is.na(age_group)]
  expect_equal(tt[stream == "sleep" & stratum_var == "age_group" & family == "groups_within_cohort" & group == "Healthy", p],
               kruskal.test(value ~ droplevels(age_group), d)$p.value)
  expect_equal(tt[stream == "sleep" & stratum_var == "age_group" & family == "cohorts_within_level" & group == "45-64", n_groups] > 0, TRUE)
})

test_that("odds ratios of the missing share by stratum level agree with a direct quasi-binomial fit and carry the omnibus F test", {
  cfg <- cfg_test(); inp <- suppressMessages(simulate_inputs(cfg, n_per_cohort = c(Insomnia = 80, Diabetes = 80, Healthy = 160, Prostate = 40), seed = 6))
  gaps <- lapply(inp$day_tables, get_gap_runs)
  pg <- rbindlist(lapply(names(gaps), function(s) with_cohort(summarise_person_gaps(inp$day_tables[[s]], gaps[[s]]), inp$persons, s)))
  m <- missing_pct_stratum_models(pg, inp$persons)
  row <- m[stratum_var == "age_group" & stream == "sleep" & cohort == "Healthy"]
  d <- merge(pg[stream == "sleep" & cohort == "Healthy", .(person_id, n_missing, n_observed)], inp$persons[, .(person_id, age_group)], by = "person_id")
  fit <- glm(cbind(n_missing, n_observed) ~ droplevels(age_group), data = d, family = quasibinomial())
  expect_equal(row$odds_ratio, unname(exp(coef(fit)[-1L]))); expect_equal(row$reference[1L], levels(droplevels(d$age_group))[1L])
  expect_equal(row$omnibus_p[1L], anova(fit, test = "F")[["Pr(>F)"]][2L]); expect_equal(row$dispersion[1L], summary(fit)$dispersion)
  expect_true(all(row$or_low < row$odds_ratio & row$odds_ratio < row$or_high))
  one <- m[stratum_var == "age_group" & stream == "sleep", .(p = omnibus_p[1L], ph = omnibus_p_holm[1L]), by = cohort]
  expect_equal(one$ph, p.adjust(one$p, "holm"))
  expect_false(any(m[stratum_var == "sex" & cohort == "Prostate"]$stream == "never"))                  # a one-sex cohort simply has no sex rows
  expect_equal(nrow(m[stratum_var == "sex" & cohort == "Prostate"]), 0L)
})
