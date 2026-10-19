wear_fixture <- function(seed = 7) {
  s <- small_study(-0.5, seed = seed)
  d <- s$inp$day_tables$activity
  pg <- summarise_person_gaps(d, get_gap_runs(d))
  cls <- classify_wear_consistency(wear_coverage_table(pg), s$cfg)
  c(s, list(pg = pg, cls = cls))
}

test_that("wear consistency = observed days / days between first and last record, classified by a rule fixed in the config", {
  w <- wear_fixture()
  cov <- wear_coverage_table(w$pg)
  expect_equal(cov$coverage, w$pg$n_observed / w$pg$window_days)
  expect_equal(cov$coverage, 1 - w$pg$pct_missing / 100)                          # the complement of the % missing of the missingness summaries
  expect_true(all(cov$coverage > 0 & cov$coverage <= 1))
  cr <- attr(w$cls, "criteria")
  expect_equal(cr$method, "median"); expect_equal(cr$threshold, median(cov$coverage))
  expect_equal(as.character(w$cls$stratum), ifelse(w$cls$coverage >= cr$threshold, "Consistent", "Inconsistent"))   # >= the median = consistent
  expect_equal(c(cr$n_consistent, cr$n_inconsistent), as.integer(table(w$cls$stratum)))
  expect_equal(levels(w$cls$stratum), c("Consistent", "Inconsistent"))
  expect_match(cr$rule, "Consistent: observed-day coverage >=")
  # a fixed threshold overrides the median; ties at the cut-off belong to the consistent group; classification is deterministic
  cfg2 <- w$cfg; cfg2$wear$method <- "fixed"; cfg2$wear$threshold <- 0.95
  c2 <- classify_wear_consistency(cov, cfg2)
  expect_equal(attr(c2, "criteria")$threshold, 0.95); expect_equal(as.character(c2$stratum), ifelse(cov$coverage >= 0.95, "Consistent", "Inconsistent"))
  expect_identical(classify_wear_consistency(cov, w$cfg)$stratum, classify_wear_consistency(cov, w$cfg)$stratum)
  tie <- data.table(person_id = 1:6, window_days = 10L, n_observed = c(10L, 10L, 10L, 10L, 9L, 8L), n_gaps = 0L, gap_max = 0L, longest_streak = 10L)
  ct <- classify_wear_consistency(wear_coverage_table(cbind(tie, pct_missing = 0)), w$cfg)
  expect_equal(as.character(ct$stratum), c(rep("Consistent", 4), "Inconsistent", "Inconsistent"))   # median = 1: ties at the median are consistent
})

test_that("descriptive tables report the strata, their sizes and how their composition differs (standardized mean differences)", {
  w <- wear_fixture(); d <- wear_consistency_describe(w$cls, w$inp$persons, "activity", w$cfg)
  expect_setequal(d$summary$group, c("All", "Consistent", "Inconsistent"))
  expect_equal(d$summary[group == "All", n_people], nrow(w$cls))
  expect_equal(d$summary[group %in% WEAR_LEVELS, sum(n_people)], nrow(w$cls))
  expect_true(all(c("coverage_median", "gaps_per_100_observed_days", "age_mean", "pct_female", "pct_Healthy") %in% names(d$summary)))
  expect_gt(d$summary[group == "Inconsistent", gaps_per_100_observed_days], d$summary[group == "Consistent", gaps_per_100_observed_days])   # by construction
  # SMD: planted imbalance
  p <- copy(w$inp$persons); a <- w$cls$person_id[w$cls$stratum == "Consistent"]
  p[, age_first := ifelse(person_id %in% a, 70, 40) + rnorm(.N)]
  d2 <- wear_consistency_describe(w$cls, p, "activity", w$cfg)
  expect_gt(d2$composition[variable == "Age at first record (years)", smd], 5); expect_true(d2$composition[variable == "Age at first record (years)", imbalanced])
  expect_equal(smd_binary(c(1, 1, 1, 0), factor(c("Consistent", "Consistent", "Inconsistent", "Inconsistent"), levels = WEAR_LEVELS)), (1 - 0.5) / sqrt((0 + 0.25) / 2))
  expect_equal(smd_continuous(c(1, 3, 2, 4), factor(c("Consistent", "Consistent", "Inconsistent", "Inconsistent"), levels = WEAR_LEVELS)), (2 - 3) / sqrt((2 + 2) / 2))
  expect_equal(nrow(d$by_person), nrow(w$cls)); expect_true(all(c("stream", "person_id", "stratum", "coverage") %in% names(d$by_person)))
})

test_that("strata are compared formally: difference of two independent estimates, ratio of odds ratios, never by which stratum is significant", {
  cmp <- compare_strata(log(0.5), 0.1, log(0.8), 0.15)
  expect_equal(cmp$difference, log(0.8) - log(0.5)); expect_equal(cmp$difference_se, sqrt(0.1^2 + 0.15^2))
  expect_equal(cmp$z, cmp$difference / cmp$difference_se); expect_equal(cmp$p_wald, 2 * pnorm(-abs(cmp$z)))
  expect_equal(cmp$ratio_of_or, 0.8 / 0.5); expect_equal(cmp$ratio_low, exp(cmp$difference_low))
  swap <- compare_strata(log(0.8), 0.15, log(0.5), 0.1)
  expect_equal(swap$p_wald, cmp$p_wald); expect_equal(swap$difference, -cmp$difference)       # reference choice changes the sign, not the test
  # one significant and one non-significant stratum do NOT make a significant difference
  a <- compare_strata(log(0.5), 0.1, log(0.8), 0.4)                                            # stratum 1: z = -6.9; stratum 2: z = -0.56
  expect_gt(a$p_wald, 0.05)
  expect_null(compare_strata(0.1, 0.1, 0.2, 0.1, log_scale = FALSE)$ratio_of_or)              # linear model: a difference, no ratio
})

test_that("stratified onset data keep only the people of the stratum and carry the stratum for the interaction test", {
  w <- wear_fixture(); map <- w$cls[, .(person_id, stratum)]; spec <- as.list(primary_spec(w$cfg)); spec$predictors <- list(spec$predictors[[1]])
  full <- build_onset_data(w$b, spec, w$cfg)
  both <- build_onset_data(w$b, spec, w$cfg, list(strata = map))
  expect_equal(nrow(both), nrow(full)); expect_equal(levels(both$stratum), c("Consistent", "Inconsistent"))
  one <- build_onset_data(w$b, spec, w$cfg, list(strata = map, stratum = "Consistent"))
  expect_true(all(one$person_id %in% map$person_id[map$stratum == "Consistent"])); expect_equal(uniqueN(one$stratum), 1L)
  expect_equal(nrow(one) + nrow(build_onset_data(w$b, spec, w$cfg, list(strata = map, stratum = "Inconsistent"))), nrow(full))   # a partition: nothing lost, nothing twice
  r <- fit_onset_spec(w$b, primary_spec(w$cfg), w$cfg, variant = list(strata = map, stratum = "Inconsistent"))
  expect_equal(r$info$wear_stratum, "Inconsistent"); expect_equal(r$info$n_persons, uniqueN(w$b$rs$person_id[w$b$rs$person_id %in% map$person_id[map$stratum == "Inconsistent"] & complete.cases(w$b$rs[, .(w_act_steps)])]))
})

test_that("the interaction test finds a planted difference between strata and not a common effect", {
  mk <- function(beta, offset, seed) {
    cfg <- cfg_test(); inp <- simulate_inputs(cfg, n_per_cohort = c(Insomnia = 40, Diabetes = 40, Healthy = 90, Prostate = 30), beta = beta, seed = seed)
    inp$day_tables <- lapply(inp$day_tables, function(d) { d <- copy(d); d[, person_id := person_id + offset]; setkey(d, person_id, date); d })
    inp$persons <- copy(inp$persons)[, person_id := person_id + offset]; inp
  }
  run <- function(beta_a, beta_b) {
    A <- mk(beta_a, 0L, 11); B <- mk(beta_b, 100000L, 12); cfg <- cfg_test()
    tabs <- list(activity = rbind(A$day_tables$activity, B$day_tables$activity)); setkey(tabs$activity, person_id, date)
    persons <- rbind(A$persons, B$persons); setkey(persons, person_id)
    b <- suppressMessages(prepare_stream_bundle("activity", tabs, persons, cfg))
    map <- data.table(person_id = persons$person_id, stratum = factor(ifelse(persons$person_id < 100000L, "Consistent", "Inconsistent"), levels = WEAR_LEVELS))
    it <- test_group_interaction(b, primary_spec(cfg), cfg, by = "stratum", variant = list(strata = map))
    ra <- fit_onset_spec(b, primary_spec(cfg), cfg, variant = list(strata = map, stratum = "Consistent"))$coef[term == "w_act_steps"]
    rb <- fit_onset_spec(b, primary_spec(cfg), cfg, variant = list(strata = map, stratum = "Inconsistent"))$coef[term == "w_act_steps"]
    list(it = it, ra = ra, rb = rb)
  }
  diff <- run(-0.9, 0); same <- run(-0.5, -0.5)
  expect_equal(diff$it$status, "ok"); expect_equal(diff$it$df, 1L); expect_lt(diff$it$p, 0.001)
  expect_lt(diff$ra$or, 0.55); expect_gt(diff$rb$or, 0.8)                                       # OR per +1 SD higher steps: strong vs none
  expect_gt(same$it$p, 0.01)                                                                     # the same effect in both strata
  cmp <- compare_strata(same$ra$estimate, same$ra$se, same$rb$estimate, same$rb$se)
  expect_gt(cmp$p_wald, 0.01)
})

test_that("the wear stage fits the primary models within each stratum as checkpointed units, resumes, and reports sparse strata", {
  w <- wear_fixture(); ck <- tempfile(); st <- as_day_store(w$inp$day_tables, w$inp$persons); dk <- data_key(st, "activity")
  built <- 0; bundle <- lazy(function() { built <<- built + 1; w$b })
  gap <- lazy(function() build_gap_level_data(w$inp$day_tables$activity, w$inp$day_tables$activity, w$inp$persons, w$cfg))
  grid <- make_onset_grid(w$cfg)
  s1 <- suppressWarnings(suppressMessages(run_wear_stage("activity", bundle, gap, grid, w$cls, w$cfg, ck, dk)))
  expect_true(all(s1 == "computed")); expect_equal(built, 1)
  expect_true(all(c("wear_onset: activity|joint|k1|adj | Consistent", "wear_onset: activity|joint|k1|adj | Inconsistent",
                    "wear_extra: interaction | activity|joint|k1|adj", "wear_gapB: B_back2_adj | Consistent", "wear_gapC: C_adj | Inconsistent") %in% names(s1)))
  s2 <- suppressMessages(run_wear_stage("activity", lazy(function() stop("no")), lazy(function() stop("no")), grid, w$cls, w$cfg, ck, dk))
  expect_true(all(s2 == "resumed"))
  r <- collect_wear_results(ck, w$cfg, list(activity = dk))
  expect_setequal(r$wear_strata_onset_coefficients$wear_stratum, WEAR_LEVELS)
  expect_equal(r$wear_strata_interaction$by, "stratum"); expect_equal(r$wear_strata_interaction$df, 1L)
  cmp <- r$wear_strata_comparison
  expect_true(all(c("model", "stream", "effect_consistent", "effect_inconsistent", "ratio_of_or", "p_wald", "p_difference") %in% names(cmp)))
  o <- cmp[model == "onset (primary)"]
  expect_equal(o$p_difference, r$wear_strata_interaction$p)                                       # the formal test of the difference is the interaction LRT
  expect_equal(o$ratio_of_or, o$effect_inconsistent / o$effect_consistent, tolerance = 1e-9)
  expect_true(nrow(r$wear_strata_onset_collinearity) > 0L)
  q <- wear_strata_quality(r, w$cfg)
  expect_true(all(c("n_events", "n_coefficients", "events_per_coefficient", "sparse", "usable") %in% names(q)))
  cfg_sp <- w$cfg; cfg_sp$wear$epv_min <- 1e6
  expect_true(all(wear_strata_quality(r, cfg_sp)$sparse)); expect_false(any(wear_strata_quality(r, cfg_sp)$usable))
  # the classification rule is part of every stratified result: another rule never reuses these fits
  cfg_fx <- w$cfg; cfg_fx$wear$method <- "fixed"
  expect_equal(nrow(collect_wear_results(ck, cfg_fx, list(activity = dk))$wear_strata_onset_coefficients), 0L)
  # strata that are too small are described but not modelled
  cfg_small <- w$cfg; cfg_small$wear$min_persons_stratum <- 1e6
  s3 <- suppressMessages(run_wear_stage("activity", bundle, gap, grid, w$cls, cfg_small, tempfile(), dk))
  expect_true(all(s3 == "too few people"))
})
