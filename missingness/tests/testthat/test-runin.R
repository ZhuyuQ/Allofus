# Run-in analysis (runin.R): the measure, the rows that are modelled, and the variants of the model.

# three people with a known run-in (28 days): person 1 misses days 3 and 10 of the run-in, person 2 has a 20-day follow-up only, person 3 misses nothing
runin_days_fixture <- function() {
  d0 <- as.Date("2023-03-01")
  one <- function(id, drop, len) data.table(person_id = id, date = d0 + setdiff(0:(len - 1L), drop - 1L))
  rbind(one(1L, c(3L, 10L), 60L), one(2L, integer(), 20L), one(3L, integer(), 45L))
}

test_that("missing days of the run-in: counted among the first N days only, NA when the follow-up is shorter, the first record is day 1", {
  rt <- runin_person_table(runin_days_fixture(), 28L)
  expect_equal(rt[order(person_id), runin_missing], c(2L, NA, 0L))
  # observed days AFTER the run-in never change the measure (no look-ahead)
  d <- runin_days_fixture(); d2 <- d[!(person_id == 1L & date > as.Date("2023-03-01") + 27L)]
  expect_equal(runin_person_table(d2, 28L)[person_id == 1L, runin_missing], 2L)
  # a gap that starts inside the run-in and runs on past it counts only its days inside the run-in
  d3 <- d[!(person_id == 3L & date >= as.Date("2023-03-01") + 26L & date <= as.Date("2023-03-01") + 40L)]
  expect_equal(runin_person_table(d3, 28L)[person_id == 3L, runin_missing], 2L)
  expect_equal(runin_person_table(d, 14L)[order(person_id), runin_missing], c(2L, 0L, 0L))            # N is a setting
})

test_that("categories and the centred log scale", {
  cfg <- make_config()
  expect_equal(runin_category_labels(cfg), c("0", "1", "2-3", ">=4"))
  expect_equal(as.character(runin_category(c(0L, 1L, 2L, 3L, 4L, 12L), cfg)), c("0", "1", "2-3", "2-3", ">=4", ">=4"))
  expect_equal(levels(runin_category(0L, cfg)), c("0", "1", "2-3", ">=4"))
  dt <- data.table(person_id = rep(1:4, c(5, 3, 8, 1)), runin_missing = rep(c(0L, 1L, 3L, 9L), c(5, 3, 8, 1)))
  add_runin_columns(dt, 28L, cfg)
  expect_equal(mean(unique(dt[, .(person_id, cov0)])$cov0), 0, tolerance = 1e-12)                      # centred over PEOPLE, not over rows
  expect_equal(dt$cov0 - dt$cov0[1L], log2(1 + dt$runin_missing) - log2(1), tolerance = 1e-12)           # per doubling of (1 + missing days)
  expect_equal(unname(attr(dt, "cov0_center")[["log"]]), mean(log2(1 + c(0, 1, 3, 9))), tolerance = 1e-12)
  expect_equal(mean(unique(dt[, .(person_id, cov0_lin)])$cov0_lin), 0, tolerance = 1e-12)
  expect_equal(runin_formula_terms("interaction", "w_act_steps", dt), c("cov0", "w_act_steps:cov0"))
  expect_equal(runin_formula_terms("reference", "w_act_steps", dt), character())
  expect_equal(runin_formula_terms("categorical", "w_act_steps", dt), "cov0_cat")
  expect_error(runin_formula_terms("nope", "w", dt), "unknown run-in variant")
  dt2 <- copy(dt)[, cov0 := 1]; expect_equal(runin_formula_terms("main", "w", dt2), character())            # a term without variation is left out
})

test_that("the onset data of the run-in analysis hold only days AFTER the run-in, the same rows in every variant", {
  s <- small_study(-0.5); cfg <- s$cfg; b <- s$b
  expect_true("runin_missing" %in% names(b$rs))
  spec <- primary_spec(cfg)
  dref <- build_onset_data(b, spec, cfg, list(runin_days = 28L, runin_terms = "reference"))
  dint <- build_onset_data(b, spec, cfg, list(runin_days = 28L, runin_terms = "interaction"))
  expect_equal(nrow(dref), nrow(dint))
  expect_true(all(c("cov0", "cov0_lin", "cov0_cat") %in% names(dint)))
  all_rows <- build_onset_data(b, spec, cfg)
  expect_lt(nrow(dref), nrow(all_rows))
  # every row of the run-in data is a day after the run-in of a person with a complete run-in
  n_after <- sum(b$rs$study_day > 28L & !is.na(b$rs$runin_missing))
  expect_lte(nrow(dref), n_after); expect_gt(nrow(dref), 0.95 * n_after)
  expect_equal(mean(unique(dint[, .(person_id, cov0)])$cov0), 0, tolerance = 1e-10)
  expect_false("runin_missing" %in% names(build_onset_data(b, spec, cfg)))                              # the other analyses never see it
})

test_that("onset variants: same rows, the interaction is named by its terms, the LRT and the effect at a given run-in follow from the stored pieces", {
  s <- small_study(-0.5); cfg <- s$cfg; b <- s$b; spec <- primary_spec(cfg)
  fits <- lapply(setNames(RUNIN_ONSET_VARIANTS, RUNIN_ONSET_VARIANTS), function(v)
    suppressWarnings(fit_onset_spec(b, spec, cfg, variant = list(runin_days = 28L, runin_terms = v))))
  expect_true(all(vapply(fits, function(f) f$info$status == "ok", TRUE)))
  expect_equal(uniqueN(vapply(fits, function(f) f$info$n_rows, 1)), 1L)
  expect_false("cov0" %in% fits$reference$coef$term); expect_true("cov0" %in% fits$main$coef$term)
  expect_true(all(c("cov0", "w_act_steps:cov0") %in% fits$interaction$coef$term))
  expect_true(any(grepl("^cov0_cat", fits$categorical$coef$term))); expect_true(all(c("cov0_lin", "w_act_steps:cov0_lin") %in% fits$linear_interaction$coef$term))
  expect_equal(fits$interaction$info$ll_df - fits$main$info$ll_df, 1)                                 # nested: one more parameter
  expect_gte(fits$interaction$info$loglik, fits$main$info$loglik - 1e-6)
  # the stored pieces reproduce the coefficients
  i <- fits$interaction
  expect_equal(i$info$b_key, i$coef[term == "w_act_steps", estimate]); expect_equal(i$info$b_int, i$coef[term == "w_act_steps:cov0", estimate])
  expect_equal(i$info$v_int, i$coef[term == "w_act_steps:cov0", se]^2, tolerance = 1e-6)
  # a variant that changes only the run-in terms leaves the other estimates close
  expect_equal(fits$main$coef[term == "w_act_steps", estimate], fits$reference$coef[term == "w_act_steps", estimate], tolerance = 0.1)
  # the interaction is nearly uncorrelated with its main terms because the run-in is centred
  gv <- fits$interaction$gvif
  expect_lt(gv[term == "w_act_steps:cov0", vif_equiv], 5); expect_lt(gv[term == "cov0", vif_equiv], 5)
})

test_that("a planted wear trait is found by the run-in term (H1) and absent when nothing is planted", {
  cfg <- cfg_test()
  inp <- suppressMessages(simulate_inputs(cfg, n_per_cohort = c(Insomnia = 100, Diabetes = 100, Healthy = 200, Prostate = 50), seed = 21))
  thin <- function(strength, seed) {
    set.seed(seed)
    u <- rnorm(uniqueN(inp$persons$person_id)); names(u) <- inp$persons$person_id
    lapply(inp$day_tables["activity"], function(d) {
      p_miss <- plogis(-2.2 + strength * u[as.character(d$person_id)])                           # the person's own share of missing days: the SAME trait drives the run-in and the later gaps
      first <- !duplicated(d$person_id)
      d[first | stats::runif(nrow(d)) > p_miss]
    })
  }
  fit_cov0 <- function(strength, seed) {
    dt <- thin(strength, seed)
    b <- suppressMessages(prepare_stream_bundle("activity", c(dt, inp$day_tables[c("sleep", "hr")]), inp$persons, cfg, trailing = FALSE))
    spec <- primary_spec(cfg)
    r <- suppressWarnings(fit_onset_spec(b, spec, cfg, variant = list(runin_days = 28L, runin_terms = "main")))
    r$coef[term == "cov0", .(estimate, p)]
  }
  planted <- fit_cov0(1.0, 3); null <- fit_cov0(0.0, 3)
  expect_gt(planted$estimate, 0.1); expect_lt(planted$p, 1e-4)
  expect_lt(abs(null$estimate), 0.15); expect_gt(null$p, 0.001)
})

test_that("gap-level data carry the run-in; models B and C use only gaps after the run-in, the same gaps in every variant", {
  s <- small_study(-0.5); cfg <- s$cfg
  days <- s$inp$day_tables$activity
  g <- build_gap_level_data(days, days, s$inp$persons, cfg)
  expect_true(all(c("runin_missing", "prev_study_day") %in% names(g)))
  ref <- fit_gap_dose_response(g, cfg, "back2_dev", adjusted = TRUE, id = "B_back2_adj", runin = "reference")
  mn  <- fit_gap_dose_response(g, cfg, "back2_dev", adjusted = TRUE, id = "B_back2_adj", runin = "main")
  int <- fit_gap_dose_response(g, cfg, "back2_dev", adjusted = TRUE, id = "B_back2_adj", runin = "interaction")
  all_g <- fit_gap_dose_response(g, cfg, "back2_dev", adjusted = TRUE, id = "B_back2_adj")
  expect_equal(c(ref$coef$n_gaps, mn$coef$n_gaps), rep(int$coef$n_gaps, 2L)); expect_lt(ref$coef$n_gaps, all_g$coef$n_gaps)
  expect_true("cov0" %in% mn$coef_all$term); expect_true(all(c("cov0", "log2_gap:cov0") %in% int$coef_all$term))
  expect_true(all(c("b_int", "se_int", "cov_key_int", "cov0_center") %in% names(int$coef)))
  expect_false("cov0" %in% ref$coef_all$term)
  expect_error(fit_gap_dose_response(g, cfg, "back2_dev", adjusted = TRUE, id = "x", runin = "nope"), "unknown run-in variant")
  cr <- suppressWarnings(fit_gap_length_ordinal(g, cfg, adjusted = TRUE, id = "C_adj", runin = "interaction"))
  if (cr$coef$status %in% c("ok", "ok_bootstrap_se")) {
    expect_true(all(c("cov0", "w_pre_day_steps:cov0") %in% cr$coef_all$term)); expect_null(cr$pred)           # no predicted probabilities for these fits
    expect_equal(unique(cr$coef_all$runin_variant), "interaction")
  } else expect_true(cr$coef$status %in% c("po_rejected", "not_estimable"))
  # models fitted with the run-in have their own checkpoint keys
  expect_false(identical(ckpt_cfg_part(cfg, "runin_gapB"), ckpt_cfg_part(cfg, "gapB")))
  cfg2 <- cfg; cfg2$runin$days <- 14L
  expect_false(identical(ckpt_cfg_part(cfg, "runin_onset"), ckpt_cfg_part(cfg2, "runin_onset")))
  expect_identical(ckpt_cfg_part(cfg, "onset"), ckpt_cfg_part(cfg2, "onset"))                            # the other analyses do not depend on the run-in length
})

test_that("tests, shape, effect by missing days and quality tables are built from the stored fits; primary and sensitivity rows are kept apart", {
  s <- small_study(-0.5); cfg <- s$cfg; b <- s$b; spec <- primary_spec(cfg)
  fits <- lapply(setNames(RUNIN_ONSET_VARIANTS, RUNIN_ONSET_VARIANTS), function(v)
    suppressWarnings(fit_onset_spec(b, spec, cfg, variant = list(runin_days = 28L, runin_terms = v))))
  res <- list(runin_onset_coefficients = rbindlist(lapply(fits, `[[`, "coef"), fill = TRUE), runin_onset_models = rbindlist(lapply(fits, `[[`, "info"), fill = TRUE))
  tt <- runin_tests(res, cfg)
  prim <- tt[analysis == "primary"]
  expect_setequal(prim$hypothesis, c("H1: run-in missing days -> later gap onset", "H2: steps effect depends on run-in"))
  expect_true(all(prim$form == "log2 coding, with interaction") && all(is.na(prim$p_lrt)))                              # the primary tests are Wald tests of the interaction model
  it <- fits$interaction$coef
  expect_equal(prim[term == "cov0", estimate], it[term == "cov0", estimate]); expect_equal(prim[term == "w_act_steps:cov0", p_wald], it[term == "w_act_steps:cov0", p])
  sens <- tt[analysis == "sensitivity"]
  lr <- sens[grepl("LRT vs the main-effect model", form)]
  chisq <- 2 * (fits$interaction$info$loglik - fits$main$info$loglik)
  expect_equal(lr$chisq_lrt, max(0, chisq)); expect_equal(lr$p_lrt, stats::pchisq(max(0, chisq), 1, lower.tail = FALSE)); expect_equal(lr$p, lr$p_lrt)
  expect_true(any(grepl("categories", sens$form)) && any(grepl("linear", sens$form)) && any(grepl("main effect only", sens$form)))
  expect_equal(sens[grepl("categories", form), df_lrt], nlevels(runin_category(0:10, cfg)) - 1)
  expect_equal(prim$p_holm, prim$p)                                                                                   # one stream: Holm changes nothing
  # only the primary variants: no sensitivity rows
  res0 <- list(runin_onset_coefficients = res$runin_onset_coefficients[runin_variant %in% RUNIN_PRIMARY_VARIANTS], runin_onset_models = res$runin_onset_models[runin_variant %in% RUNIN_PRIMARY_VARIANTS])
  expect_true(all(runin_tests(res0, cfg)$analysis %in% c("primary", "sensitivity")))
  expect_equal(nrow(runin_tests(res0, cfg)[analysis == "primary"]), 2L)
  sh <- runin_shape(res, cfg)
  expect_equal(sh[variant == "reference", delta_aic_vs_reference], 0); expect_equal(sh[variant == "main", delta_aic_vs_main], 0)
  expect_equal(sh[variant == "categorical", df_vs_reference], nlevels(runin_category(0:10, cfg)) - 1)
  eb <- runin_effect_by_missing(res, cfg)
  expect_equal(sort(eb$missing_days), c(0L, 1L, 3L, 7L))
  i <- fits$interaction$info; L <- log2(1 + 3) - i$cov0_center
  expect_equal(eb[missing_days == 3L, estimate], i$b_key + i$b_int * L)
  expect_equal(eb[missing_days == 3L, se], sqrt(i$v_key + L^2 * i$v_int + 2 * L * i$cov_key_int))
  pe <- runin_primary_effect(res, cfg)
  expect_setequal(pe$variant, RUNIN_ONSET_VARIANTS); expect_true(all(pe$model == "onset"))
  q <- runin_quality(res, cfg); expect_true(all(c("events_per_coefficient", "sparse", "usable") %in% names(q)))
  # diagnostics of the run-in models are their own tables; one figure per model is flagged
  cl <- build_runin_collinearity(list(runin_onset_collinearity = rbindlist(lapply(fits, `[[`, "gvif"), fill = TRUE)))
  expect_true(all(cl[figure == TRUE, variant] == "interaction")); expect_equal(unique(cl[figure == TRUE, figure_id]), "onset")
  expect_true(all(c("main", "categorical", "linear_interaction") %in% cl$variant) && !any(cl[variant != "interaction", figure]))
  expect_equal(RUNIN_PRIMARY_VARIANTS, c("reference", "interaction")); expect_false(any(RUNIN_PRIMARY_VARIANTS %in% c(RUNIN_SENS_ONSET, RUNIN_SENS_GAP)))
})
