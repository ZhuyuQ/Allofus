# ==============================================================================
# runin.R
# ------------------------------------------------------------------------------
# RUN-IN ANALYSIS: does the wear of a person in the FIRST days of follow-up (the run-in) relate to the gaps that come LATER?
#
# Why not whole-period coverage / strata of it: coverage = observed days / days between the first and last record contains
# the gaps that are the outcome (a person with more gaps has lower coverage by construction), so stratifying or adjusting on it
# conditions on the outcome. The run-in measure is computed from days that all lie BEFORE the days that are modelled:
#
#   run-in            study days 1..N of the person in the modelled stream (N = cfg$runin$days, 28)
#   missing days m    N - observed days among them (0 .. N-1; the first record is observed by definition). NA if the person's
#                     follow-up is shorter than N days.
#   modelled days     onset model: prior days d with study_day > N (the outcome is a gap starting on d + 1);
#                     models B and C: gaps whose last observed day before the gap has study_day > N.
#   cov0              log2(1 + m), centred at its mean over the people analysed (one coefficient = per DOUBLING of 1 + m; the
#                     centring keeps the main effect of the steps term the effect at the average run-in)
#   cov0_cat          0 | 1 | 2-3 | >=4 missing days, reference 0 (shape check; onset only)
#   cov0_lin          100 * m / N / 10, centred: per +10 percentage points of the run-in missing (sensitivity; onset only)
#
# Hypotheses (stated before looking at the results):
#   H1  (wear trait -> later gaps)   more missing days in the run-in -> a higher odds that a gap starts later (cov0 main effect > 0).
#   H2  (effect modification)        the association between a higher-than-usual steps day and a gap tomorrow depends on the run-in
#                                    (key term x cov0). Models B and C: the same question for the gap-length term (B) and the
#                                    prior-day steps term (C). H2 is the confirmatory test (Holm over the streams).
#   H3  (robustness)                 the primary odds ratio of steps is the same without and with the run-in terms.
# Variants of the SAME model on the SAME rows (only rows after the run-in):
#   PRIMARY analysis (stage "runin"):
#   reference            no run-in term                     (onset, B, C)
#   interaction          + cov0 + key term x cov0           (onset, B, C)   <- the run-in and its interaction: H1 (cov0) and H2 (the product)
#   SENSITIVITY analysis (stage "runin_sensitivity"):
#   main                 + cov0                             (onset, B, C; the null model of the likelihood-ratio test of H2)
#   categorical          + cov0_cat                          (onset; shape check against the log scale by AIC / LRT)
#   linear_interaction   + cov0_lin + key term x cov0_lin   (onset; sensitivity to the functional form)
# Limits: the measure is noisy (28 days; attenuates), a gap in the run-in makes any_prior_gap = 1 for the rest of follow-up (the two
# terms overlap by construction: GVIF is reported), and the 70% observed-days selection of the cohort ties early to later coverage.
# ==============================================================================

# PRIMARY run-in analysis: the model without a run-in term and the model with the run-in AND its interaction with the key term.
# SENSITIVITY: every other form (run-in main effect only, categories, linear coding); they have their own stage, tables and figures.
RUNIN_PRIMARY_VARIANTS  <- c("reference", "interaction")
RUNIN_SENS_ONSET        <- c("main", "categorical", "linear_interaction")
RUNIN_SENS_GAP          <- "main"
RUNIN_ONSET_VARIANTS    <- c("reference", "main", "interaction", "categorical", "linear_interaction")      # display order
RUNIN_GAP_VARIANTS      <- c("reference", "main", "interaction")
RUNIN_LABELS <- c(reference = "Without run-in term", main = "+ run-in missing days", interaction = "+ run-in and its interaction",
                  categorical = "Run-in as categories", linear_interaction = "Run-in linear (per 10 points) + interaction")
RUNIN_SCHEMA <- 1L          # version of the stored run-in results (part of every checkpoint key)

# ---- the run-in measure ------------------------------------------------------------------------------------------------

# Per person: first / last record, and the missing days of the run-in (NA if the follow-up is shorter than n days).
#   days: day table of the stream (person_id, date; one row per observed day)
runin_person_table <- function(days, n) {
  d <- days[, .(person_id, date)]
  out <- d[, .(first_date = min(date), last_date = max(date), n_in = sum(date <= min(date) + n - 1L)), by = person_id]
  out[, runin_missing := fifelse(as.integer(last_date - first_date) + 1L >= n, as.integer(n - n_in), NA_integer_)]
  out[, n_in := NULL]
  out[]
}

# Labels of the missing-days categories from the upper bounds: c(0, 1, 3) -> "0", "1", "2-3", ">=4".
runin_category_labels <- function(cfg) {
  b <- cfg$runin$cat_breaks
  lo <- c(b[1L], b[-length(b)] + 1L)
  c(ifelse(lo == b, as.character(b), paste0(lo, "-", b)), paste0(">=", b[length(b)] + 1L))
}
runin_category <- function(m, cfg) {
  b <- cfg$runin$cat_breaks
  factor(runin_category_labels(cfg)[findInterval(m, b + 1L) + 1L], levels = runin_category_labels(cfg))
}

# Adds cov0 (centred log2(1 + m)), cov0_lin (centred per +10 points) and cov0_cat to a table with person_id and runin_missing, BY
# REFERENCE. The centres are person-level means over the people in the table; they are kept in attr(, "cov0_center").
add_runin_columns <- function(dt, n_days, cfg) {
  pm <- unique(dt[, .(person_id, m = runin_missing)])
  c_log <- mean(log2(1 + pm$m)); c_lin <- mean(100 * pm$m / n_days / 10)
  dt[, `:=`(cov0 = log2(1 + runin_missing) - c_log, cov0_lin = 100 * runin_missing / n_days / 10 - c_lin,
            cov0_cat = droplevels(runin_category(runin_missing, cfg)))]
  setattr(dt, "cov0_center", c(log = c_log, lin = c_lin))
  invisible(dt)
}

# Formula terms of a run-in variant for the key term of the model (onset: the primary within-person term). A term without variation
# in the data is left out (the variant then equals the reference; the coefficient table shows it).
runin_formula_terms <- function(variant, key_term, dat) {
  terms <- switch(variant, reference = character(), main = "cov0", interaction = c("cov0", paste0(key_term, ":cov0")),
                  categorical = "cov0_cat", linear_interaction = c("cov0_lin", paste0(key_term, ":cov0_lin")),
                  stop("unknown run-in variant: ", variant, call. = FALSE))
  col <- switch(variant, categorical = "cov0_cat", linear_interaction = "cov0_lin", "cov0")
  if (length(terms) && col %in% names(dat) && data.table::uniqueN(dat[[col]]) < 2L) character() else terms
}

# Log-likelihood, AIC and the pieces of the key-term x run-in contrast of one fitted model (one row), for the likelihood-ratio tests
# and the effect of the key term at a given number of run-in missing days.
runin_fit_summary <- function(fit, V, ftxt, key_term, center) {
  ll <- stats::logLik(fit)
  fe <- if (inherits(fit, "glmmTMB")) glmmTMB::fixef(fit)$cond else lme4::fixef(fit)
  V <- as.matrix(V)
  it <- grep(paste0("^", key_term, ":cov0(_lin)?$"), names(fe), value = TRUE)
  has_int <- length(it) == 1L && key_term %in% names(fe) && all(c(key_term, it) %in% rownames(V))
  data.table(loglik = as.numeric(ll), ll_df = attr(ll, "df"), aic = -2 * as.numeric(ll) + 2 * attr(ll, "df"),
             cov0_center = unname(center[["log"]]), cov0_center_lin = unname(center[["lin"]]),
             b_key = if (has_int) unname(fe[[key_term]]) else NA_real_, v_key = if (has_int) V[key_term, key_term] else NA_real_,
             b_int = if (has_int) unname(fe[[it]]) else NA_real_, v_int = if (has_int) V[it, it] else NA_real_,
             cov_key_int = if (has_int) V[key_term, it] else NA_real_, int_term = if (has_int) it else NA_character_)
}

# ---- descriptive statistics (no model) -----------------------------------------------------------------------------------

# What the run-in looks like and how it relates to the gaps that come later, for ONE stream.
#   bundle: prepare_stream_bundle() result; k: gap definition of the outcome (primary); term: the within-person steps term
# Returns list(flow, distribution, rate, steps_by_group, association), every table with a `stream` column.
runin_describe <- function(bundle, cfg, k = cfg$onset$primary_k, term = cfg$primary_term) {
  n_days <- cfg$runin$days; ycol <- paste0("y_ge", k); s <- bundle$stream
  rs <- bundle$rs
  pp_all <- unique(rs[, .(person_id, cohort, runin_missing)])
  after <- rs[rs$study_day > n_days & !is.na(rs$runin_missing) & !is.na(rs[[ycol]]), c("person_id", "cohort", "runin_missing", ycol, term), with = FALSE]
  setnames(after, c(ycol, term), c("y", "z"))
  after[, category := runin_category(runin_missing, cfg)]
  co_levels <- c(cfg$cohorts, "All cohorts")

  # sample flow: who is in the run-in analysis
  analysed <- after[, .(n_rows = .N, n_events = sum(y)), by = .(person_id, cohort, runin_missing, category)]
  flow <- rbindlist(lapply(c(cfg$cohorts, "All cohorts"), function(co) {
    sel <- function(x) if (co == "All cohorts") x else x[as.character(cohort) == co]
    a <- sel(pp_all); b <- sel(analysed)
    data.table(stream = s, cohort = co, n_people = nrow(a), n_complete_runin = sum(!is.na(a$runin_missing)), n_analysed = nrow(b),
               n_rows_after_runin = sum(b$n_rows), n_events_after_runin = sum(b$n_events))
  }))

  # distribution of the missing days in the run-in
  dist <- rbindlist(lapply(co_levels, function(co) {
    x <- if (co == "All cohorts") analysed else analysed[as.character(cohort) == co]
    if (nrow(x) == 0L) return(NULL)
    cnt <- x[, .N, by = category]
    out <- data.table(stream = s, cohort = co, category = factor(runin_category_labels(cfg), levels = runin_category_labels(cfg)))
    out[, n_people := cnt$N[match(as.character(category), as.character(cnt$category))]]
    out[is.na(n_people), n_people := 0L]
    out[, `:=`(pct = 100 * n_people / sum(n_people), total_people = sum(n_people),
               median_missing_days_in_category = vapply(as.character(category), function(cl) { v <- x$runin_missing[as.character(x$category) == cl]; if (length(v)) as.numeric(stats::median(v)) else NA_real_ }, numeric(1)),
               mean_missing_days = mean(x$runin_missing), pct_no_missing_day = 100 * mean(x$runin_missing == 0L))]
    out
  }))

  # H1, descriptive: the daily rate of gap onsets AFTER the run-in by run-in category (cluster-robust interval, people are the clusters)
  rate <- rbindlist(lapply(co_levels, function(co) {
    x <- if (co == "All cohorts") after else after[as.character(cohort) == co]
    rbindlist(lapply(runin_category_labels(cfg), function(cl) {
      xx <- x[as.character(category) == cl]
      if (nrow(xx) < 2L) return(NULL)
      ci <- cluster_mean_ci(xx$y, xx$person_id)
      data.table(stream = s, cohort = co, category = cl, n_people = ci$n_persons, n_rows = ci$n_gaps, n_events = sum(xx$y),
                 rate = 100 * ci$mean, rate_low = 100 * max(0, ci$lo), rate_high = 100 * ci$hi)
    }))
  }))

  # H2, descriptive: the observed gap rate by bin of the prior-day steps z-score, separately for each run-in category and cohort
  nb <- cfg$runin$deciles
  zz <- after[!is.na(z)]
  br <- if (nrow(zz)) unique(stats::quantile(zz$z, probs = seq(0, 1, length.out = nb + 1L), names = FALSE)) else numeric()
  steps <- if (length(br) >= 3L) {
    zz[, bin := findInterval(z, br, rightmost.closed = TRUE, all.inside = TRUE)]
    zz[, {
      ci <- cluster_mean_ci(y, person_id)
      .(mean_z = mean(z), n_rows = .N, n_events = sum(y), n_people = ci$n_persons, rate = 100 * ci$mean, rate_low = 100 * max(0, ci$lo), rate_high = 100 * ci$hi)
    }, by = .(cohort = as.character(cohort), category = as.character(category), bin)][, stream := s][order(cohort, category, bin)]
  } else data.table()

  # person-level rank correlation: missing days in the run-in vs the person's gap-onset rate after it (descriptive, Holm over the rows)
  assoc <- rbindlist(lapply(co_levels, function(co) {
    x <- (if (co == "All cohorts") analysed else analysed[as.character(cohort) == co])[n_rows >= 30L]
    if (nrow(x) < 20L || data.table::uniqueN(x$runin_missing) < 3L) return(NULL)
    r <- x$n_events / x$n_rows
    ct <- suppressWarnings(stats::cor.test(x$runin_missing, r, method = "spearman", exact = FALSE))
    ci <- spearman_ci(unname(ct$estimate), nrow(x))
    data.table(stream = s, cohort = co, n_people = nrow(x), rho = unname(ct$estimate), rho_low = ci[1L], rho_high = ci[2L], p = ct$p.value)
  }))
  if (nrow(assoc)) assoc[, p_holm := stats::p.adjust(p, "holm")]
  list(flow = flow, distribution = dist, rate = rate, steps_by_group = steps, association = assoc)
}

# ---- the units: description and models of ONE stream ----------------------------------------------------------------------

# Fits (or resumes) every run-in unit of one stream. Each fit is its own checkpoint, so one failure never costs the others.
#   bundle   lazy(prepare_stream_bundle(...))          gap_data  lazy(build_gap_level_data(...)) or NULL
#   with_models FALSE: only the description unit;  sensitivity TRUE: also the other forms (main, categorical, linear)
run_runin_stage <- function(s, bundle, gap_data, grid, cfg, ckpt_dir, dkey, with_models = TRUE, sensitivity = FALSE) {
  status <- character()
  n_days <- cfg$runin$days
  unit <- function(kind, id, compute, spec) {
    u <- run_unit(ckpt_dir, list(kind = kind, stream = s, id = id),
                  ckpt_key(cfg, kind, dkey, c(spec, list(schema = RUNIN_SCHEMA)), cfg$model$sample_persons), compute, cfg)
    status[paste(kind, id, sep = ": ")] <<- u$status
  }
  prim <- grid[stream == s & role == "primary"]
  if (nrow(prim) != 1L) return(status)
  unit("runin_extra", "describe", function() runin_describe(bundle$get(), cfg, k = prim$k, term = cfg$primary_term),
       list(k = prim$k, term = cfg$primary_term, what = "describe"))
  if (!with_models) return(status)

  v_on  <- c(RUNIN_PRIMARY_VARIANTS, if (sensitivity) RUNIN_SENS_ONSET)
  v_gap <- c(RUNIN_PRIMARY_VARIANTS, if (sensitivity) RUNIN_SENS_GAP)
  for (v in v_on) local({
    vv <- v
    unit("runin_onset", paste(prim$spec_id, vv, sep = " | "), function() {
      r <- fit_onset_spec(bundle$get(), prim, cfg, variant = list(runin_days = n_days, runin_terms = vv)); r$fit <- NULL; r
    }, c(spec_key(prim), list(runin = vv)))
  })
  if (!is.null(gap_data)) {
    for (v in v_gap) local({
      vv <- v
      unit("runin_gapB", paste("B_back2_adj", vv, sep = " | "), function()
        isolated(cfg, fit_gap_dose_response(gap_data$get(), cfg, "back2_dev", adjusted = TRUE, adjust_pre = FALSE, id = "B_back2_adj", runin = vv)),
        list(model = "B_back2_adj", runin = vv, gap_schema = GAP_SCHEMA))
      unit("runin_gapC", paste("C_adj", vv, sep = " | "), function()
        isolated(cfg, fit_gap_length_ordinal(gap_data$get(), cfg, adjusted = TRUE, id = "C_adj", runin = vv)),
        list(model = "C_adj", runin = vv, gap_schema = GAP_SCHEMA))
    })
  }
  status
}

# ---- reading the results -----------------------------------------------------------------------------------------------------

# Everything finished, all streams: coefficients, model information, diagnostics, and the derived tables
# (runin_tests, runin_shape, runin_effect_by_missing, runin_primary_effect, runin_quality).
collect_runin_results <- function(ckpt_dir, cfg, data_keys) {
  grab <- function(cu, part) rbindlist(Map(function(v, st) {
    x <- v[[part]]; if (is.null(x) || !is.data.frame(x) || nrow(x) == 0L) return(NULL)
    if (!"stream" %in% names(x)) x <- cbind(stream = st, x)
    x[]
  }, cu$values, cu$units$stream), fill = TRUE)
  on <- collect_units(ckpt_dir, "runin_onset", cfg, data_keys); ex <- collect_units(ckpt_dir, "runin_extra", cfg, data_keys)
  gb <- collect_units(ckpt_dir, "runin_gapB", cfg, data_keys); gc <- collect_units(ckpt_dir, "runin_gapC", cfg, data_keys)
  out <- list(
    runin_flow = grab(ex, "flow"), runin_distribution = grab(ex, "distribution"), runin_rate = grab(ex, "rate"),
    runin_steps_by_group = grab(ex, "steps_by_group"), runin_association = grab(ex, "association"),
    runin_onset_coefficients = grab(on, "coef"), runin_onset_models = grab(on, "info"),
    runin_onset_collinearity = grab(on, "gvif"), runin_onset_term_tests = grab(on, "wald"),
    runin_gapB = grab(gb, "coef"), runin_gapB_all_terms = grab(gb, "coef_all"),
    runin_gapC = grab(gc, "coef"), runin_gapC_all_terms = grab(gc, "coef_all"),
    runin_gap_collinearity = rbindlist(list(grab(gb, "gvif"), grab(gc, "gvif")), fill = TRUE),
    runin_gap_term_tests = rbindlist(list(grab(gb, "wald"), grab(gc, "wald")), fill = TRUE),
    runin_failures = rbindlist(list(ex$failures, on$failures, gb$failures, gc$failures), fill = TRUE))
  out$runin_primary_effect     <- runin_primary_effect(out, cfg)
  out$runin_tests              <- runin_tests(out, cfg)
  out$runin_shape              <- runin_shape(out, cfg)
  out$runin_effect_by_missing  <- runin_effect_by_missing(out, cfg)
  out$runin_quality            <- runin_quality(out, cfg)
  out
}

# The key term of each model, in every variant: onset = the primary within-person term (OR per +1 SD higher), model C = the prior-day
# steps term (OR of a longer gap per +1 SD), model B = the gap-length term (SD change in activity per doubling of the gap).
runin_primary_effect <- function(res, cfg) {
  one <- function(tab, term_name, model_name, is_or) {          # (the gap tables have a column `model`: the argument must not share its name)
    if (!nz(tab) || !"runin_variant" %in% names(tab)) return(NULL)
    x <- tab[term == term_name & !is.na(estimate)]
    if (nrow(x) == 0L) return(NULL)
    if (is_or) x[, .(model = model_name, stream, variant = runin_variant, term = term_name, effect = exp(estimate), low = exp(conf_low), high = exp(conf_high), p, status, measure = "odds ratio")]
    else       x[, .(model = model_name, stream, variant = runin_variant, term = term_name, effect = estimate, low = conf_low, high = conf_high, p, status, measure = "SD change per doubling")]
  }
  out <- rbindlist(list(one(res$runin_onset_coefficients, cfg$primary_term, "onset", TRUE), one(res$runin_gapB, "log2_gap", "model B", FALSE),
                        one(res$runin_gapC, "w_pre_day_steps", "model C", TRUE)), fill = TRUE)
  if (nrow(out)) out[, variant := factor(variant, levels = RUNIN_ONSET_VARIANTS)]
  out[]
}

# The run-in terms as hypothesis tests.
#   analysis "primary": the model with the run-in and its interaction. H1 = cov0 (the run-in effect at the average steps, since the steps term is
#     centred), H2 = key term x cov0; Wald z. Holm over the streams within a model.
#   analysis "sensitivity" (needs stage runin_sensitivity): the likelihood-ratio test of H2 (interaction vs the main-effect model, same rows),
#     H1 from the main-effect model (Wald and LRT against the reference), the linear coding, and the shape of the coding (categories vs reference, LRT).
#   Models B and C: Wald z (the B fit is REML).
runin_tests <- function(res, cfg) {
  rows <- list()
  pick <- function(tab, st, v, term_name) { x <- tab[stream == st & runin_variant == v & term == term_name]; if (nrow(x)) x[1L] }
  lrt <- function(mods, st, big, small) {
    a <- mods[stream == st & runin_variant == big & !is.na(loglik)]; b <- mods[stream == st & runin_variant == small & !is.na(loglik)]
    if (nrow(a) != 1L || nrow(b) != 1L) return(NULL)
    chisq <- max(0, 2 * (a$loglik - b$loglik)); df <- a$ll_df - b$ll_df
    list(chisq = chisq, df = df, p = if (df > 0) stats::pchisq(chisq, df, lower.tail = FALSE) else NA_real_)
  }
  add <- function(model, st, analysis, hyp, form, term_label_, x, is_or, test = NULL, per = "doubling of (1 + missing days)") {
    if (is.null(x) && is.null(test)) return(invisible(NULL))
    rows[[length(rows) + 1L]] <<- data.table(model = model, stream = st, analysis = analysis, hypothesis = hyp, form = form, term = term_label_,
      effect_measure = if (is.null(x)) NA_character_ else paste(if (is_or) "odds ratio per" else "SD change per", per),
      estimate = x$estimate %||% NA_real_, se = x$se %||% NA_real_, effect = if (is.null(x)) NA_real_ else if (is_or) exp(x$estimate) else x$estimate,
      low = if (is.null(x)) NA_real_ else if (is_or) exp(x$conf_low) else x$conf_low, high = if (is.null(x)) NA_real_ else if (is_or) exp(x$conf_high) else x$conf_high,
      p_wald = x$p %||% NA_real_, chisq_lrt = test$chisq %||% NA_real_, df_lrt = test$df %||% NA_real_, p_lrt = test$p %||% NA_real_, status = x$status %||% "ok")
  }
  key <- c(onset = cfg$primary_term, `model B` = "log2_gap", `model C` = "w_pre_day_steps")
  H1 <- "H1: run-in missing days -> later gap onset"; H2 <- "H2: key term depends on run-in"
  if (nz(res$runin_onset_coefficients)) {
    cf <- res$runin_onset_coefficients; md <- res$runin_onset_models
    for (st in unique(cf$stream)) {
      it <- paste0(key[["onset"]], ":cov0")
      add("onset", st, "primary", H1, "log2 coding, with interaction", "cov0", pick(cf, st, "interaction", "cov0"), TRUE)
      add("onset", st, "primary", H2, "log2 coding, with interaction", it, pick(cf, st, "interaction", it), TRUE)
      add("onset", st, "sensitivity", H1, "log2 coding, main effect only", "cov0", pick(cf, st, "main", "cov0"), TRUE, lrt(md, st, "main", "reference"))
      l2 <- lrt(md, st, "interaction", "main")                                                              # needs the main-effect model (sensitivity stage)
      if (!is.null(l2)) add("onset", st, "sensitivity", H2, "log2 coding, LRT vs the main-effect model", it, pick(cf, st, "interaction", it), TRUE, l2)
      itl <- paste0(key[["onset"]], ":cov0_lin")
      add("onset", st, "sensitivity", H1, "linear coding (per 10 points), with interaction", "cov0_lin", pick(cf, st, "linear_interaction", "cov0_lin"), TRUE, per = "+10 percentage points missing in the run-in")
      add("onset", st, "sensitivity", H2, "linear coding (per 10 points), with interaction", itl, pick(cf, st, "linear_interaction", itl), TRUE, per = "+10 percentage points missing in the run-in")
      sh <- lrt(md, st, "categorical", "reference")
      if (!is.null(sh)) add("onset", st, "sensitivity", "Shape: categories vs reference", "categories (0 | 1 | 2-3 | >=4), LRT", "cov0_cat", NULL, TRUE, sh)
    }
  }
  for (spec in list(list(tab = res$runin_gapB_all_terms, model = "model B", is_or = FALSE), list(tab = res$runin_gapC_all_terms, model = "model C", is_or = TRUE))) {
    cf <- spec$tab
    if (!nz(cf) || !"runin_variant" %in% names(cf)) next
    it <- paste0(key[[spec$model]], ":cov0")
    for (st in unique(cf$stream)) {
      add(spec$model, st, "primary", H1, "log2 coding, with interaction", "cov0", pick(cf, st, "interaction", "cov0"), spec$is_or)
      add(spec$model, st, "primary", H2, "log2 coding, with interaction", it, pick(cf, st, "interaction", it), spec$is_or)
      add(spec$model, st, "sensitivity", H1, "log2 coding, main effect only", "cov0", pick(cf, st, "main", "cov0"), spec$is_or)
    }
  }
  out <- rbindlist(rows, fill = TRUE)
  if (nrow(out) == 0L) return(out)
  out[, p := fifelse(!is.na(p_lrt), p_lrt, p_wald)]
  out[!is.na(p), p_holm := stats::p.adjust(p, "holm"), by = .(model, analysis, hypothesis, form)]
  out[model == "model B", hypothesis := sub("key term", "gap-length effect", hypothesis, fixed = TRUE)]
  out[model != "model B", hypothesis := sub("key term", "steps effect", hypothesis, fixed = TRUE)]
  out[]
}

# Functional form of the run-in (onset only): log-likelihood, AIC and the likelihood-ratio tests of each coding against the reference.
runin_shape <- function(res, cfg) {
  md <- res$runin_onset_models
  if (!nz(md) || !"loglik" %in% names(md)) return(data.table())
  x <- md[!is.na(loglik), .(stream, variant = runin_variant, n_rows, n_events, loglik, ll_df, aic)]
  x[, `:=`(delta_aic_vs_main = aic - aic[variant == "main"][1L], delta_aic_vs_reference = aic - aic[variant == "reference"][1L]), by = stream]
  x[, `:=`(chisq_vs_reference = pmax(0, 2 * (loglik - loglik[variant == "reference"][1L])), df_vs_reference = ll_df - ll_df[variant == "reference"][1L]), by = stream]
  x[variant == "reference", `:=`(chisq_vs_reference = NA_real_, df_vs_reference = NA_real_)]
  x[, p_vs_reference := fifelse(!is.na(df_vs_reference) & df_vs_reference > 0, stats::pchisq(chisq_vs_reference, df_vs_reference, lower.tail = FALSE), NA_real_)]
  x[, variant := factor(variant, levels = RUNIN_ONSET_VARIANTS)]
  x[order(stream, variant)][]
}

# Effect of the key term at given numbers of missing days in the run-in, from the interaction variant:
#   log effect(m) = b_key + b_int * (log2(1 + m) - centre), variance v_key + L^2 v_int + 2 L cov.
runin_effect_by_missing <- function(res, cfg, at = c(0L, 1L, 3L, 7L)) {
  rows <- list()
  from_info <- function(info, model, is_or) {
    if (!nz(info) || !all(c("b_key", "cov0_center") %in% names(info))) return(NULL)
    x <- info[runin_variant == "interaction" & !is.na(b_int)]
    if (nrow(x) == 0L) return(NULL)
    rbindlist(lapply(seq_len(nrow(x)), function(i) rbindlist(lapply(at, function(m) {
      L <- log2(1 + m) - x$cov0_center[i]
      est <- x$b_key[i] + x$b_int[i] * L; se <- sqrt(x$v_key[i] + L^2 * x$v_int[i] + 2 * L * x$cov_key_int[i])
      data.table(model = model, stream = x$stream[i], missing_days = m, estimate = est, se = se,
                 effect = if (is_or) exp(est) else est, low = if (is_or) exp(est - 1.96 * se) else est - 1.96 * se,
                 high = if (is_or) exp(est + 1.96 * se) else est + 1.96 * se, measure = if (is_or) "odds ratio per +1 SD" else "SD change per doubling of the gap")
    }))))
  }
  gap_info <- function(coef, model, is_or) {
    if (!nz(coef) || !"b_int" %in% names(coef)) return(NULL)
    x <- coef[runin_variant == "interaction" & !is.na(b_int)]
    if (nrow(x) == 0L) return(NULL)
    x <- x[, .(stream, b_key = estimate, v_key = se^2, b_int, v_int = se_int^2, cov_key_int, cov0_center)]
    from_info(cbind(x, runin_variant = "interaction"), model, is_or)
  }
  out <- rbindlist(list(from_info(res$runin_onset_models, "onset", TRUE), gap_info(res$runin_gapB, "model B", FALSE), gap_info(res$runin_gapC, "model C", TRUE)), fill = TRUE)
  out[]
}

# Sparse-data and convergence flags of every run-in fit: events (onset) or gaps (B, C) per estimated coefficient.
runin_quality <- function(res, cfg) {
  rows <- list()
  if (nz(res$runin_onset_models)) {
    m <- res$runin_onset_models
    npar <- res$runin_onset_coefficients[term != "(Intercept)", .(n_coefficients = .N), by = .(stream, runin_variant)]
    m <- merge(m, npar, by = c("stream", "runin_variant"), all.x = TRUE)
    rows[[1L]] <- m[, .(model = "onset", stream, variant = runin_variant, n_people = n_persons, n_events, n_coefficients, events_per_coefficient = n_events / n_coefficients,
                        status, singular = if ("singular" %in% names(m)) singular else NA)]
  }
  for (nm in c("runin_gapB", "runin_gapC")) {
    m <- res[[nm]]; ta <- res[[paste0(nm, "_all_terms")]]
    if (!nz(m) || !nz(ta)) next
    npar <- ta[, .(n_coefficients = .N), by = .(stream, model_id, runin_variant)]
    m <- merge(m, npar, by = c("stream", "model_id", "runin_variant"), all.x = TRUE)
    rows[[length(rows) + 1L]] <- m[, .(model = if (nm == "runin_gapB") "model B" else "model C", stream, variant = runin_variant, n_people = n_persons, n_events = n_gaps,
                                       n_coefficients, events_per_coefficient = n_gaps / n_coefficients, status, singular = if ("singular" %in% names(m)) singular else NA)]
  }
  out <- rbindlist(rows, fill = TRUE)
  if (nrow(out) == 0L) return(out)
  out[, sparse := events_per_coefficient < cfg$runin$epv_min]
  out[, usable := status %in% c("ok", "ok_bootstrap_se") & !sparse]
  out[]
}
