# ==============================================================================
# sensitivity.R
# ------------------------------------------------------------------------------
# TRUE sensitivity analyses of the primary onset model: the SAME question
# (does a lower-than-usual prior day raise the odds that a gap starts?) answered
# under different assumptions. Each returns one row for the primary term, so the
# results can be read side by side and plotted as a forest plot.
#
#   main                  the reference fit
#   trailing_baseline     "usual" = mean of the previous 14 observed days (no look-ahead)
#   exclude_partial_wear  drop prior days that look like partial wear (device on/off)
#   first_gap_only        only days up to each person's first gap (no gap history at all)
#   main_sleep_only       sleep only: a night counts as observed only if main sleep exists
#   conditional_logit     within-person comparison, no random-effect assumption
#   cluster_robust_logit  pooled logistic model, person-clustered standard errors
#   spline_vs_linear      likelihood-ratio test: is the within effect non-linear?
#   random_slope          likelihood-ratio test: does the effect differ between people?
#
# (The stricter gap definitions k = 2, 3, 7 are model D, in the onset grid.)
# ==============================================================================

# Rows of one sensitivity analysis: the tested term (`tested = TRUE`) AND every other within-person
# term of the same model (`w_` / `tw_`), so one analysis can be drawn for all the stream's predictors.
# The tested term alone is what `sensitivity.tsv` reports; the full set is `sensitivity_all_terms.tsv`.
extract_term_row <- function(res, analysis, which_term) {
  info <- res$info
  if (is.null(res$coef)) {
    return(data.table(analysis = analysis, term = which_term, estimate = NA_real_, se = NA_real_,
                      conf_low = NA_real_, conf_high = NA_real_, p = NA_real_,
                      or = NA_real_, or_low = NA_real_, or_high = NA_real_,
                      n_events = info$n_events, n_persons_event = info$n_persons_event,
                      status = info$status, note = info$note, tested = TRUE))
  }
  rows <- res$coef[grepl("^t?w_", term) | term == which_term]
  rows[, .(analysis = analysis, term, estimate, se, conf_low, conf_high, p, or, or_low, or_high,
           n_events, n_persons_event, status = info$status, note = info$note, tested = term == which_term)]
}

# Same, for an estimator that returns its own coefficient table.
term_row_from_table <- function(tab, analysis, which_term, ec, status = "ok", note = "") {
  if (nrow(tab[term == which_term]) == 0L) {
    return(data.table(analysis = analysis, term = which_term, estimate = NA_real_, se = NA_real_,
                      conf_low = NA_real_, conf_high = NA_real_, p = NA_real_, or = NA_real_,
                      or_low = NA_real_, or_high = NA_real_, n_events = ec$n_events,
                      n_persons_event = ec$n_persons_event, status = "term_missing", note = note, tested = TRUE))
  }
  rows <- add_odds_ratios(copy(tab[grepl("^t?w_", term) | term == which_term]))
  rows[, .(analysis = analysis, term, estimate, se, conf_low, conf_high, p, or, or_low, or_high,
           n_events = ec$n_events, n_persons_event = ec$n_persons_event, status = status, note = note,
           tested = term == which_term)]
}

# ---- alternative estimators on the same model data ----------------------------------------

# Conditional logistic regression: each person is compared only with themselves.
# Time-invariant terms (cohort, age, sex, between-person levels) cancel out.
# Uses the Efron approximation to the exact conditional likelihood (the exact method
# overflows with many days per person); this is accurate because events are rare
# relative to the days at risk within each person.
fit_conditional_logit <- function(dat, spec, ycol) {
  preds <- spec$predictors[[1L]]
  keep_ids <- dat[get(ycol) == 1L, unique(person_id)]
  d <- dat[person_id %in% keep_ids]
  tv <- if (spec$adjusted) c("dow_t", "splines::ns(study_day, 3)", "prev_gap_len_l",
                             "any_prior_gap", "days_since_return_l") else character()
  tv <- tv[tv != "dow_t" | uniqueN(d$dow_t) > 1L]
  f <- stats::as.formula(paste0(ycol, " ~ ", paste(c(paste0("w_", preds), tv, "survival::strata(person_id)"),
                                                   collapse = " + ")))
  fit <- with_attached("survival", survival::clogit(f, data = d, method = "efron"))
  cf <- summary(fit)$coefficients
  data.table(term = rownames(cf), estimate = cf[, "coef"], se = cf[, "se(coef)"], p = cf[, ncol(cf)])[
    , `:=`(conf_low = estimate - 1.96 * se, conf_high = estimate + 1.96 * se)]
}

# Pooled logistic regression (no random effect) with person-clustered robust SEs.
fit_cluster_robust_logit <- function(dat, formula_text) {
  f_fixed <- sub(" \\+ \\(1 \\| person_id\\)$", "", formula_text)
  fit <- stats::glm(stats::as.formula(f_fixed), data = dat, family = stats::binomial())
  V   <- sandwich::vcovCL(fit, cluster = dat$person_id)
  est <- stats::coef(fit); se <- sqrt(diag(V))
  p   <- 2 * stats::pnorm(abs(est / se), lower.tail = FALSE)
  data.table(term = names(est), estimate = est, se = se, p = p,
             conf_low = est - 1.96 * se, conf_high = est + 1.96 * se)
}

# Likelihood-ratio test between two fitted models.
lrt <- function(m0, m1, boundary = FALSE) {
  chisq <- max(0, 2 * (as.numeric(stats::logLik(m1)) - as.numeric(stats::logLik(m0))))
  df <- attr(stats::logLik(m1), "df") - attr(stats::logLik(m0), "df")
  p  <- stats::pchisq(chisq, df, lower.tail = FALSE)
  if (boundary) p <- 0.5 * p      # variance on the boundary: 50:50 mixture of chi-squares
  list(chisq = chisq, df = df, p = p)
}

# Likelihood-ratio test from two stored log-likelihoods (see fit_loglik).
lrt_ll <- function(l0, l1, boundary = FALSE) {
  chisq <- max(0, 2 * (l1$ll - l0$ll)); df <- l1$df - l0$df
  p <- stats::pchisq(chisq, df, lower.tail = FALSE)
  if (boundary) p <- 0.5 * p
  list(chisq = chisq, df = df, p = p)
}

# ---- the sensitivity suite --------------------------------------------------------------------

# Everything the alternative estimators share: the model data of the primary spec, its event counts
# and its formula text. Built inside the unit that needs it (a worker), never kept between units.
sensitivity_context <- function(bundle, spec_l, cfg) {
  dat  <- build_onset_data(bundle, spec_l, cfg)
  ycol <- attr(dat, "ycol")
  ec   <- event_counts(dat, ycol)
  drop_cols <- constant_adjuster_cols(dat, spec_l$adjusted)
  list(dat = dat, ycol = ycol, ec = ec,
       ftxt = build_onset_formula(spec_l$predictors[[1L]], spec_l$adjusted, spec_l$k, drop_cols = drop_cols))
}

# One-row table for a unit that cannot be estimated (too few events).
not_estimable_row <- function(analysis, term, ec) {
  data.table(analysis = analysis, term = term, estimate = NA_real_, se = NA_real_, conf_low = NA_real_,
             conf_high = NA_real_, p = NA_real_, or = NA_real_, or_low = NA_real_, or_high = NA_real_,
             n_events = ec$n_events, n_persons_event = ec$n_persons_event,
             status = "not_estimable", note = "too few events", tested = TRUE)
}

# A likelihood-ratio test reported in the same row layout (the test is in `p`).
lrt_row <- function(analysis, term, ec, t, meaning) {
  data.table(analysis = analysis, term = term, estimate = NA_real_, se = NA_real_, conf_low = NA_real_,
             conf_high = NA_real_, p = t$p, or = NA_real_, or_low = NA_real_, or_high = NA_real_,
             n_events = ec$n_events, n_persons_event = ec$n_persons_event, status = "ok",
             note = sprintf("LRT chisq=%.2f, df=%d (%s)", t$chisq, t$df, meaning), tested = TRUE)
}

# Version of the stored sensitivity rows (2: every within-person term is kept, not only the tested one).
# It is part of the checkpoint key, so rows stored by an older version are recomputed, never mixed in.
SENSITIVITY_SCHEMA <- 2L

# The sensitivity analyses of the PRIMARY spec of one stream, each an independent checkpointed unit
# (so a failed or slow analysis never costs the others).
#   bundle : lazy(prepare_stream_bundle(..., trailing = TRUE))   alt : lazy bundle of the main-sleep-only
#   rule (sleep only; NULL otherwise)      prim : the primary row of the onset grid
# Returns the run_unit() status of each analysis.
run_sensitivity_stage <- function(s, bundle, alt, prim, cfg, ckpt_dir, dkey, term = cfg$primary_term) {
  spec_l <- as.list(prim); spec_l$predictors <- list(spec_l$predictors[[1L]])
  tw_term <- sub("^w_", "tw_", term)
  status <- character()
  unit <- function(id, compute, spec_extra = NULL) {
    u <- run_unit(ckpt_dir, list(kind = "sensitivity", stream = s, id = id),
                  ckpt_key(cfg, "sensitivity", dkey, c(spec_key(prim), spec_extra, list(schema = SENSITIVITY_SCHEMA)), cfg$model$sample_persons),
                  compute, cfg)
    status[id] <<- u$status
  }
  refit <- function(id, variant, which_term = term, the_bundle = bundle)
    unit(id, function() extract_term_row(fit_onset_spec(the_bundle$get(), prim, cfg, variant), id, which_term))

  # "main" is the primary onset model itself: it shares that unit's checkpoint (no second fit)
  cfg_p <- cfg; cfg_p$model$sample_persons <- effective_sample(cfg, prim$role)
  u <- run_unit(ckpt_dir, list(kind = "onset", stream = s, id = prim$spec_id),
                ckpt_key(cfg_p, "onset", dkey, spec_key(prim), cfg_p$model$sample_persons),
                function() { r <- fit_onset_spec(bundle$get(), prim, cfg_p); r$fit <- NULL; r }, cfg)
  status["main"] <- u$status
  if (!is.null(u$value)) unit("main", function() extract_term_row(u$value, "main", term))

  refit("trailing_baseline", list(baseline = "trailing"), tw_term)
  refit("exclude_partial_wear", list(drop_partial = TRUE))
  refit("first_gap_only", list(first_gap = TRUE))
  if (s == "sleep" && !is.null(alt)) refit("main_sleep_only", list(), term, alt)

  # estimators that reuse the model data: each builds it inside its own worker and frees it again
  with_context <- function(id, body)
    unit(id, function() {
      b <- bundle$get()                           # built (once) here in the parent, shared with the worker by fork
      isolated(cfg, {
        cx <- sensitivity_context(b, spec_l, cfg)
        if (!enough_events(cx$ec, cfg)) not_estimable_row(id, term, cx$ec) else body(cx)
      })
    })
  with_context("conditional_logit", function(cx)
    term_row_from_table(fit_conditional_logit(cx$dat, spec_l, cx$ycol), "conditional_logit", term, cx$ec))
  with_context("cluster_robust_logit", function(cx)
    term_row_from_table(fit_cluster_robust_logit(cx$dat, cx$ftxt), "cluster_robust_logit", term, cx$ec))
  with_context("spline_vs_linear", function(cx) {
    f_spline <- sub(paste0("\\b", term, "\\b"), paste0("splines::ns(", term, ", 3)"), cx$ftxt)
    lrt_row("spline_vs_linear", term, cx$ec, lrt_ll(fit_loglik(cx$ftxt, cx$dat, cfg), fit_loglik(f_spline, cx$dat, cfg)),
            "small p = non-linear")
  })
  with_context("random_slope", function(cx) {
    f_slope <- sub("\\(1 \\| person_id\\)$", paste0("(1 | person_id) + (0 + ", term, " | person_id)"), cx$ftxt)
    lrt_row("random_slope", term, cx$ec, lrt_ll(fit_loglik(cx$ftxt, cx$dat, cfg), fit_loglik(f_slope, cx$dat, cfg), boundary = TRUE),
            "small p = effect differs by person")
  })
  status
}

# All finished sensitivity rows (all streams) + failures.
collect_sensitivity_results <- function(ckpt_dir, cfg, data_keys) {
  cu <- collect_units(ckpt_dir, "sensitivity", cfg, data_keys)
  tab <- rbindlist(Map(function(v, st) cbind(stream = st, v), cu$values, cu$units$stream), fill = TRUE)
  if (nrow(tab) == 0L) return(list(sensitivity = tab, sensitivity_all_terms = tab, sensitivity_failures = cu$failures))
  # `sensitivity`: the tested term of every analysis (as before); `sensitivity_all_terms`: all within-person terms
  list(sensitivity = tab[tested %in% TRUE][, tested := NULL][],
       sensitivity_all_terms = tab[!analysis %in% c("spline_vs_linear", "random_slope")][, tested := NULL][],   # estimates only (the other two are tests)
       sensitivity_failures = cu$failures)
}

# ---- engine check -----------------------------------------------------------------------------------

# Compare the fast default engine (lme4, nAGQ = 0) with glmmTMB on a seeded person
# sample, for the primary term. Reports the difference so the shortcut is documented.
verify_engine_agreement <- function(bundle, spec, cfg, n_persons = 1500L, which_term = cfg$primary_term) {
  cfg_s <- cfg; cfg_s$model$sample_persons <- n_persons
  fits <- lapply(c("lme4", "glmmTMB"), function(eng) {
    cfg_s$model$engine <- eng
    r <- fit_onset_spec(bundle, spec, cfg_s)
    cbind(engine = eng, r$coef[term == which_term, .(term, estimate, se)])
  })
  out <- rbindlist(fits)
  out[, diff_vs_lme4 := estimate - estimate[engine == "lme4"]]
  out[]
}
