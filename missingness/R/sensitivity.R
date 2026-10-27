# ==============================================================================
# sensitivity.R
# ------------------------------------------------------------------------------
# TRUE sensitivity analyses of the primary onset model: the SAME question
# (is a higher-than-usual prior day followed by lower odds that a gap starts?) answered
# under different assumptions. Each returns one row for the primary term, so the
# results can be read side by side and plotted as a forest plot.
#
#   main                  the reference fit
#   full_period_standardization  the person's "usual" (and the between-person reference) from the WHOLE observation period
#   trailing_baseline     "usual" = mean of the previous 14 observed days (no look-ahead)
#   exclude_partial_wear  drop prior days that look like partial wear (device on/off)
#   first_gap_only        only days up to each person's first gap (no gap history at all)
#   cci_removed / cci_as_groups / cci_within_cohort
#                         the CCI score (strongly tied to the cohorts) left out / entered as groups (reference = the most common group) /
#                         entered as the score minus the cohort mean (the part that is not cohort membership)
#   day_before_gap / *_same_rows / both_days
#                         the PREDICTOR-DAY check (see below): the primary model reads the predictors from the last FULL day d - 1; these read them from the
#                         day right before the gap (d, which may be a partial day), on all rows / on the rows of the primary model, or both days in one model
#   main_sleep_only       sleep only: a night counts as observed only if main sleep exists
#   conditional_logit     within-person comparison, no random-effect assumption
#   cluster_robust_logit  pooled logistic model, person-clustered standard errors
#   spline_vs_linear      likelihood-ratio test: is the within effect non-linear?
#   random_slope          likelihood-ratio test: does the effect differ between people?
#
# (The stricter gap definitions k = 2, 3, 7 are model D, in the onset grid.)
# ==============================================================================

# Rows of one sensitivity analysis: the tested term (`tested = TRUE`) AND every other term of the same model
# (predictors and adjustment covariates; not the intercept and not the spline basis, whose coefficients have no
# reading one by one - they are tested jointly, see term_tests), so one analysis can be drawn for all terms.
# The tested term alone is what `sensitivity.tsv` reports; the full set is `sensitivity_all_terms.tsv`.
# The multicollinearity table of the fit travels along as attr(, "gvif").
is_reportable_term <- function(term) term != "(Intercept)" & !grepl("^splines::ns", term)

extract_term_row <- function(res, analysis, which_term) {
  info <- res$info
  if (is.null(res$coef)) {
    return(data.table(analysis = analysis, term = which_term, estimate = NA_real_, se = NA_real_,
                      conf_low = NA_real_, conf_high = NA_real_, p = NA_real_,
                      or = NA_real_, or_low = NA_real_, or_high = NA_real_,
                      n_events = info$n_events, n_persons_event = info$n_persons_event,
                      status = info$status, note = info$note, tested = TRUE))
  }
  rows <- res$coef[is_reportable_term(term)]
  out <- rows[, .(analysis = analysis, term, estimate, se, conf_low, conf_high, p, or, or_low, or_high,
                  n_events, n_persons_event, status = info$status, note = info$note, tested = term == which_term)]
  if (!is.null(res$gvif)) setattr(out, "gvif", cbind(analysis = analysis, res$gvif[, .(term, df, gvif, gvif_adj, vif_equiv, flag)]))
  out
}

# Same, for an estimator that returns its own coefficient table.
term_row_from_table <- function(tab, analysis, which_term, ec, status = "ok", note = "") {
  if (nrow(tab[term == which_term]) == 0L) {
    return(data.table(analysis = analysis, term = which_term, estimate = NA_real_, se = NA_real_,
                      conf_low = NA_real_, conf_high = NA_real_, p = NA_real_, or = NA_real_,
                      or_low = NA_real_, or_high = NA_real_, n_events = ec$n_events,
                      n_persons_event = ec$n_persons_event, status = "term_missing", note = note, tested = TRUE))
  }
  rows <- add_odds_ratios(copy(tab[is_reportable_term(term)]))
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
  f <- stats::as.formula(paste0(ycol, " ~ ", paste(c(paste0("w_", preds), tv, "strata(person_id)"),
                                                   collapse = " + ")))
  # `strata()` must be written WITHOUT the "survival::" prefix: only the plain name is recognised as the stratification
  # special (with the prefix the person ids would enter as ~1000 dummy covariates and no conditional likelihood would be used)
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

# ---- the predictor-day check ---------------------------------------------------------------------------------------
#
# The onset model asks, for every day d with a record, whether a gap starts on d + 1. Day d is the LAST day with a record before the gap; the device may
# have been taken off during it, so its steps can be low for a technical reason, which creates the association "low steps -> gap tomorrow" without any
# change in behaviour. The PRIMARY model therefore reads the predictors from the last FULL day, d - 1 (the calendar day before d; as model B uses the first
# FULL day back, r + 1, and not the partial day of the return). It needs a record on d - 1: a d that is the first record or the first day after an
# earlier gap has no d - 1 and is not in the primary model. The sensitivity analyses read the predictors from the day right before the gap instead:
#
#   day_before_gap             the earlier primary model: predictors of day d, every row with complete predictors on d (a few more rows than the primary model)
#   last_full_day_same_rows /  the predictors of d - 1 / of d on the rows that have complete predictors on BOTH days (same rows and outcomes, only the
#   day_before_gap_same_rows   predictor day differs: the clean comparison)
#   both_days                  day d and the primary term on day d - 1 in ONE model (d - 1 given d), the same rows
# Reading. A lag of one more day lowers the association even for a real behavioural effect (activity on consecutive days is correlated, rho), so the
# two single-day estimates alone cannot separate "partial-day artifact" from "real effect". Three quantities do: (1) the joint model and the Wald test of
# the difference b_d - b_(d-1) (the covariance of the two estimates is in the fit); (2) the observed mean steps z-score on day d and on day d - 1 by whether
# a gap follows, and above all the CHANGE from d - 1 to d: a partial last day shows as a drop at d that is larger before a gap; (3) the ratio of the
# estimates against rho. Association, not proof: day d - 1 can itself be partial, and day d can be low for a behavioural reason (illness) that also ends
# the wear.

# Pieces of the joint model: the two estimates, their variances and covariance (for the test of their difference).
both_days_fit_summary <- function(fe, V, term) {
  lag <- paste0("lag1_", term); V <- as.matrix(V)
  if (!all(c(term, lag) %in% names(fe)) || !all(c(term, lag) %in% rownames(V)))
    return(data.table(b_cur = NA_real_, v_cur = NA_real_, b_lag = NA_real_, v_lag = NA_real_, cov_cur_lag = NA_real_))
  data.table(b_cur = unname(fe[[term]]), v_cur = V[term, term], b_lag = unname(fe[[lag]]), v_lag = V[lag, lag], cov_cur_lag = V[term, lag])
}

# Cluster-robust (person) comparison of x between rows with y = 1 and y = 0: difference of means, its standard error (CR1 sandwich of the regression
# x ~ y, closed form so that it works on millions of rows) and the two group means with their own cluster-robust intervals.
day_contrast <- function(x, y, id) {
  ok <- !is.na(x) & !is.na(y); x <- x[ok]; y <- y[ok]; id <- id[ok]
  n <- length(x); n1 <- sum(y == 1L); n0 <- n - n1
  if (n1 < 2L || n0 < 2L) return(NULL)
  m1 <- mean(x[y == 1L]); m0 <- mean(x[y == 0L]); b <- m1 - m0
  e <- x - ifelse(y == 1L, m1, m0)
  S <- rowsum(cbind(e, e * y), id)                              # cluster scores of (intercept, y)
  meat <- crossprod(S); bread <- solve(matrix(c(n, n1, n1, n1), 2L))
  g <- nrow(S); V <- bread %*% meat %*% bread * (g / (g - 1)) * ((n - 1) / (n - 2))
  se <- sqrt(V[2L, 2L])
  ci1 <- cluster_mean_ci(x[y == 1L], id[y == 1L]); ci0 <- cluster_mean_ci(x[y == 0L], id[y == 0L])
  list(mean_gap = m1, mean_no_gap = m0, difference = b, se = se, low = b - 1.96 * se, high = b + 1.96 * se,
       p = 2 * stats::pnorm(abs(b / se), lower.tail = FALSE), ci_gap = ci1, ci_no_gap = ci0)
}

# The observed pattern behind the check, on the rows of the primary model (identical in every mode): the mean within-person steps z-score of the day
# right before the gap (d) and of the last full day (d - 1) by whether a gap starts on d + 1, and their change. Returns list(pattern, descriptive) with a
# `stream` column.
day_describe <- function(bundle, spec, cfg, term = cfg$primary_term) {
  spec <- as.list(spec); spec$predictors <- list(spec$predictors[[1L]])
  dat <- build_onset_data(bundle, spec, cfg, list(day = "both"))
  ycol <- attr(dat, "ycol"); lagc <- paste0("lag1_", term)
  d <- data.table(person_id = dat$person_id, y = dat[[ycol]], z0 = dat[[term]], z1 = dat[[lagc]]); rm(dat)
  d[, dz := z0 - z1]
  s <- bundle$stream; rho <- suppressWarnings(stats::cor(d$z0, d$z1))
  qs <- list(z1 = DAY_LABELS[["d1"]], z0 = DAY_LABELS[["d0"]], dz = DAY_LABELS[["dz"]])
  pat <- list(); des <- list()
  for (q in names(qs)) {
    r <- day_contrast(d[[q]], d$y, d$person_id)
    if (is.null(r)) next
    des[[q]] <- data.table(stream = s, quantity = qs[[q]], mean_gap = r$mean_gap, mean_no_gap = r$mean_no_gap, difference = r$difference, se = r$se,
                           low = r$low, high = r$high, p = r$p, n_rows = nrow(d), n_events = sum(d$y == 1L), n_people = uniqueN(d$person_id), rho_d_dminus1 = rho)
    if (q != "dz") for (g in c("gap", "no_gap")) {
      ci <- if (g == "gap") r$ci_gap else r$ci_no_gap
      pat[[paste(q, g)]] <- data.table(stream = s, day = qs[[q]], group = if (g == "gap") "A gap starts the next day" else "No gap the next day",
                                       mean_z = ci$mean, low = ci$lo, high = ci$hi, n_rows = ci$n_gaps, n_people = ci$n_persons)
    }
  }
  list(pattern = rbindlist(pat), descriptive = rbindlist(des))
}

# Labels of the two days (shared by the table and the figures).
DAY_LABELS <- c(d1 = "Last full day (d - 1)", d0 = "Day before the gap (d)", dz = "Change from d - 1 to d")
DAY_MODEL_LEVELS <- c("Primary model: last full day d - 1", "Last full day d - 1, rows with both days", "Day before the gap d, rows with both days",
                      "Day before the gap d, all rows", "Joint model: last full day d - 1", "Joint model: day before the gap d", "Difference d minus d - 1 (joint model)")

# The table of the check: the primary model, the day-before-gap fits, and the formal difference of the joint model, per stream.
#   sens: the stored sensitivity rows (all analyses); info: the `info` rows of the joint fits (both_days_fit_summary + loglik...)
day_check_table <- function(sens, info, desc, cfg) {
  term <- cfg$primary_term; lag <- paste0("lag1_", term)
  if (!nz(sens) || !any(sens$analysis %in% c("day_before_gap", "day_before_gap_same_rows", "last_full_day_same_rows", "both_days"))) return(data.table())
  take <- function(analysis_, term_, label) {
    x <- sens[analysis == analysis_ & term == term_ & !is.na(estimate)]
    if (!nrow(x)) return(NULL)
    x[, .(stream, model = label, term = term_, estimate, se, or, or_low, or_high, p, n_events, n_persons_event, status)]
  }
  out <- rbindlist(list(
    take("main", term, DAY_MODEL_LEVELS[1L]),
    take("last_full_day_same_rows", term, DAY_MODEL_LEVELS[2L]),
    take("day_before_gap_same_rows", term, DAY_MODEL_LEVELS[3L]),
    take("day_before_gap", term, DAY_MODEL_LEVELS[4L]),
    take("both_days", lag, DAY_MODEL_LEVELS[5L]),
    take("both_days", term, DAY_MODEL_LEVELS[6L])), fill = TRUE)
  if (nrow(out) == 0L) return(out)
  out[, model := as.character(model)]
  if (nz(info) && "b_cur" %in% names(info)) {
    ct <- info[!is.na(b_cur), {
      dlt <- b_cur - b_lag; se <- sqrt(pmax(v_cur + v_lag - 2 * cov_cur_lag, 0))
      .(model = DAY_MODEL_LEVELS[7L], term = "b_d - b_(d-1)", estimate = dlt, se = se, or = exp(dlt), or_low = exp(dlt - 1.96 * se), or_high = exp(dlt + 1.96 * se),
        p = 2 * stats::pnorm(abs(dlt / se), lower.tail = FALSE), n_events = n_events, n_persons_event = n_persons_event, status = status)
    }, by = stream]
    if (nrow(ct)) {
      ct[, p_holm := stats::p.adjust(p, "holm")]
      lr <- merge(out[model == DAY_MODEL_LEVELS[3L], .(stream, b_day_d_alone = estimate)], out[model == DAY_MODEL_LEVELS[2L], .(stream, b_full_day_alone = estimate)], by = "stream")
      lr[, lagged_over_current := b_full_day_alone / b_day_d_alone]
      ct <- merge(ct, lr[, .(stream, lagged_over_current)], by = "stream", all.x = TRUE)
      if (nz(desc)) ct <- merge(ct, unique(desc[, .(stream, rho_d_dminus1)]), by = "stream", all.x = TRUE)
      out <- rbind(out, ct, fill = TRUE)
    }
  }
  out[, model := factor(model, levels = DAY_MODEL_LEVELS)]
  out[order(stream, as.integer(model))][]
}

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
SENSITIVITY_SCHEMA <- 3L

# The sensitivity analyses of the PRIMARY spec of one stream, each an independent checkpointed unit
# (so a failed or slow analysis never costs the others).
#   bundle : lazy(prepare_stream_bundle(..., trailing = TRUE))   alt : lazy bundle of the main-sleep-only
#   rule (sleep only; NULL otherwise)      prim : the primary row of the onset grid
#   alt_full : lazy bundle with standardization = "full_period" (NULL: that analysis is not run)
# Returns the run_unit() status of each analysis.
run_sensitivity_stage <- function(s, bundle, alt, prim, cfg, ckpt_dir, dkey, term = cfg$primary_term, alt_full = NULL) {
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
  # CCI is strongly tied to the cohorts (the score counts diabetes and cancer): the primary model uses the score; these three refits show whether the
  # within-person estimates depend on how (or whether) it is adjusted for
  if (isTRUE(cfg$model$adjust_cci)) {
    refit("cci_removed", list(cci_mode = "none"))
    refit("cci_as_groups", list(cci_mode = "groups"))
    refit("cci_within_cohort", list(cci_mode = "within_cohort"))
  }
  if (s == "sleep" && !is.null(alt)) refit("main_sleep_only", list(), term, alt)
  # the predictor-day check (see above): the primary model reads day d - 1; here the predictors of the day right before the gap (d), and both days in one model
  refit("day_before_gap", list(day = "day_before_gap"))
  refit("last_full_day_same_rows", list(day = "last_full_day", same_rows = TRUE))
  refit("day_before_gap_same_rows", list(day = "day_before_gap", same_rows = TRUE))
  unit("both_days", function() {
    r <- fit_onset_spec(bundle$get(), prim, cfg, list(day = "both"))
    out <- extract_term_row(r, "both_days", term)
    setattr(out, "day_info", r$info)                        # the two estimates, their variances and covariance: the formal difference
    out
  })
  unit("day_description", function() {
    b <- bundle$get()
    isolated(cfg, { x <- day_describe(b, prim, cfg); setattr(x, "kind", "day_description"); x })
  })

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
  # same model, but the person's usual (and the between-person reference) from the WHOLE observation period (look-ahead).
  # Run last: the other bundles are released first, so two person-day tables are never alive together.
  if (!is.null(alt_full)) {
    bundle$release(); if (!is.null(alt)) alt$release()
    refit("full_period_standardization", list(), term, alt_full)
  }
  status
}

# All finished sensitivity rows (all streams) + failures.
collect_sensitivity_results <- function(ckpt_dir, cfg, data_keys) {
  cu <- collect_units(ckpt_dir, "sensitivity", cfg, data_keys)
  is_desc <- vapply(cu$values, function(v) identical(attr(v, "kind"), "day_description"), logical(1))     # the predictor-day pattern is not a row of the model table
  desc_vals <- cu$values[is_desc]; desc_st <- cu$units$stream[is_desc]
  cu$values <- cu$values[!is_desc]; cu$units <- cu$units[!is_desc]
  tab <- rbindlist(Map(function(v, st) cbind(stream = st, v), cu$values, cu$units$stream), fill = TRUE)
  gv  <- rbindlist(Map(function(v, st) { g <- attr(v, "gvif"); if (is.null(g)) NULL else cbind(stream = st, g) }, cu$values, cu$units$stream), fill = TRUE)
  fd_info <- rbindlist(Map(function(v, st) { g <- attr(v, "day_info"); if (is.null(g)) NULL else cbind(stream = st, g) }, cu$values, cu$units$stream), fill = TRUE)
  fd_pat  <- rbindlist(lapply(desc_vals, `[[`, "pattern"), fill = TRUE); fd_des <- rbindlist(lapply(desc_vals, `[[`, "descriptive"), fill = TRUE)
  if (nrow(tab) == 0L) return(list(sensitivity = tab, sensitivity_all_terms = tab, sensitivity_collinearity = gv, sensitivity_failures = cu$failures))
  # `sensitivity`: the tested term of every analysis (as before); `sensitivity_all_terms`: all within-person terms
  list(sensitivity = tab[tested %in% TRUE][, tested := NULL][],
       sensitivity_all_terms = tab[!analysis %in% c("spline_vs_linear", "random_slope")][, tested := NULL][],   # estimates only (the other two are tests)
       sensitivity_collinearity = gv, sensitivity_failures = cu$failures,
       day_check = day_check_table(tab, fd_info, fd_des, cfg), day_pattern = fd_pat, day_descriptive = fd_des)
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
