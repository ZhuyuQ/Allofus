# ==============================================================================
# onset_model.R
# ------------------------------------------------------------------------------
# Mixed-effects logistic model of GAP ONSET (the primary analysis, plus models A and D).
#
#   One row  = one observed day d (the "prior day"), except each person's last day.
#   Outcome  = 1 if a gap starts on day d+1 (y_ge<k>: ... and lasts >= k days).
#   Predictors = the prior day's activity / sleep / heart-rate values, split into
#                w_  (how unusual for THIS person) and b_ (how high/low the person is).
#   (1 | person_id) gives each person their own baseline tendency to have gaps.
#
# Role of each model (see the design notes):
#   primary : k = 1, joint predictors, adjusted -> test of cfg$primary_term
#   model A : gap-history covariates inside the adjusted formula
#   model D : the same joint adjusted model with k = 2, 3, 7
# ==============================================================================

# ---- specification grid ------------------------------------------------------------

# One row per model to fit. `predictors` is a list column of predictor names.
make_onset_grid <- function(cfg) {
  rows <- list()
  add <- function(stream, set_name, preds, k, adjusted, role) {
    rows[[length(rows) + 1L]] <<- data.table(
      stream = stream, set_name = set_name, predictors = list(preds),
      k = k, adjusted = adjusted, role = role)
  }
  act <- cfg$activity_predictors
  for (stream in names(cfg$own_predictors)) {
    own <- cfg$own_predictors[[stream]]
    sets <- if (stream == "activity") {
      c(setNames(as.list(act), paste0("uni_", act)), list(joint = act))
    } else {
      c(setNames(as.list(own), paste0("own_uni_", own)), list(own_multi = own),
        setNames(as.list(act), paste0("act_uni_", act)), list(act_multi = act),
        list(joint = c(own, act)))
    }
    for (nm in names(sets)) {
      for (adj in c(FALSE, TRUE)) {
        role <- if (nm == "joint" && adj) "primary" else "exploratory"
        add(stream, nm, sets[[nm]], cfg$onset$primary_k, adj, role)
      }
    }
    # model D: stricter gap definitions, joint adjusted model only
    for (k in cfg$onset$secondary_k) add(stream, "joint", sets$joint, k, TRUE, "secondary_threshold")
  }
  g <- rbindlist(rows)
  g[, spec_id := sprintf("%s|%s|k%d|%s", stream, set_name, k, fifelse(adjusted, "adj", "unadj"))]
  setcolorder(g, "spec_id")
  g[]
}

# Adjustment terms: formula term -> data column it depends on (so constant columns
# can be dropped automatically, e.g. `sex2` inside the all-male Prostate cohort).
# Calendar year is NOT a covariate. Device family (replaces it), CCI score and its missing flag are; a term whose column is not in the
# data (switched off in cfg$model, or removed by a sensitivity variant) is dropped like a constant one.
adjuster_terms <- function() {
  data.table(
    term = c("cohort", "age10", "sex2", "dow_t", "splines::ns(study_day, 3)", "device", "cci", "cci_missing",
             "prev_gap_len_l", "any_prior_gap", "days_since_return_l"),
    col  = c("cohort", "age10", "sex2", "dow_t", "study_day", "device", "cci", "cci_missing",
             "prev_gap_len_l", "any_prior_gap", "days_since_return_l"))
}

# Build the model formula text.
#   baseline : "whole" -> w_ terms, "trailing" -> tw_ terms
#   drop_cols: adjuster columns to leave out (constant in the data)
build_onset_formula <- function(preds, adjusted, k, baseline = "whole", drop_cols = character(),
                                random = "(1 | person_id)", extra_terms = character()) {
  pre   <- if (baseline == "trailing") "tw_" else "w_"
  terms <- c(paste0(pre, preds), paste0("b_", preds))
  if (adjusted) {
    adj <- adjuster_terms()[!col %in% drop_cols]
    terms <- c(terms, adj$term)
  }
  terms <- c(terms, extra_terms)
  paste0("y_ge", k, " ~ ", paste(terms, collapse = " + "), if (nzchar(random)) paste0(" + ", random))
}

# ---- model data -----------------------------------------------------------------------

# Select the rows and columns for one specification.
#
# variant fields (all optional):
#   baseline      "whole" | "trailing"
#   risk_set      "common" (default): complete cases over ALL predictors of the stream, so
#                 univariate / multivariate / joint models use identical rows;
#                 "spec": complete cases over the spec's own predictors only
#   drop_partial  TRUE -> drop days flagged as partial wear (own stream or activity)
#   first_gap     TRUE -> keep only days up to each person's first gap onset
#   cohort        a cohort label -> fit within that cohort only
#   runin_days    N -> the run-in analysis (runin.R): only prior days AFTER the first N days of follow-up (study_day > N) of people with
#                 a complete run-in; adds cov0 (log2(1 + missing days in the run-in), centred), cov0_lin and cov0_cat
#   day           WHICH DAY the predictors are read from (default cfg$onset$predictor_day = "last_full_day"). Let d be the last day with a record before a
#                 possible gap (the gap starts on d + 1). The device may be taken off during d, so d can be a partial day.
#                   "last_full_day"  (primary): the predictors of the calendar day d - 1, the last FULL day. Rows whose day d - 1 has no record (d is
#                                    the first record, or the first day after an earlier gap) or has incomplete predictors are not in the model.
#                   "day_before_gap" (sensitivity): the predictors of day d itself (the primary model of the earlier versions)
#                   "both":          day d, plus the primary term on day d - 1 as `lag1_<term>` in the same model (rows with both days)
#   same_rows     TRUE -> only the rows that have complete predictors on BOTH days d and d - 1 (the rows of the "both" model): with day = "last_full_day"
#                 or "day_before_gap" the clean comparison of the two days (same rows and outcomes, only the predictor day differs)
#   runin_terms   which run-in terms the formula gets: "reference" | "main" | "interaction" | "categorical" | "linear_interaction"
build_onset_data <- function(bundle, spec, cfg, variant = list()) {
  baseline <- variant$baseline %||% "whole"
  pre   <- if (baseline == "trailing") "tw_" else "w_"
  preds <- spec$predictors[[1L]]
  ycol  <- paste0("y_ge", spec$k)
  rs    <- bundle$rs

  if (baseline == "trailing" && !isTRUE(bundle$trailing)) {
    stop("this bundle has no trailing baselines: build it with prepare_stream_bundle(..., trailing = TRUE)", call. = FALSE)
  }
  all_preds <- unique(c(cfg$own_predictors[[bundle$stream]], cfg$activity_predictors))
  cc_preds  <- if ((variant$risk_set %||% "common") == "common") all_preds else preds
  cc_cols   <- c(paste0(pre, cc_preds), paste0("b_", cc_preds))
  if (spec$adjusted) cc_cols <- c(cc_cols, "age10", "sex2", if (is.null(variant$cohort)) "cohort")

  # rows to keep, decided with logical vectors (no copy of the big table)
  day <- variant$day %||% cfg$onset$predictor_day %||% "last_full_day"
  if (!day %in% c("last_full_day", "day_before_gap", "both")) stop("unknown predictor day: ", day, call. = FALSE)
  use_lag  <- day != "day_before_gap" || isTRUE(variant$same_rows)       # is the row of the calendar day d - 1 needed?
  cc_d <- complete.cases(rs[, ..cc_cols])                                  # complete predictors on the row's own day d
  base <- rep(TRUE, nrow(rs))
  if (isTRUE(variant$first_gap)) base <- base & first_gap_rows(rs)
  if (!is.null(variant$cohort))  base <- base & (rs$cohort == variant$cohort)
  lag_idx <- NULL
  if (use_lag) {
    # row of the calendar day BEFORE each row's day d: NA when that day has no record (d is the first record, or the first day after a gap)
    lag_idx <- rep(NA_integer_, nrow(rs)); b_ix <- which(base)
    lag_idx[b_ix] <- rs[.(rs$person_id[b_ix], rs$date[b_ix] - 1L), which = TRUE, on = c("person_id", "date")]
    rm(b_ix)
  }
  cc_l <- if (use_lag) cc_d[lag_idx] %in% TRUE                              # complete predictors on day d - 1
  keep <- base & switch(day, last_full_day = if (isTRUE(variant$same_rows)) cc_d & cc_l else cc_l, day_before_gap = if (use_lag) cc_d & cc_l else cc_d, both = cc_d & cc_l)
  rm(base, cc_d)
  if (isTRUE(variant$drop_partial)) {                                      # a partial day among the days the predictors are read from
    pw_rows <- function(ix) (rs$pw_own[ix] %in% TRUE) | (rs$pw_act[ix] %in% TRUE)
    if (day != "last_full_day") keep <- keep & !pw_rows(seq_len(nrow(rs)))
    if (use_lag)                keep <- keep & !pw_rows(lag_idx)
  }
  if (!is.null(variant$runin_days)) {                                     # only days AFTER the run-in: the days the predictors are read from
    sd_pred <- if (day == "last_full_day" || use_lag) rs$study_day[lag_idx] else rs$study_day
    keep <- keep & (sd_pred > variant$runin_days) & !is.na(rs$runin_missing); rm(sd_pred)
  }

  # only the columns this model uses (the bundle has ~30; copying all of them dominated memory)
  # CCI variant (sensitivity): "score" (default, per +1 point + missing flag), "none", "groups" (cci_grp), "within_cohort" (score minus the cohort mean)
  cci_mode <- variant$cci_mode %||% "score"
  adj_cols <- c("cohort", "age10", "sex2", "dow_t", "study_day", "device", "prev_gap_len", "any_prior_gap", "days_since_return",
                switch(cci_mode, score = c("cci", "cci_missing"), none = character(), groups = "cci_grp", within_cohort = c("cci_wc", "cci_missing")))
  adj_cols <- intersect(adj_cols, names(rs))                               # device / cci columns exist only when cfg$model switches them on
  pred_cols <- c(paste0(pre, preds), paste0("b_", preds))
  out_cols <- unique(c("person_id", ycol, pred_cols,
                       if (spec$adjusted) adj_cols, if (!is.null(variant$runin_days)) "runin_missing"))
  out <- rs[keep, ..out_cols]
  # the adjusters (weekday, time in study, gap history) stay those of day d, the day before the gap: weekday of d is one-to-one with the weekday of d - 1,
  # and the gap history of d is that of d - 1 (no gap lies between them)
  if (day == "last_full_day") {
    li <- lag_idx[keep]
    for (cl in pred_cols) set(out, j = cl, value = rs[[cl]][li])
  } else if (day == "both") {
    li <- lag_idx[keep]
    set(out, j = paste0("lag1_", cfg$primary_term), value = rs[[cfg$primary_term]][li])
  }
  rm(lag_idx, keep, cc_l)

  if (!is.null(cfg$model$sample_persons)) {
    ids <- sample_person_ids(out$person_id, cfg$model$sample_persons, cfg$seed)
    out <- out[person_id %in% ids]
  }
  if (spec$adjusted) {
    # day counts on the log2 scale, uncapped: a coefficient is "per doubling" of (1 + days)
    out[, `:=`(prev_gap_len_l = log2(1 + prev_gap_len), days_since_return_l = log2(1 + days_since_return))]
    out[, c("prev_gap_len", "days_since_return") := NULL]
    if (cci_mode == "within_cohort") { out[, cci := cci_wc]; out[, cci_wc := NULL] }
  }
  if (!is.null(variant$runin_days)) add_runin_columns(out, variant$runin_days, cfg)     # after the person sample, so the centring refers to the people analysed
  attr(out, "ycol") <- ycol
  out[]
}

# The values of `cols` of the risk-set table on the day the model reads its predictors from (cfg$onset$predictor_day), row by row: the row's own day d for
# "day_before_gap", the calendar day d - 1 for "last_full_day" (NA when that day has no record). For the descriptive tables that must describe the model's rows.
rs_at_predictor_day <- function(rs, cols, cfg) {
  if (!identical(cfg$onset$predictor_day %||% "last_full_day", "last_full_day")) return(rs[, ..cols])
  li <- rs[.(rs$person_id, rs$date - 1L), which = TRUE, on = c("person_id", "date")]
  rs[li, ..cols]
}

# Are there enough events to fit? (reported, not an error)
event_counts <- function(dat, ycol) {
  y <- dat[[ycol]]
  list(n_rows = length(y), n_persons = uniqueN(dat$person_id), n_events = sum(y),
       n_persons_event = uniqueN(dat$person_id[y == 1L]))
}
enough_events <- function(ec, cfg) {
  ec$n_events >= cfg$model$min_events && ec$n_persons_event >= cfg$model$min_persons_with_events
}

# Columns of the adjusters that are constant in `dat` (to be dropped from the formula).
constant_adjuster_cols <- function(dat, adjusted) {
  if (!adjusted) return(character())
  cols <- adjuster_terms()$col
  absent <- cols[!cols %in% names(dat)]                       # switched off or removed by a sensitivity variant
  present <- intersect(cols, names(dat))
  c(absent, present[vapply(present, function(cl) uniqueN(dat[[cl]]) < 2L, logical(1))])
}

# ---- fitting ----------------------------------------------------------------------------

# Fit a binomial mixed model with the configured engine.
fit_binary_mixed <- function(formula_text, dat, cfg) {
  f <- stats::as.formula(formula_text)
  if (cfg$model$engine == "lme4") {
    lme4::glmer(f, data = dat, family = stats::binomial(), nAGQ = cfg$model$nAGQ,
                control = lme4::glmerControl(optimizer = "bobyqa",
                                             optCtrl = list(maxfun = 2e5)))
  } else {
    glmmTMB::glmmTMB(f, data = dat, family = stats::binomial())
  }
}

# Fit a model, keep only its log-likelihood and degrees of freedom, and free the model at once.
# Comparing two models this way needs the memory of ONE fit, not two.
fit_loglik <- function(formula_text, dat, cfg) {
  fit <- fit_binary_mixed(formula_text, dat, cfg)
  ll <- stats::logLik(fit)
  out <- list(ll = as.numeric(ll), df = attr(ll, "df"))
  rm(fit, ll); invisible(gc())
  out
}

# Did an lme4 optimiser converge? lme4 stores the "boundary (singular) fit" message next to real
# convergence problems; a singular fit (person variance ~ 0) did converge and is reported separately
# (`singular`), so only the other messages count here.
lme4_converged <- function(fit) {
  msgs <- fit@optinfo$conv$lme4$messages
  isTRUE(fit@optinfo$conv$opt == 0) && !any(!grepl("singular", msgs, fixed = TRUE))
}

# Did the optimiser converge? (a model that did not should not be interpreted)
fit_converged <- function(fit) {
  if (inherits(fit, "glmmTMB")) isTRUE(fit$fit$convergence == 0) && isTRUE(fit$sdr$pdHess)
  else lme4_converged(fit)
}

# Fixed-effect table with Wald 95% intervals.
tidy_fixed <- function(fit) {
  cf <- if (inherits(fit, "glmmTMB")) summary(fit)$coefficients$cond else stats::coef(summary(fit))
  out <- data.table(term = rownames(cf), estimate = cf[, 1], se = cf[, 2], p = cf[, 4])
  out[, `:=`(conf_low = estimate - 1.96 * se, conf_high = estimate + 1.96 * se)]
  out[]
}

# Variance of the person random intercept (latent scale).
person_variance <- function(fit) {
  vc <- if (inherits(fit, "glmmTMB")) glmmTMB::VarCorr(fit)$cond$person_id else lme4::VarCorr(fit)$person_id
  as.numeric(vc[1L, 1L])
}

# Median odds ratio: the typical odds ratio between two random people.
median_odds_ratio <- function(sigma2) exp(sqrt(2 * sigma2) * stats::qnorm(0.75))

# Add odds ratios, all in the NORMAL direction (a higher value of the predictor is a positive step; OR = exp(estimate)):
#   w_ / tw_ terms -> OR per +1 SD HIGHER than the person's usual (OR < 1: higher activity, fewer gaps)
#   b_ terms       -> OR per +1 SD higher person-level mean
#   other          -> OR per +1 unit (or against the reference level for factors)
add_odds_ratios <- function(tab) {
  tab[, `:=`(or      = exp(estimate), or_low = exp(conf_low), or_high = exp(conf_high),
             or_unit = fifelse(grepl("^(lag1_)?(t?w|b)_", term), "per +1 SD", "per +1 unit"))]
  tab[]
}

# Fit one specification and return a tidy result (never throws for a model failure:
# the failure is recorded in `status`).
#   status: "ok" | "not_estimable" (too few events) | "not_converged" | "failed"
fit_onset_spec_inner <- function(bundle, spec, cfg, variant = list(), keep_fit = FALSE) {
  spec <- as.list(spec); spec$predictors <- list(spec$predictors[[1L]])
  dat  <- build_onset_data(bundle, spec, cfg, variant)
  ycol <- attr(dat, "ycol")
  ec   <- event_counts(dat, ycol)
  meta <- data.table(
    spec_id = spec$spec_id %||% NA_character_, stream = spec$stream, set_name = spec$set_name,
    k = spec$k, adjusted = spec$adjusted, role = spec$role %||% NA_character_,
    baseline = variant$baseline %||% "whole", cohort_stratum = variant$cohort %||% "all",
    standardization = bundle$standardization %||% cfg$standardization$method,
    n_rows = ec$n_rows, n_persons = ec$n_persons, n_events = ec$n_events,
    n_persons_event = ec$n_persons_event)
  if (!is.null(variant$runin_days)) meta[, `:=`(runin_variant = variant$runin_terms %||% "reference", runin_days = variant$runin_days)]
  meta[, predictor_day := variant$day %||% cfg$onset$predictor_day %||% "last_full_day"]
  if (isTRUE(variant$same_rows)) meta[, predictor_day := paste0(predictor_day, "_same_rows")]

  if (!enough_events(ec, cfg)) {
    return(list(coef = NULL, info = cbind(meta, status = "not_estimable", sigma2 = NA_real_,
                                          mor = NA_real_, note = "too few events"), fit = NULL))
  }
  drop_cols <- constant_adjuster_cols(dat, spec$adjusted)
  extra <- c(if (spec$adjusted && identical(variant$cci_mode, "groups") && "cci_grp" %in% names(dat) && uniqueN(dat$cci_grp) > 1L) "cci_grp",
             if (!is.null(variant$runin_days)) runin_formula_terms(variant$runin_terms %||% "reference", cfg$primary_term, dat),
             if (identical(variant$day, "both") && paste0("lag1_", cfg$primary_term) %in% names(dat)) paste0("lag1_", cfg$primary_term))
  ftxt <- build_onset_formula(spec$predictors[[1L]], spec$adjusted, spec$k,
                              baseline = variant$baseline %||% "whole", drop_cols = drop_cols, extra_terms = extra)

  notes <- character()
  fit <- tryCatch(
    withCallingHandlers(suppressMessages(fit_binary_mixed(ftxt, dat, cfg)),     # lme4's "boundary (singular)" message is reported in `singular`
                        warning = function(w) { notes <<- c(notes, conditionMessage(w))
                                                invokeRestart("muffleWarning") }),
    error = function(e) e)
  if (inherits(fit, "error")) {
    return(list(coef = NULL, info = cbind(meta, status = "failed", sigma2 = NA_real_,
                                          mor = NA_real_, note = conditionMessage(fit)), fit = NULL))
  }
  ok <- fit_converged(fit)
  s2 <- person_variance(fit)
  coef <- add_odds_ratios(tidy_fixed(fit))
  coef <- cbind(meta[rep(1L, nrow(coef))], coef)
  coef[, status := if (ok) "ok" else "not_converged"]
  # covariance of the intercept and the primary term: needed for the interval of the absolute-risk curve
  V <- if (inherits(fit, "glmmTMB")) stats::vcov(fit)$cond else as.matrix(stats::vcov(fit))
  cov_b0_term <- if (all(c("(Intercept)", cfg$primary_term) %in% rownames(V))) V["(Intercept)", cfg$primary_term] else NA_real_
  sing <- if (inherits(fit, "glmerMod")) lme4::isSingular(fit) else NA
  info <- cbind(meta, status = if (ok) "ok" else "not_converged", sigma2 = s2,
                mor = median_odds_ratio(s2), singular = sing, cov_b0_term = cov_b0_term,
                note = paste(c(if (isTRUE(sing)) "singular fit (person variance ~ 0)", utils::head(unique(notes), 3L)), collapse = " | "))
  # multicollinearity (GVIF per term) and joint Wald tests per term, from the fixed-effect covariance of this very fit
  dg <- tryCatch(fixed_effect_diagnostics(stats::setNames(coef$estimate, coef$term), V, ftxt), error = function(e) NULL)
  id_cols <- meta[1L, c("spec_id", "stream", "set_name", "k", "adjusted", "role", "baseline", "cohort_stratum", "standardization", intersect(c("runin_variant", "runin_days", "predictor_day"), names(meta))), with = FALSE]
  tag_meta <- function(tab) if (is.null(tab) || nrow(tab) == 0L) NULL else cbind(id_cols[rep(1L, nrow(tab))], tab)
  if (!is.null(variant$runin_days)) info <- cbind(info, runin_fit_summary(fit, V, ftxt, cfg$primary_term, attr(dat, "cov0_center")))
  if (identical(variant$day, "both")) info <- cbind(info, both_days_fit_summary(stats::setNames(coef$estimate, coef$term), V, cfg$primary_term))
  res <- list(coef = coef, info = info, fit = if (keep_fit) fit else NULL,
              gvif = if (!is.null(dg)) tag_meta(dg$gvif), wald = if (!is.null(dg)) tag_meta(dg$wald))
  rm(dat, fit); invisible(gc())          # give the memory back before the next model
  res
}

# Rough memory a fit needs (GB): rows x model columns x bytes per cell (cfg$resources$bytes_per_cell,
# measured for lme4::glmer on this model, including its working copies).
estimate_fit_memory_gb <- function(n_rows, n_model_cols, bytes_per_cell = 56) {
  n_rows * n_model_cols * bytes_per_cell / 1e9
}

# The person sample that applies to a spec (NULL = everybody): the exploratory sample for
# exploratory specs when one is set, otherwise cfg$model$sample_persons.
effective_sample <- function(cfg, role) {
  if (identical(role, "exploratory") && !is.null(cfg$model$exploratory_sample_persons))
    cfg$model$exploratory_sample_persons else cfg$model$sample_persons
}

# Stop ONE model, with a diagnosis, when it cannot fit in the memory that is free (the unit is then
# recorded as failed and the other models still run). Nothing is subsampled automatically; the advice
# names the explicit, documented alternative. Switch off with cfg$resources$memory_guard <- FALSE.
check_fit_memory <- function(bundle, spec, cfg, variant = list(), extra_cols = 0, safety = 1) {
  if (!isTRUE(cfg$resources$memory_guard)) return(invisible(TRUE))
  n_pred <- length(spec$predictors[[1L]])
  # adjusters: cohort 3, age 1, sex 1, weekday 6, study-time spline 3, gap history 3, device family (levels - 1), CCI 2 (an upper bound)
  n_adj <- 17 + (if (!is.null(variant$runin_days)) 4L else 0L) + (if ("device" %in% names(bundle$rs)) max(0L, nlevels(bundle$rs$device) - 1L) else 0L) + (if ("cci" %in% names(bundle$rs)) 2L else 0L)
  n_cols <- 2 * n_pred + (if (isTRUE(spec$adjusted)) n_adj else 0) + 3 + extra_cols
  n_rows <- if (!is.null(variant$cohort)) sum(bundle$rs$cohort == variant$cohort, na.rm = TRUE) else nrow(bundle$rs)
  if (!is.null(variant$runin_days)) n_rows <- sum(bundle$rs$study_day > variant$runin_days, na.rm = TRUE)       # the run-in analysis only models the days after the run-in
  samp   <- effective_sample(cfg, spec$role)
  if (!is.null(samp) && !is.null(bundle$n_persons)) n_rows <- n_rows * min(1, samp / bundle$n_persons)
  need <- safety * estimate_fit_memory_gb(n_rows, n_cols, cfg$resources$bytes_per_cell)
  headroom <- memory_headroom_gb(cfg)
  if (need > headroom && need <= cfg$resources$memory_budget_gb) headroom <- wait_for_headroom(cfg, need)
  if (need > headroom) {
    keep <- max(500L, floor((bundle$n_persons %||% 5000L) * headroom / need * 0.7))
    stop(sprintf(paste0("this model needs about %.0f GB but only %.0f GB can be used (free memory minus the %.0f GB reserve, ",
                        "budget %.0f GB). Nothing was subsampled. Options: restart the kernel and run only the setup and ",
                        "load cells; or set cfg$model$sample_persons <- %d to fit a seeded sample of people (same model, ",
                        "wider intervals); or cfg$resources$memory_guard <- FALSE to try anyway."),
                 need, headroom, cfg$resources$memory_reserve_gb, cfg$resources$memory_budget_gb, keep), call. = FALSE)
  }
  invisible(TRUE)
}

# Public entry point: the fit runs in a worker process (see run_isolated) unless the model
# object itself is wanted (keep_fit = TRUE).
fit_onset_spec <- function(bundle, spec, cfg, variant = list(), keep_fit = FALSE) {
  check_fit_memory(bundle, as.list(spec), cfg, variant)
  if (keep_fit) return(fit_onset_spec_inner(bundle, spec, cfg, variant, keep_fit = TRUE))
  isolated(cfg, fit_onset_spec_inner(bundle, spec, cfg, variant, FALSE))
}

# ---- reading the result in absolute terms ----------------------------------------------------

# Marginal (population-averaged) gap risk on the next day, integrating over the
# person random intercept. `coef` is a fitted coefficient table, `term` the within
# term; all other covariates are held at 0 / reference.
# at = within-person z-scores, e.g. c(-1, 0, 1) = -1 SD, 0 (usual), +1 SD (higher z = more steps).
marginal_risk <- function(coef, term, sigma2, at = c(-1, 0, 1), n_grid = 201L,
                          cov_b0_term = NA_real_, n_draw = 2000L, seed = 1L) {
  which_term <- term                                    # `term` is also a column of `coef`
  b0 <- coef[term == "(Intercept)", estimate]; s0 <- coef[term == "(Intercept)", se]
  b1 <- coef[term == which_term, estimate];    s1 <- coef[term == which_term, se]
  u  <- stats::qnorm(stats::ppoints(n_grid)) * sqrt(sigma2)
  risk <- function(b0, b1, a) mean(stats::plogis(b0 + b1 * a + u))
  out <- data.table(within_z = at, risk = vapply(at, function(a) risk(b0, b1, a), numeric(1)),
                    risk_low = NA_real_, risk_high = NA_real_)
  if (is.finite(cov_b0_term)) {                         # interval from the uncertainty of (b0, b1); sigma2 held fixed
    V <- matrix(c(s0^2, cov_b0_term, cov_b0_term, s1^2), 2L)
    L <- tryCatch(chol(V), error = function(e) NULL)
    if (!is.null(L)) {
      draws <- with_seed(seed, matrix(stats::rnorm(2L * n_draw), n_draw) %*% L)
      draws <- sweep(draws, 2L, c(b0, b1), "+")
      for (i in seq_along(at)) {
        r <- vapply(seq_len(n_draw), function(j) risk(draws[j, 1L], draws[j, 2L], at[i]), numeric(1))
        q <- stats::quantile(r, c(0.025, 0.975), names = FALSE)
        set(out, i, c("risk_low", "risk_high"), as.list(q))
      }
    }
  }
  out[]
}

# ---- cohort differences in the primary effect ----------------------------------------------------

# Likelihood-ratio test: does the within-person effect (`term`) differ between the groups of `by`?
#   by = "cohort"  (default): the model with vs without the term x cohort interaction; cohort is already an adjuster.
#   by = any other factor column of the model data: its main effect is added to both models, so the null model is "the same
#   effect in every group, groups may differ in their overall gap risk".
# One chi-square test on df = (number of groups - 1); it is the formal comparison of the groups, not a comparison of
# which group's own estimate is significant.
test_group_interaction_inner <- function(bundle, spec, cfg, term = cfg$primary_term, by = "cohort", variant = list()) {
  spec <- data.table::as.data.table(spec)
  dat  <- build_onset_data(bundle, spec, cfg, variant)
  ycol <- attr(dat, "ycol")
  if (!enough_events(event_counts(dat, ycol), cfg) || uniqueN(dat[[by]]) < 2L) {
    return(data.table(term = term, by = by, chisq = NA_real_, df = NA_integer_, p = NA_real_, status = "not_estimable"))
  }
  drop_cols <- constant_adjuster_cols(dat, spec$adjusted)
  f0 <- sub(" \\+ \\(1 \\| person_id\\)$", "", build_onset_formula(spec$predictors[[1L]], spec$adjusted, spec$k, drop_cols = drop_cols))
  if (by != "cohort") f0 <- paste0(f0, " + ", by)
  l0 <- fit_loglik(paste0(f0, " + (1 | person_id)"), dat, cfg)
  l1 <- fit_loglik(paste0(f0, " + ", term, ":", by, " + (1 | person_id)"), dat, cfg)     # one model in memory at a time
  chisq <- max(0, 2 * (l1$ll - l0$ll))
  df <- l1$df - l0$df
  data.table(term = term, by = by, chisq = chisq, df = df,
             p = stats::pchisq(chisq, df, lower.tail = FALSE), status = "ok")
}

test_group_interaction <- function(bundle, spec, cfg, term = cfg$primary_term, by = "cohort", variant = list()) {
  # two fits run one after the other in the same worker, whose memory is not returned in between (measured: ~15% more)
  check_fit_memory(bundle, as.list(spec), cfg, variant, extra_cols = 3L, safety = 1.25)
  isolated(cfg, test_group_interaction_inner(bundle, spec, cfg, term, by, variant))
}

test_cohort_interaction <- function(bundle, spec, cfg, term = cfg$primary_term) {
  r <- test_group_interaction(bundle, spec, cfg, term, by = "cohort")
  r[, by := NULL][]
}

# ---- collinearity check ---------------------------------------------------------------------------

# Variance inflation factor of each predictor (1 / (1 - R^2) from regressing it on the others).
vif_table <- function(dat, cols) {
  cols <- cols[cols %in% names(dat)]
  if (length(cols) < 2L) return(data.table(term = cols, vif = NA_real_))
  R <- stats::cor(dat[, ..cols], use = "pairwise.complete.obs")   # only a small matrix, no big copies
  v <- tryCatch(diag(solve(R)), error = function(e) rep(NA_real_, length(cols)))
  data.table(term = cols, vif = unname(v))
}

# ---- the onset stage for ONE stream, as independent checkpointed units ------------------------------

# Spec identity used in checkpoint keys (what makes one model different from another).
spec_key <- function(g_row) list(stream = g_row$stream, set_name = g_row$set_name, predictors = g_row$predictors[[1L]],
                                 k = g_row$k, adjusted = g_row$adjusted, role = g_row$role)

# Observed (model-free) gap rate by decile of the within-person z-score of the primary term on the predictor day (the last full day d - 1), within each cohort:
# deciles are cut over all rows at risk, so a decile holds the same days in every cohort. For each decile and cohort: mean z,
# rows, gaps that started, people, the rate and its 95% interval that respects several days per person (ratio estimator,
# person-clustered sandwich variance; see cluster_mean_ci). Unadjusted: it shows the shape of the association that the
# odds ratio summarises, not the adjusted effect.
onset_deciles <- function(bundle, cfg, k = cfg$onset$primary_k, term = cfg$primary_term, n_bins = 10L) {
  ycol <- paste0("y_ge", k)
  rs <- cbind(bundle$rs[, c("person_id", "cohort", ycol), with = FALSE], rs_at_predictor_day(bundle$rs, term, cfg))     # the term on the model's predictor day
  setnames(rs, c(term, ycol), c("z", "y"))
  rs <- rs[!is.na(z) & !is.na(y)]
  if (nrow(rs) == 0L) return(data.table())
  br <- unique(stats::quantile(rs$z, probs = seq(0, 1, length.out = n_bins + 1L), names = FALSE))
  if (length(br) < 3L) return(data.table())
  rs[, decile := findInterval(z, br, rightmost.closed = TRUE, all.inside = TRUE)]
  out <- rs[, {
    ci <- cluster_mean_ci(y, person_id)
    .(mean_z = mean(z), n_rows = .N, n_events = sum(y), n_persons = ci$n_persons,
      rate = ci$mean, rate_low = max(0, ci$lo), rate_high = min(1, ci$hi))
  }, by = .(cohort, decile)]
  out[order(cohort, decile)][]
}

# Difference of two cluster-sampled means (people are the clusters): mean(x1) - mean(x0) with the linearised cluster-robust variance
#   Var = g / (g - 1) * sum_i (u1_i - u0_i)^2,   u_ji = (S_ji - m_j n_ji) / N_j     (S = person sum, n = person rows, N = all rows of the group)
# The two groups come from the SAME people (a person has days followed by a gap and days that are not), so their covariance is included.
cluster_diff_ci <- function(x1, id1, x0, id0, level = 0.95) {
  if (length(x1) < 2L || length(x0) < 2L) return(list(difference = NA_real_, low = NA_real_, high = NA_real_, se = NA_real_, p = NA_real_))
  m1 <- mean(x1); m0 <- mean(x0); N1 <- length(x1); N0 <- length(x0)
  S1 <- rowsum(x1, id1); n1 <- rowsum(rep(1, N1), id1); S0 <- rowsum(x0, id0); n0 <- rowsum(rep(1, N0), id0)
  ids <- union(rownames(S1), rownames(S0)); g <- length(ids)
  u1 <- numeric(g); u0 <- numeric(g); names(u1) <- names(u0) <- ids
  u1[rownames(S1)] <- (S1[, 1L] - m1 * n1[, 1L]) / N1; u0[rownames(S0)] <- (S0[, 1L] - m0 * n0[, 1L]) / N0
  se <- sqrt(g / (g - 1) * sum((u1 - u0)^2)); d <- m1 - m0; q <- stats::qnorm(1 - (1 - level) / 2)
  list(difference = d, low = d - q * se, high = d + q * se, se = se, p = 2 * stats::pnorm(abs(d / se), lower.tail = FALSE))
}

# OBSERVED pre-gap trajectory of the primary term (the onset model's counterpart of the observed post-gap activity of model B): the within-person z-score
# (against the person's usual BEFORE that day, as in the model) on the day d - j, j = 0 .. 7, where d is the last day with a record before a gap
# (day 0; the gap starts on d + 1). Events = days d followed by a gap, by the length of that gap (1, 2, 3, 4-7, 8-14, 15+ days, all gaps; the length is known
# only afterwards, it describes the events and is never a predictor); controls = every other day d, on which no gap starts. Means with person-clustered
# intervals, and the difference events - controls at each day with its cluster-robust interval (the same people contribute both). A day d - j without a
# record has no z-score: the mean at that day uses the events that have one (n is reported). Descriptive and unadjusted; one value per (stream, group, day).
onset_pre_gap_trajectory <- function(bundle, cfg, k = cfg$onset$primary_k, term = cfg$primary_term, days_before = 0:7) {
  ycol <- paste0("y_ge", k); s <- bundle$stream
  base <- bundle$rs[, c("person_id", "date", "gap_len_next", ycol, term), with = FALSE]
  setkey(base, person_id, date)
  y <- base[[ycol]]; pid <- base$person_id
  grp <- as.character(gap_length_bin(pmax(base$gap_len_next, 1L), cfg))
  groups <- c(levels(gap_length_bin(1L, cfg)), "all gaps")
  ev_all <- which(y == 1L); ct <- which(y == 0L)
  ev_grp <- c(split(ev_all, grp[ev_all]), list(`all gaps` = ev_all))
  out <- list()
  for (j in days_before) {
    zj <- if (j == 0L) base[[term]] else base[.(pid, base$date - j), term, with = FALSE, on = c("person_id", "date")][[1L]]
    okc <- ct[!is.na(zj[ct])]
    cc <- if (length(okc) >= 2L) cluster_mean_ci(zj[okc], pid[okc]) else NULL
    for (g in groups) {
      ei <- ev_grp[[g]]; ei <- ei[!is.na(zj[ei])]
      if (length(ei) < 2L || is.null(cc)) next
      ce <- cluster_mean_ci(zj[ei], pid[ei]); df <- cluster_diff_ci(zj[ei], pid[ei], zj[okc], pid[okc])
      out[[length(out) + 1L]] <- data.table(group = g, days_before = j, mean = ce$mean, lo = ce$lo, hi = ce$hi, n_gaps = ce$n_gaps, n_persons = ce$n_persons,
                                            median_gap_length = as.numeric(stats::median(base$gap_len_next[ei])),
                                            mean_control = cc$mean, lo_control = cc$lo, hi_control = cc$hi, n_control_rows = cc$n_gaps, n_control_persons = cc$n_persons,
                                            difference = df$difference, diff_low = df$low, diff_high = df$high, p = df$p)
    }
  }
  res <- rbindlist(out)
  if (nrow(res)) res[, group := factor(group, levels = groups)]
  res[]
}

# Fit (or resume) every spec of the grid for one stream, plus the per-stream extras. Each model is
# one unit with its own checkpoint, so a crash or a failed model never costs the others.
#   s      : the stream name
#   bundle : lazy(prepare_stream_bundle(...)) - built only if some unit really has to run, then released
#   dkey   : data_key(store, stream), the input fingerprints the checkpoints are tied to
# Returns the run_unit() status of every unit (character vector); results are read with collect_onset_results().
run_onset_stage <- function(s, bundle, grid, cfg, ckpt_dir, dkey) {
  g <- grid[stream == s]
  status <- character()
  for (i in seq_len(nrow(g))) {
    cfg_i <- cfg; cfg_i$model$sample_persons <- effective_sample(cfg, g$role[i])
    u <- run_unit(ckpt_dir, list(kind = "onset", stream = s, id = g$spec_id[i]),
                  ckpt_key(cfg_i, "onset", dkey, spec_key(g[i]), cfg_i$model$sample_persons),
                  function() { r <- fit_onset_spec(bundle$get(), g[i], cfg_i); r$fit <- NULL; r }, cfg)
    status[g$spec_id[i]] <- u$status
  }

  # per-stream extras, each its own unit
  prim <- g[role == "primary"]
  allp <- unique(c(cfg$own_predictors[[s]], cfg$activity_predictors))
  extra <- function(id, compute, spec = NULL) {
    u <- run_unit(ckpt_dir, list(kind = "onset_extra", stream = s, id = id),
                  ckpt_key(cfg, "onset_extra", dkey, spec, cfg$model$sample_persons), compute, cfg)
    status[paste0("extra:", id)] <<- u$status
  }
  extra("vif", function() cbind(stream = s, vif_table(bundle$get()$rs, paste0("w_", allp))))
  if (nrow(prim) == 1L) {
    extra("cohort_interaction", function() cbind(stream = s, test_cohort_interaction(bundle$get(), prim, cfg)), spec_key(prim))
    for (co in cfg$cohorts) {
      extra(paste0("cohort_", co), function() {
        r <- fit_onset_spec(bundle$get(), prim, cfg, variant = list(cohort = co))
        if (is.null(r$coef)) cbind(r$info, term = cfg$primary_term) else r$coef[term == cfg$primary_term]
      }, spec_key(prim))
    }
    # observed gap rate by decile of the primary term: a model-free check of the odds ratio (one extra figure per stream)
    extra("deciles", function() cbind(stream = s, onset_deciles(bundle$get(), cfg, k = prim$k)), c(spec_key(prim), list(schema = 2L)))
    # observed steps on the days before a gap (the counterpart of the observed post-gap activity of model B): one more figure per stream
    extra("pre_gap_trajectory", function() cbind(stream = s, onset_pre_gap_trajectory(bundle$get(), cfg, k = prim$k)), c(spec_key(prim), list(schema = 1L)))
    # the absolute-risk curve comes from the unadjusted joint model, whatever roles this call fits
    unadj <- make_onset_grid(cfg)[stream == s & set_name == "joint" & k == cfg$onset$primary_k & !adjusted][, role := "risk_curve"]
    if (nrow(unadj) == 1L) {
      extra("risk_unadjusted", function() {
        r <- fit_onset_spec(bundle$get(), unadj, cfg)
        if (is.null(r$coef)) return(data.table(stream = s, status = r$info$status, note = r$info$note))
        cbind(stream = s, marginal_risk(r$coef, cfg$primary_term, r$info$sigma2, cov_b0_term = r$info$cov_b0_term))
      }, spec_key(unadj))
    }
  }
  status
}

# Gather everything finished so far (all streams, all roles) from the checkpoint folder. Only results
# made with the CURRENT settings and the CURRENT input data are used, so each call shows the complete picture.
# `failures` lists units whose last attempt failed (they are re-run on the next call).
collect_onset_results <- function(ckpt_dir, cfg, data_keys) {
  sp <- collect_units(ckpt_dir, "onset", cfg, data_keys)
  ex <- collect_units(ckpt_dir, "onset_extra", cfg, data_keys)
  pick <- function(ids) rbindlist(ex$values[ex$units$id %in% ids], fill = TRUE)
  out <- list(onset_coefficients = rbindlist(lapply(sp$values, `[[`, "coef"), fill = TRUE),
              onset_models       = rbindlist(lapply(sp$values, `[[`, "info"), fill = TRUE),
              onset_collinearity = rbindlist(lapply(sp$values, `[[`, "gvif"), fill = TRUE),
              onset_term_tests   = rbindlist(lapply(sp$values, `[[`, "wald"), fill = TRUE),
              vif = pick("vif"), cohort_interaction = pick("cohort_interaction"),
              cohort_specific = pick(paste0("cohort_", cfg$cohorts)), risk_unadjusted = pick("risk_unadjusted"),
              onset_deciles = pick("deciles"), onset_pre_gap_trajectory = pick("pre_gap_trajectory"),
              onset_failures = rbindlist(list(sp$failures, ex$failures), fill = TRUE))
  if (nrow(out$onset_pre_gap_trajectory) > 0L)                                  # Holm over the 8 days, within stream and gap-length group
    out$onset_pre_gap_trajectory[, p_holm := stats::p.adjust(p, "holm"), by = .(stream, group)]
  prim <- if (nrow(out$onset_coefficients) > 0L) out$onset_coefficients[role == "primary" & term == cfg$primary_term]
          else out$onset_coefficients
  if (nrow(prim) > 0L) {
    prim[, p_holm := stats::p.adjust(p, "holm")]
    out$primary_test <- prim[, .(stream, term, estimate, or, or_low, or_high, p, p_holm, n_events, status)]
  }
  out
}
