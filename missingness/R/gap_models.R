# ==============================================================================
# gap_models.R
# ------------------------------------------------------------------------------
# Models B and C: they use gap LENGTH, which is only known after a gap ends, so
# gap length is never a predictor of whether a gap starts (that would use the
# answer to predict the answer). One row per GAP, not per day.
#
#   B  "Is activity lower after longer gaps?"
#        outcome   : activity on the first FULL day after return (within-person z-score)
#        predictor : log2(gap length)  -> effect per DOUBLING of the gap
#        The first day back is a sensitivity (the watch is often put on mid-day).
#
#   C  "Once a gap starts, does low activity the day before predict a longer gap?"
#        outcome   : gap length group  1 | 2 | 3 | 4-7 | 8-14 | 15+
#        predictor : prior-day activity (within-person z-score)
#        ordinal mixed model - run ONLY if the proportional-odds assumption holds.
#
# Both are conditional on a gap having happened (so they describe gap-prone
# people) and are exploratory / secondary.
#
# Memory: these models need one row per gap (about 10^5 rows, a few MB). They are built from the
# stream's observed dates and ONE activity predictor, never from the person-day risk set, and every
# model is a separate checkpointed unit that is fitted in its own worker process.
# ==============================================================================

# Label the gap lengths into the ordered groups used by model C.
gap_length_bin <- function(len, cfg) {
  br <- c(cfg$gap_models$length_bins, Inf)       # 1,2,3,4,8,15,Inf
  cut(len, breaks = br, right = FALSE, ordered_result = TRUE,
      labels = c("1", "2", "3", "4-7", "8-14", "15+"))
}

# Merge length groups that contain fewer than `min_n` gaps into the next SHORTER group
# (so the ordinal model never sees an almost-empty category). Labels show what was merged,
# e.g. "4-7+8-14".
collapse_sparse_bins <- function(f, min_n) {
  lv <- levels(f); cnt <- as.integer(table(f)[lv]); keep <- rep(TRUE, length(lv))
  for (i in rev(seq_along(lv)[-1L])) {
    if (cnt[i] < min_n) { cnt[i - 1L] <- cnt[i - 1L] + cnt[i]; cnt[i] <- 0L; keep[i] <- FALSE }
  }
  target <- vapply(seq_along(lv), function(i) { j <- i; while (!keep[j]) j <- j - 1L; j }, integer(1))
  new_lv <- vapply(which(keep), function(j) paste(lv[target == j], collapse = "+"), character(1))
  factor(new_lv[match(target[as.integer(f)], which(keep))], levels = new_lv, ordered = TRUE)
}

# ---- the gap table ------------------------------------------------------------------------------------

# One row per gap with the activity z-scores around it, built from
#   days      the stream's day table (person_id, date)
#   act_days  the activity day table (person_id, date + the predictor's source column)
#   pred      the activity predictor (a name from cfg$predictors)
#   standardization  "pre_gap" (primary) or "full_period" (sensitivity), see predictors.R
# Columns (z = SD of the person's own day-to-day variation; + = more steps than usual):
#   w_pre_day_steps : z of the last observed day before the gap, against the person's usual BEFORE that day
#   w_pre_day_mvpa  : the same for MVPA minutes (log1p) on that day: MVPA of the last observed day before the gap against the
#                     person's usual MVPA BEFORE that day (identical timing and standardization as w_pre_day_steps); NA if the
#                     activity table has no MVPA column. It is the MVPA covariate of the `*_mvpa` units of models B and C.
#   post_k0 .. post_k7 (pre_gap only): z of the day r + k, r = first observed day after the gap (the "return day", k = 0),
#                     against the person's usual THROUGH the last observed day before the gap (data observed before the gap);
#                     NA when day r + k has no activity record
#   back1_dev       : first day back (= post_k0; sensitivity for B)
#   back2_dev       : first FULL day back = r + 1 (= post_k1; primary outcome of B)
# With "full_period" back1_dev / back2_dev are z-scores against the person's whole-period mean (look-ahead).
# Gaps without a usual (fewer than baseline_min_n observed days before them) have NA z-scores.
# Only the columns the models use are kept (no dates, no unused predictors).
POST_GAP_DAYS <- 0:7
build_gap_level_data <- function(days, act_days, persons, cfg, pred = "act_steps",
                                 standardization = cfg$standardization$method %||% "pre_gap") {
  gaps <- get_gap_runs(days)
  pre_gap <- standardization == "pre_gap"
  mvpa <- "act_mvpa"
  has_mvpa <- cfg$predictors[[mvpa]]$col %in% names(act_days)
  a <- predictor_table(act_days, c(pred, if (has_mvpa) mvpa), cfg, trailing = FALSE, keep_x = pre_gap,
                       standardization = standardization, keep_baseline = pre_gap)
  sd_w <- attr(a, "scales")[predictor == pred, sd_within]                 # the steps scale, as before (each predictor has its own constants)
  wcol <- paste0("w_", pred)
  a_mvpa <- if (has_mvpa) a[, c("person_id", "date", paste0("w_", mvpa)), with = FALSE]
  a <- a[, c("person_id", "date", wcol, if (pre_gap) c(paste0("x_", pred), paste0("m_", pred))), with = FALSE]
  setnames(a, c(wcol, if (pre_gap) c(paste0("x_", pred), paste0("m_", pred))), c("w", if (pre_gap) c("x", "m")))
  setkey(a, person_id, date)

  gaps[, w_pre_day_steps := a[.(gaps$person_id, gaps$prev_obs_date), w]]
  if (has_mvpa) {                                                          # MVPA of the SAME day (the last observed day before the gap), same rule
    setnames(a_mvpa, paste0("w_", mvpa), "w"); setkey(a_mvpa, person_id, date)
    gaps[, w_pre_day_mvpa := a_mvpa[.(gaps$person_id, gaps$prev_obs_date), w]]
    rm(a_mvpa)
  } else gaps[, w_pre_day_mvpa := NA_real_]
  if (pre_gap) {
    m_g <- a[.(gaps$person_id, gaps$prev_obs_date), m]               # usual through the last day before the gap
    for (k in POST_GAP_DAYS)
      set(gaps, j = paste0("post_k", k), value = (a[.(gaps$person_id, gaps$next_obs_date + k), x] - m_g) / sd_w)
    gaps[, `:=`(back1_dev = post_k0, back2_dev = post_k1)]
  } else {
    gaps[, `:=`(back1_dev = a[.(gaps$person_id, gaps$next_obs_date), w],
                back2_dev = a[.(gaps$person_id, gaps$next_obs_date + 1L), w])]
  }
  rm(a)
  gaps[, `:=`(log2_gap = log2(gap_length),
              length_bin = collapse_sparse_bins(gap_length_bin(gap_length, cfg), cfg$gap_models$min_bin_n))]
  gaps[, c("gap_index", "prev_obs_date", "gap_start", "gap_end", "next_obs_date") := NULL]

  cov <- person_covariates(persons, cfg)                 # cohort, age10, sex2, CCI and device family (covariates.R)
  merge(gaps, cov, by = "person_id", all.x = TRUE, sort = FALSE)[]
}

# Covariate text for the adjusted gap-level models; constant ones are dropped.
gap_adjusters <- function(gdat) {
  cand <- intersect(c("cohort", "age10", "sex2", "device", "cci", "cci_missing"), names(gdat))     # device / CCI only when cfg$model switches them on
  cand[vapply(cand, function(cl) uniqueN(gdat[[cl]][!is.na(gdat[[cl]])]) > 1L, logical(1))]
}

# Mean of y with a 95% interval that respects the clustering of gaps within people
# (ratio estimator, sandwich variance over persons).
cluster_mean_ci <- function(y, id) {
  n <- length(y); m <- mean(y)
  per <- rowsum(y - m, id)[, 1L]
  g <- length(per)
  se <- if (g > 1L) sqrt(g / (g - 1) * sum(per^2)) / n else NA_real_
  list(mean = m, lo = m - 1.96 * se, hi = m + 1.96 * se, n_gaps = n, n_persons = g)
}

# Observed mean outcome by gap-length group (what the model curve of B is compared with).
observed_by_gap_length <- function(d, outcome, cfg) {
  d <- copy(d)[, grp := gap_length_bin(gap_length, cfg)]
  d[, {
    ci <- cluster_mean_ci(get(outcome), person_id)
    .(x = as.numeric(stats::median(gap_length)), mean = ci$mean, lo = ci$lo, hi = ci$hi,
      n_gaps = ci$n_gaps, n_persons = ci$n_persons)
  }, by = grp][order(x)][]
}

# ---- Model B ---------------------------------------------------------------------------------------

empty_B_row <- function(base, status, note = "") {
  cbind(base, term = "log2_gap", estimate = NA_real_, se = NA_real_, conf_low = NA_real_, conf_high = NA_real_,
        p = NA_real_, singular = NA, converged = NA, sigma_person = NA_real_, sigma_resid = NA_real_,
        status = status, note = note)
}

# Fitted mean outcome at gap lengths 1..max, with a 95% interval from the fixed-effect covariance.
# The model is linear, so this is also the population-average mean: covariates are held at their
# mean in the analysed gaps (the share of each cohort / sex, the mean age).
gap_curve_B <- function(fit, cfg) {
  X <- lme4::getME(fit, "X"); xbar <- colMeans(X)
  beta <- lme4::fixef(fit); V <- as.matrix(stats::vcov(fit))
  L <- seq_len(cfg$predict$gap_curve_max)
  rbindlist(lapply(L, function(l) {
    x <- xbar; x["log2_gap"] <- log2(l)
    est <- sum(x * beta[names(x)]); se <- sqrt(drop(t(x) %*% V[names(x), names(x)] %*% x))
    data.table(gap_length = l, est = est, lo = est - 1.96 * se, hi = est + 1.96 * se)
  }))
}

# outcome: "back2_dev" (first full day, primary) or "back1_dev" (first day back).
# adjust_pre = TRUE also adjusts for the pre-gap deviation ("change from before to after").
# adjust_mvpa = TRUE also adjusts for MVPA on the last observed day BEFORE the gap (w_pre_day_mvpa: the same day and the same
#   pre-gap standardization as w_pre_day_steps, never the days after the gap, which are the outcome); gaps without it are
#   dropped (n_gaps_no_mvpa counts them).
# Returns list(coef, curve, observed, info); status "failed..." is never reported as a fit.
fit_gap_dose_response <- function(gdat, cfg, outcome = "back2_dev", adjusted = FALSE, adjust_pre = FALSE,
                                  id = NA_character_, adjust_mvpa = FALSE) {
  d <- gdat[!is.na(get(outcome)) & !is.na(log2_gap)]
  if (adjust_pre) d <- d[!is.na(w_pre_day_steps)]
  if (adjusted) d <- d[complete.cases(d[, .(age10, sex2, cohort)])]
  n_no_mvpa <- if (adjust_mvpa) sum(is.na(d$w_pre_day_mvpa)) else 0L
  if (adjust_mvpa) d <- d[!is.na(w_pre_day_mvpa)]
  base <- data.table(model = "B", model_id = id, outcome = outcome, adjusted = adjusted, adjust_pre = adjust_pre,
                     adjust_mvpa = adjust_mvpa, n_gaps_no_mvpa = n_no_mvpa,
                     n_gaps = nrow(d), n_persons = uniqueN(d$person_id),
                     outcome_mean = mean(d[[outcome]]), outcome_sd = stats::sd(d[[outcome]]),    # outcome distribution (z-scores: ~0 and ~1)
                     gaps_per_person_max = if (nrow(d)) max(table(d$person_id)) else NA_integer_)
  finish <- function(coef, coef_all = NULL) list(coef = coef, curve = NULL, observed = NULL, coef_all = coef_all,
                                                  gvif = NULL, wald = NULL, info = coef[, .(model_id, status, note)])
  if (nrow(d) < cfg$gap_models$min_gaps) return(finish(empty_B_row(base, "not_estimable", "too few gaps")))

  rhs <- c("log2_gap", if (adjust_pre) "w_pre_day_steps", if (adjust_mvpa) "w_pre_day_mvpa", if (adjusted) gap_adjusters(d), "(1 | person_id)")
  notes <- character()
  fit <- tryCatch(withCallingHandlers(
    suppressMessages(lme4::lmer(stats::as.formula(paste(outcome, "~", paste(rhs, collapse = " + "))), data = d)),   # singularity is reported in `singular`
    warning = function(w) { notes <<- c(notes, conditionMessage(w)); invokeRestart("muffleWarning") }),
    error = function(e) e)
  if (inherits(fit, "error")) return(finish(empty_B_row(base, paste("failed:", conditionMessage(fit)))))

  conv <- lme4_converged(fit)
  vc <- as.data.frame(lme4::VarCorr(fit))
  coef <- cbind(base, fixed_effect_row(fit, "log2_gap"), singular = lme4::isSingular(fit), converged = conv,
                sigma_person = vc$sdcor[vc$grp == "person_id"], sigma_resid = vc$sdcor[vc$grp == "Residual"],
                status = if (conv) "ok" else "not_converged",
                note = paste(c(if (lme4::isSingular(fit)) "singular fit (person variance ~ 0)", utils::head(unique(notes), 3L)),
                             collapse = " | "))
  curve <- cbind(model_id = id, outcome = outcome, adjusted = adjusted, adjust_pre = adjust_pre, gap_curve_B(fit, cfg))
  obs   <- cbind(model_id = id, outcome = outcome, observed_by_gap_length(d, outcome, cfg))
  # every fixed effect of the model (not the intercept), for the forest plot of all predictors
  all_fx <- fixed_effects_table(fit)
  all_fx <- cbind(model_id = id, outcome = outcome, adjusted = adjusted, adjust_pre = adjust_pre, status = coef$status,
                  all_fx[term != "(Intercept)"])
  dg <- fixed_effect_diagnostics(lme4::fixef(fit), as.matrix(stats::vcov(fit)), paste(outcome, "~", paste(setdiff(rhs, "(1 | person_id)"), collapse = " + ")))
  list(coef = coef, curve = curve, observed = obs, coef_all = all_fx,
       gvif = tag_diagnostics(dg$gvif, model = "B", model_id = id), wald = tag_diagnostics(dg$wald, model = "B", model_id = id),
       info = coef[, .(model_id, status, note)])
}

# All fixed effects of an lmer fit with Wald 95% intervals (normal approximation, as in fixed_effect_row).
fixed_effects_table <- function(fit) {
  cf <- coef(summary(fit))
  data.table(term = rownames(cf), estimate = cf[, 1], se = cf[, 2],
             conf_low = cf[, 1] - 1.96 * cf[, 2], conf_high = cf[, 1] + 1.96 * cf[, 2],
             p = 2 * stats::pnorm(abs(cf[, 1] / cf[, 2]), lower.tail = FALSE))
}

fixed_effect_row <- function(fit, term) {
  cf <- coef(summary(fit))
  data.table(term = term, estimate = cf[term, 1], se = cf[term, 2],
             conf_low = cf[term, 1] - 1.96 * cf[term, 2], conf_high = cf[term, 1] + 1.96 * cf[term, 2],
             p = 2 * stats::pnorm(abs(cf[term, 1] / cf[term, 2]), lower.tail = FALSE))
}

# ---- Model C -----------------------------------------------------------------------------------------

# Test of proportional odds: one odds ratio for the prior-day z-score across all
# cut points (no-random-effect version; a small p means the assumption fails).
check_proportional_odds <- function(d, cfg, adj) {
  rhs <- paste(adj, collapse = " + ")
  f_po  <- stats::as.formula(paste("length_bin ~ w_pre_day_steps", if (nzchar(rhs)) paste("+", rhs)))
  f_loc <- stats::as.formula(paste("length_bin ~", if (nzchar(rhs)) rhs else "1"))
  tryCatch({
    m_po  <- ordinal::clm(f_po, data = d)
    m_nom <- ordinal::clm(f_loc, nominal = ~ w_pre_day_steps, data = d)
    p <- stats::anova(m_po, m_nom)[["Pr(>Chisq)"]][2L]
    data.table(po_p = p, po_ok = is.finite(p) && p > cfg$gap_models$po_alpha, note = "")
  }, error = function(e) data.table(po_p = NA_real_, po_ok = FALSE, note = conditionMessage(e)))
}

# Ordinal model without random effect; the draws are the whole parameter vector (thresholds +
# coefficients) refitted on bootstrap samples of PEOPLE. Replicate b uses seed cfg$seed + b, so the
# result does not depend on how many workers run it. Workers are capped by free memory.
ordinal_cluster_bootstrap <- function(d, f_fixed, cfg) {
  f <- stats::as.formula(f_fixed)
  par_hat <- stats::coef(ordinal::clm(f, data = d))
  rows_by_person <- split(seq_len(nrow(d)), d$person_id)
  ids <- names(rows_by_person)
  one <- function(b) {
    pick <- with_seed(cfg$seed + b, sample(ids, length(ids), replace = TRUE))
    fit <- tryCatch(ordinal::clm(f, data = d[unlist(rows_by_person[pick], use.names = FALSE)]),
                    error = function(e) NULL)
    if (is.null(fit) || !identical(names(stats::coef(fit)), names(par_hat))) rep(NA_real_, length(par_hat))
    else unname(stats::coef(fit))
  }
  cores <- if (.Platform$OS.type == "unix") n_workers(cfg, per_worker_gb = 0.5, requested = cfg$resources$light_workers) else 1L
  draws <- do.call(rbind, parallel::mclapply(seq_len(cfg$boot$B), one, mc.cores = cores))
  colnames(draws) <- names(par_hat)
  ok <- stats::complete.cases(draws)
  if (sum(ok) < 0.8 * cfg$boot$B) stop(sprintf("only %d of %d bootstrap fits succeeded", sum(ok), cfg$boot$B), call. = FALSE)
  list(par_hat = par_hat, draws = draws[ok, , drop = FALSE], se = apply(draws[ok, , drop = FALSE], 2L, stats::sd))
}

# Probability of each length group at given prior-day z-scores, averaged over (a) the covariate mix of
# the analysed gaps (a seeded sample of rows) and (b) the person random intercept (sigma2 = 0 for the
# model without one). These are POPULATION-AVERAGE probabilities, not those of one typical person.
#   par   : named vector (thresholds + coefficients)    Xadj : adjuster design matrix (rows = sampled gaps)
ordinal_probs <- function(par, alpha_names, Xadj, sigma2, w_values, n_u = 25L) {
  alpha <- par[alpha_names]
  eta_adj <- if (ncol(Xadj) > 0L) drop(Xadj %*% par[colnames(Xadj)]) else rep(0, nrow(Xadj))
  u <- if (sigma2 > 0) stats::qnorm(stats::ppoints(n_u)) * sqrt(sigma2) else 0
  t(vapply(w_values, function(w) {
    eta <- eta_adj + par[["w_pre_day_steps"]] * w
    acc <- 0
    for (uu in u) {
      cum <- cbind(0, stats::plogis(outer(eta + uu, alpha, function(e, a) a - e)), 1)
      acc <- acc + colMeans(cum[, -1L, drop = FALSE] - cum[, -ncol(cum), drop = FALSE])
    }
    acc / length(u)
  }, numeric(length(alpha) + 1L)))
}

# Point estimate + 95% interval (percentiles over parameter draws) of those probabilities.
ordinal_pred_table <- function(par_hat, par_draws, levels_, Xadj, sigma2, w_values, cfg) {
  alpha_names <- names(par_hat)[grepl("\\|", names(par_hat))]
  pt <- ordinal_probs(par_hat, alpha_names, Xadj, sigma2, w_values)
  dr <- lapply(seq_len(nrow(par_draws)), function(i) ordinal_probs(par_draws[i, ], alpha_names, Xadj, sigma2, w_values))
  arr <- simplify2array(dr)                                  # w x K x draws
  q <- apply(arr, c(1L, 2L), stats::quantile, probs = c(0.025, 0.975), names = FALSE)
  out <- CJ(within_z = w_values, length_group = levels_, sorted = FALSE)
  out[, `:=`(prob = as.vector(t(pt)), prob_low = as.vector(t(q[1L, , ])), prob_high = as.vector(t(q[2L, , ])))]
  out[, length_group := factor(length_group, levels = levels_, ordered = TRUE)]
  out[]
}

empty_C_row <- function(base, status, note = "", po_p = NA_real_) {
  cbind(base, term = "w_pre_day_steps", estimate = NA_real_, se = NA_real_, conf_low = NA_real_,
        conf_high = NA_real_, p = NA_real_, or = NA_real_, or_low = NA_real_, or_high = NA_real_,
        po_p = as.numeric(po_p), converged = NA, sigma_person = NA_real_, status = status, note = note)
}

# Ordinal mixed model for the gap-length group. Status "po_rejected" (and no estimate) when
# proportional odds fails - then model D is the analysis to report. Reported odds ratio: odds of a
# LONGER gap per +1 SD HIGHER than the person's usual activity. Returns list(coef, pred, info).
# adjust_mvpa = TRUE adds MVPA on the last observed day before the gap (w_pre_day_mvpa, same day and same pre-gap standardization
#   as the exposure w_pre_day_steps) as a covariate; it enters the proportional-odds check and the predicted probabilities as an
#   ordinary covariate (only the effect of the steps term is tested for proportionality). Gaps without it are dropped.
fit_gap_length_ordinal <- function(gdat, cfg, adjusted = FALSE, id = NA_character_, adjust_mvpa = FALSE) {
  d <- gdat[!is.na(w_pre_day_steps) & !is.na(length_bin)]
  adj <- if (adjusted) gap_adjusters(d) else character()
  if (adjusted) d <- d[complete.cases(d[, ..adj])]
  n_no_mvpa <- if (adjust_mvpa) sum(is.na(d$w_pre_day_mvpa)) else 0L
  if (adjust_mvpa) d <- d[!is.na(w_pre_day_mvpa)]
  cov_terms <- c(if (adjust_mvpa) "w_pre_day_mvpa", adj)                           # every covariate besides the exposure
  d <- d[, c("person_id", "length_bin", "w_pre_day_steps", cov_terms), with = FALSE]      # nothing else is needed
  d[, length_bin := droplevels(length_bin)]
  counts <- table(d$length_bin)
  base <- data.table(model = "C", model_id = id, adjusted = adjusted, adjust_mvpa = adjust_mvpa, n_gaps_no_mvpa = n_no_mvpa,
                     n_gaps = nrow(d), n_persons = uniqueN(d$person_id),
                     length_groups = paste(names(counts), collapse = " | "),
                     group_counts = paste(counts, collapse = " | "))
  finish <- function(coef, pred = NULL, coef_all = NULL, dg = NULL)
    list(coef = coef, pred = pred, coef_all = coef_all,
         gvif = if (!is.null(dg)) tag_diagnostics(dg$gvif, model = "C", model_id = id), wald = if (!is.null(dg)) tag_diagnostics(dg$wald, model = "C", model_id = id),
         info = coef[, .(model_id, status, note)])
  if (nrow(d) < cfg$gap_models$min_gaps) return(finish(empty_C_row(base, "not_estimable", "too few gaps")))
  if (nlevels(d$length_bin) < 3L) return(finish(empty_C_row(base, "not_estimable", "fewer than 3 length groups")))

  po <- check_proportional_odds(d, cfg, cov_terms)
  if (!isTRUE(po$po_ok)) return(finish(empty_C_row(base, "po_rejected", paste("use model D;", po$note), po$po_p)))

  f_fixed <- paste("length_bin ~ w_pre_day_steps", if (length(cov_terms)) paste("+", paste(cov_terms, collapse = " + ")))
  notes <- character()
  fit <- tryCatch(withCallingHandlers(
    ordinal::clmm(stats::as.formula(paste(f_fixed, "+ (1 | person_id)")), data = d,
                  control = ordinal::clmm.control(useMatrix = isTRUE(cfg$resources$ordinal_sparse))),
    warning = function(w) { notes <<- c(notes, conditionMessage(w)); invokeRestart("muffleWarning") }),
    error = function(e) e)

  par_names <- NULL; est <- se <- NA_real_; status <- "ok"; note <- ""; sigma2 <- 0; conv <- NA
  mixed_ok <- FALSE
  if (!inherits(fit, "error")) {
    par_hat <- c(fit$alpha, fit$beta); par_names <- names(par_hat)
    V <- tryCatch(as.matrix(stats::vcov(fit))[par_names, par_names, drop = FALSE], error = function(e) NULL)
    conv <- isTRUE(fit$optRes$convergence == 0)
    sigma2 <- as.numeric(ordinal::VarCorr(fit)$person_id[1L, 1L])
    mixed_ok <- !is.null(V) && all(is.finite(diag(V))) && conv
    if (mixed_ok) { est <- par_hat[["w_pre_day_steps"]]; se <- sqrt(V["w_pre_day_steps", "w_pre_day_steps"]) }
    if (!conv) note <- "mixed model did not converge"
  }
  Xadj <- if (length(cov_terms)) {
    ix <- with_seed(cfg$seed, sample.int(nrow(d), min(cfg$predict$n_covariate_rows, nrow(d))))
    mm <- stats::model.matrix(stats::as.formula(paste("~", paste(cov_terms, collapse = " + "))), d[ix])[, -1L, drop = FALSE]
    mm
  } else matrix(0, 1L, 0L)

  if (mixed_ok) {
    L <- tryCatch(chol(V), error = function(e) NULL)
    if (is.null(L)) mixed_ok <- FALSE
  }
  if (mixed_ok) {
    z <- with_seed(cfg$seed, matrix(stats::rnorm(cfg$predict$n_draws * length(par_hat)), cfg$predict$n_draws))
    par_draws <- sweep(z %*% L, 2L, par_hat, "+"); colnames(par_draws) <- par_names
    se_vec <- sqrt(diag(V)); V_use <- V
  } else {
    # The random-intercept variance is on the boundary (no person-to-person difference in gap length)
    # or the Hessian is unusable, so the mixed model has no covariance. Fall back to the same ordinal model
    # without the random effect and a person-level (cluster) bootstrap for all uncertainty.
    bt <- tryCatch(ordinal_cluster_bootstrap(d, f_fixed, cfg), error = function(e) e)
    if (inherits(bt, "error")) return(finish(empty_C_row(base, paste("failed:", conditionMessage(bt)), note, po$po_p)))
    par_hat <- bt$par_hat; par_draws <- bt$draws; sigma2 <- 0; se_vec <- bt$se
    V_use <- stats::cov(bt$draws)                       # covariance of the cluster-bootstrap draws
    est <- par_hat[["w_pre_day_steps"]]; se <- bt$se[["w_pre_day_steps"]]
    status <- "ok_bootstrap_se"
    note <- paste(c(note, "random-intercept variance on boundary or Hessian unusable; no random effect, cluster-bootstrap uncertainty"), collapse = "; ")
  }
  w_values <- c(-1, 0, 1)
  pred <- ordinal_pred_table(par_hat, par_draws, levels(d$length_bin), Xadj, sigma2, w_values, cfg)
  pred <- cbind(model_id = id, adjusted = adjusted, pred)

  z <- est / se
  # clm/clmm: a positive coefficient -> longer gaps, so the odds ratio of a longer gap per +1 SD is exp(estimate)
  coef <- cbind(base, data.table(term = "w_pre_day_steps", estimate = est, se = se,
                         conf_low = est - 1.96 * se, conf_high = est + 1.96 * se,
                         p = 2 * stats::pnorm(abs(z), lower.tail = FALSE), or = exp(est),
                         or_low = exp(est - 1.96 * se), or_high = exp(est + 1.96 * se),
                         po_p = po$po_p, converged = conv, sigma_person = sqrt(sigma2), status = status,
                         note = paste(c(note, utils::head(unique(notes), 3L)), collapse = " | ")))
  # every coefficient of the model (thresholds excluded) with its odds ratio of a LONGER gap: per +1 SD above usual for the
  # prior-day term (sign flipped, as in `coef`), per +1 unit / vs the reference level for the others
  b_names <- setdiff(names(par_hat), names(par_hat)[grepl("\\|", names(par_hat))])
  all_fx <- data.table(term = b_names, estimate = unname(par_hat[b_names]), se = unname(se_vec[b_names]))
  all_fx[, `:=`(conf_low = estimate - 1.96 * se, conf_high = estimate + 1.96 * se,
                p = 2 * stats::pnorm(abs(estimate / se), lower.tail = FALSE))]
  all_fx <- add_odds_ratios(all_fx)          # the `w_pre_day_steps` term is a within-person z-score: OR per +1 SD above usual
  all_fx <- cbind(model_id = id, adjusted = adjusted, status = status, all_fx)
  dg <- tryCatch(fixed_effect_diagnostics(par_hat[b_names], V_use[b_names, b_names, drop = FALSE], f_fixed), error = function(e) NULL)
  finish(coef, pred, all_fx, dg)
}

# ---- the gap-model stage for ONE stream, as independent checkpointed units ---------------------------------

# Model B has four specifications and model C two; each is its own unit, so a failure of one never
# forces a refit of another.
# Two sensitivity units use the full-period standardization (the person's usual from the WHOLE observation period) and one
# unit summarises the OBSERVED post-gap activity (no model); the others use the primary pre-gap standardization.
gap_unit_specs <- function() list(
  B_back2_unadj   = list(model = "B", outcome = "back2_dev", adjusted = FALSE, adjust_pre = FALSE),
  B_back2_adj     = list(model = "B", outcome = "back2_dev", adjusted = TRUE,  adjust_pre = FALSE),
  B_back2_adj_pre = list(model = "B", outcome = "back2_dev", adjusted = TRUE,  adjust_pre = TRUE),
  B_back1_adj     = list(model = "B", outcome = "back1_dev", adjusted = TRUE,  adjust_pre = FALSE),
  B_back2_adj_fullperiod = list(model = "B", outcome = "back2_dev", adjusted = TRUE, adjust_pre = FALSE, standardization = "full_period"),
  B_observed_trajectory  = list(model = "BOBS"),
  # MVPA before the gap as an additional covariate (same day and standardization as the pre-gap steps); the specifications above are unchanged
  B_back2_adj_mvpa     = list(model = "B", outcome = "back2_dev", adjusted = TRUE, adjust_pre = FALSE, adjust_mvpa = TRUE),
  B_back2_adj_pre_mvpa = list(model = "B", outcome = "back2_dev", adjusted = TRUE, adjust_pre = TRUE,  adjust_mvpa = TRUE),
  C_unadj         = list(model = "C", adjusted = FALSE),
  C_adj           = list(model = "C", adjusted = TRUE),
  C_adj_fullperiod = list(model = "C", adjusted = TRUE, standardization = "full_period"),
  C_adj_mvpa       = list(model = "C", adjusted = TRUE, adjust_mvpa = TRUE))

# Version of the stored B / C results (2: every coefficient of the model is kept for the forest plots). Part of the
# checkpoint key, so results stored by an older version are recomputed, never mixed in.
GAP_SCHEMA <- 2L

# gap_data : lazy(build_gap_level_data(...)) - built only when some unit has to run.
# gap_data_full : lazy(build_gap_level_data(..., standardization = "full_period")) for the full-period sensitivity units
#   (NULL: those units are not run). Built only when one of them has to run.
# The gap table is small (one row per gap) and a unit needs well under 1 GB, so the units that still
# have to run are executed side by side (cfg$resources$light_workers, capped by free memory). Model C
# (a mixed ordinal model) is by far the slowest, so it starts first. Each worker runs run_unit() for its
# own unit, i.e. reads / writes only that unit's checkpoint.
run_gap_stage <- function(s, gap_data, cfg, ckpt_dir, dkey, gap_data_full = NULL) {
  specs <- gap_unit_specs()
  if (is.null(gap_data_full)) specs <- specs[vapply(specs, function(sp) !identical(sp$standardization, "full_period"), logical(1))]
  ids <- names(specs)
  kind_of <- function(id) if (specs[[id]]$model == "C") "gapC" else "gapB"
  unit_of <- function(id) list(kind = kind_of(id), stream = s, id = id)
  key_of  <- function(id) ckpt_key(cfg, kind_of(id), dkey, c(specs[[id]], list(schema = GAP_SCHEMA)))
  data_of <- function(sp) if (identical(sp$standardization, "full_period")) gap_data_full$get() else gap_data$get()
  one <- function(id, cfg_u) {
    sp <- specs[[id]]
    run_unit(ckpt_dir, unit_of(id), key_of(id), function() {
      g <- data_of(sp)                           # built in the parent before the workers start (shared by fork)
      if (sp$model == "B") isolated(cfg_u, fit_gap_dose_response(g, cfg_u, sp$outcome, sp$adjusted, sp$adjust_pre, id, adjust_mvpa = isTRUE(sp$adjust_mvpa)))
      else if (sp$model == "C") isolated(cfg_u, fit_gap_length_ordinal(g, cfg_u, sp$adjusted, id, adjust_mvpa = isTRUE(sp$adjust_mvpa)))
      else isolated(cfg_u, list(trajectory = cbind(model_id = id, observed_post_gap_trajectory(g, cfg_u)),
                                info = data.table(model_id = id, status = "ok", note = "")))
    }, cfg_u)$status
  }
  finished <- vapply(ids, function(id) ckpt_read(ckpt_path(ckpt_dir, unit_of(id), key_of(id)), key_of(id))$status == "valid", logical(1))
  status <- setNames(rep("resumed", length(ids)), ids)
  todo <- ids[!finished]; todo <- todo[order(!startsWith(todo, "C_"))]
  for (id in ids[finished]) one(id, cfg)                                   # loads the checkpoint (logs the skip)
  if (length(todo)) {
    for (id in todo) data_of(specs[[id]])                                  # build the needed gap tables once, here in the parent
    cfg_u <- cfg; cfg_u$resources$light_workers <- 1L                      # no workers inside a worker
    cores <- if (.Platform$OS.type == "unix") n_workers(cfg, per_worker_gb = 1, requested = cfg$resources$light_workers) else 1L
    out <- parallel::mclapply(todo, one, cfg_u = cfg_u, mc.cores = min(cores, length(todo)), mc.preschedule = FALSE)
    status[todo] <- vapply(out, function(o) if (is.character(o) && length(o) == 1L) o else "failed", "")
  }
  status
}

# Observed (model-free) activity after gaps, by the length of the gap and the day since the return.
#   time origin  : k = 0 is the RETURN DAY, the first observed day after the gap; k = 1 is the next calendar day, ... k = 7
#   window       : the 8 days r .. r + 7 (a day without an activity record contributes nothing to that k)
#   value        : steps (log1p) as a z-score against the person's usual THROUGH the last observed day before the gap
#                  (pre_gap standardization), in units of the within-person SD; 0 = the person's usual, + = more steps
#   groups       : gap length 1, 2, 3, 4-7, 8-14, 15+ days and all gaps
# For each group and k: the mean, its 95% interval and the +/- 1 SD of the mean from the person-clustered ratio estimator
# (a person with several gaps counts as one cluster), the number of gaps and of people that contribute.
observed_post_gap_trajectory <- function(gdat, cfg, days = POST_GAP_DAYS) {
  cols <- paste0("post_k", days)
  if (!all(cols %in% names(gdat))) stop("the gap table has no post-gap columns: build it with the pre_gap standardization", call. = FALSE)
  g <- gdat[, c("person_id", "gap_length", cols), with = FALSE]
  g[, grp := as.character(gap_length_bin(gap_length, cfg))]
  groups <- c(levels(gap_length_bin(1L, cfg)), "all gaps")
  out <- rbindlist(lapply(groups, function(gr) {
    d <- if (gr == "all gaps") g else g[grp == gr]
    rbindlist(lapply(days, function(k) {
      y <- d[[paste0("post_k", k)]]; ok <- !is.na(y)
      if (sum(ok) < 2L) return(NULL)
      ci <- cluster_mean_ci(y[ok], d$person_id[ok])
      data.table(gap_group = gr, days_since_return = k, mean = ci$mean, se = (ci$hi - ci$lo) / (2 * 1.96), lo = ci$lo, hi = ci$hi,
                 n_gaps = ci$n_gaps, n_persons = ci$n_persons, median_gap_length = as.numeric(stats::median(d$gap_length[ok])))
    }))
  }))
  out[, gap_group := factor(gap_group, levels = groups)]
  out[]
}

# ---- MVPA before the gap: effect estimates of the `*_mvpa` units --------------------------------------------------------

# Model pairs: the MVPA-adjusted unit and the unit it extends (same specification without MVPA).
GAP_MVPA_PAIRS <- c(B_back2_adj_mvpa = "B_back2_adj", B_back2_adj_pre_mvpa = "B_back2_adj_pre", C_adj_mvpa = "C_adj")

# One table with the MVPA coefficient and the main terms with and without MVPA, with 95% intervals:
#   B (linear, SD of the person's usual steps): change in activity on the first full day back per +1 SD of pre-gap MVPA, per doubling
#   of the gap, per +1 SD of pre-gap steps;   C (cumulative logit): odds ratio of a LONGER gap per +1 SD of pre-gap MVPA / steps.
# vif_equiv = the VIF-equivalent of that term in its model (collinearity.R), n_gaps_no_mvpa = gaps dropped because MVPA is missing.
build_gap_mvpa_table <- function(b_all, c_all, b_coef, c_coef, gvif = NULL) {
  if (!nz(b_all) && !nz(c_all)) return(data.table())
  take <- function(tab, coef, model) {
    if (!nz(tab) || !nz(coef)) return(NULL)
    ids <- c(names(GAP_MVPA_PAIRS), unname(GAP_MVPA_PAIRS))
    x <- tab[model_id %in% ids & term %in% c("log2_gap", "w_pre_day_steps", "w_pre_day_mvpa")]
    if (nrow(x) == 0L) return(NULL)
    n <- unique(coef[, .(stream, model_id, n_gaps, n_persons, n_gaps_no_mvpa = if ("n_gaps_no_mvpa" %in% names(coef)) n_gaps_no_mvpa else 0L)])
    x <- merge(x[, .(stream, model_id, term, estimate, se, conf_low, conf_high, p, status, or = if ("or" %in% names(x)) or else NA_real_,
                     or_low = if ("or_low" %in% names(x)) or_low else NA_real_, or_high = if ("or_high" %in% names(x)) or_high else NA_real_)],
               n, by = c("stream", "model_id"), all.x = TRUE)
    x[, `:=`(model = model, mvpa_adjusted = model_id %in% names(GAP_MVPA_PAIRS))]
    x[, pair := ifelse(mvpa_adjusted, model_id, names(GAP_MVPA_PAIRS)[match(model_id, GAP_MVPA_PAIRS)])]
    x
  }
  out <- rbindlist(list(take(b_all, b_coef, "B"), take(c_all, c_coef, "C")), fill = TRUE)
  if (nrow(out) == 0L) return(data.table())
  out <- out[!is.na(pair)]
  out[, effect_measure := fifelse(model == "B", "SD change in activity (per +1 SD of the term; log2_gap: per doubling)", "odds ratio of a longer gap per +1 SD")]
  if (nz(gvif)) out <- merge(out, gvif[, .(stream, model_id, term, vif_equiv)], by = c("stream", "model_id", "term"), all.x = TRUE)
  setcolorder(out, c("model", "stream", "pair", "model_id", "mvpa_adjusted", "term", "effect_measure"))
  out[order(model, stream, pair, -mvpa_adjusted, term)][]
}

# Everything finished so far, all streams: coefficient tables, curves, predicted probabilities, failures.
collect_gap_results <- function(ckpt_dir, cfg, data_keys) {
  b <- collect_units(ckpt_dir, "gapB", cfg, data_keys); cc <- collect_units(ckpt_dir, "gapC", cfg, data_keys)
  tag <- function(vals, units, part) rbindlist(Map(function(v, st) if (is.null(v[[part]])) NULL else cbind(stream = st, v[[part]]),
                                                    vals, units$stream), fill = TRUE)
  out <- list(model_B_post_gap_activity = tag(b$values, b$units, "coef"),
       model_B_curve             = tag(b$values, b$units, "curve"),
       model_B_observed          = tag(b$values, b$units, "observed"),
       model_C_gap_length        = tag(cc$values, cc$units, "coef"),
       model_B_observed_trajectory = tag(b$values, b$units, "trajectory"),
       gap_collinearity          = rbindlist(list(tag(b$values, b$units, "gvif"), tag(cc$values, cc$units, "gvif")), fill = TRUE),
       gap_term_tests            = rbindlist(list(tag(b$values, b$units, "wald"), tag(cc$values, cc$units, "wald")), fill = TRUE),
       model_C_predicted_probs   = tag(cc$values, cc$units, "pred"),
       model_B_all_terms         = tag(b$values, b$units, "coef_all"),
       model_C_all_terms         = tag(cc$values, cc$units, "coef_all"),
       gap_failures              = rbindlist(list(b$failures, cc$failures), fill = TRUE))
  out$gap_mvpa_effects <- build_gap_mvpa_table(out$model_B_all_terms, out$model_C_all_terms, out$model_B_post_gap_activity,
                                               out$model_C_gap_length, out$gap_collinearity)
  out
}
