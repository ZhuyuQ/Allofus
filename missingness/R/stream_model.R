# ==============================================================================
# stream_model.R
# ------------------------------------------------------------------------------
# Does an unusually high-step day d predict the recording status of day d + 1? A risk set, a usual level built only from EARLIER days, within- and between-person
# components and the covariates of the onset model.
#
#   day status      activity present = a record remains after the activity filter; heart rate present = a record; sleep present = a MAIN sleep record dated
#                   that day (the preceding night). Days outside a stream's own first-to-last period, or on a device without a heart-rate sensor, are excluded
#                   (the person-day table of stream_states.R: reason == "analysed").
#   risk set        days d with all three streams present (fully recorded), at least wear_hours (20) of heart rate on d, and at least min_earlier_days (7)
#                   earlier days that met both conditions (the baseline days)
#   outcome         status of day d + 1: 0 fully recorded | 1 all missing (000) | 2 activity only (100) | 3 sleep missing (110) | 4 inconsistent (every other pattern:
#                   101, 001, 011, 010; negative control). Rows followed by an excluded day, and the last day of a window, are dropped and counted.
#   exposure        x = log(1 + steps) on d; m = mean of x over the person's EARLIER qualifying days (before d: the onset model's "pre_gap" usual level);
#                   w = (x - m) / sd_w (within person), b = (m - mu_b) / sd_b (between persons); sd_w, mu_b, sd_b are constants over the risk-set rows
# Model A: one mixed logistic model per category k (outcome k against 0, random person intercept): subject-specific odds ratios per +1 SD of w, H1 (category 1, two-sided Wald),
# equivalence (two one-sided tests, 90% interval within 0.97-1.03) for categories 2-4 (Holm across the three, unadjusted p kept).
# Model B: one multinomial logit (reference 0) with person-clustered standard errors: H2 (category 1 against 2 and 1 against 3, Holm over the two) and the predicted probabilities.
# Event-count thresholds decide which claims are allowed.
# ==============================================================================

suppressPackageStartupMessages(library(data.table))

STM_LABELS <- c(`0` = "Fully recorded", `1` = "All missing", `2` = "Activity only", `3` = "Sleep missing", `4` = "Inconsistent (negative control)")
STM_WEEKDAYS <- c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")
STM_COVARIATES <- c("cohort", "age10", "sex2", "cci", "cci_missing", "device", "dow_next", "prev_spell_l", "any_prior_spell", "dsr_l", "ns1", "ns2", "ns3")

stm_config <- function(overrides = list()) {
  utils::modifyList(list(wear_hours = 20, min_earlier_days = 7L, equiv_margin = c(0.97, 1.03), min_events_effect = 1500L, min_events_equiv = 6000L,
                         w_grid = seq(-2, 2, by = 0.5), n_profile_rows = 2000L, n_beta_draws = 2000L, n_re_grid = 41L, first_days_excluded = 14L,
                         verify_people = 1500L, deciles = 10L, probably_worn_cuts = c(20, 22), multinom_max_iter = 30L, seed = 20260511L), overrides)
}

# ---- days -----------------------------------------------------------------------------------------------------------------------------------------

# One row per person-day with the day status, the qualifying-day flag, the usual level of the exposure (EARLIER qualifying days only) and the recording history; the status
# and wear of the next row (d + 1) are attached. D: the person-day table of stream_states.R (reason assigned); naps: count naps as sleep.
stm_days <- function(D, scfg, naps = FALSE) {
  X <- D[, .(person_id, date, analysed = reason == "analysed", a = act_after, h = hr_present, s = if (naps) slp_row else sleep_main, steps, dev_version,
             wear = if ("hr_wear_min" %in% names(D)) hr_wear_min else hr_min)]
  setorder(X, person_id, date)
  X[, `:=`(full = analysed & a & h & s, study_day = as.integer(date - min(date)) + 1L), by = person_id]
  X[, status := fifelse(analysed, paste0(as.integer(a), as.integer(h), as.integer(s)), NA_character_)]
  X[, `:=`(qual = full & !is.na(steps) & (wear >= 60 * scfg$wear_hours) %in% TRUE, x = log1p(steps))]
  X <- stm_history(X)
  X[, `:=`(cs = cumsum(fifelse(qual, x, 0)), cn = cumsum(qual)), by = person_id]
  X[, `:=`(n_before = shift(cn, fill = 0L), m_before = shift(cs, fill = 0) / shift(cn, fill = 0L)), by = person_id]
  X[n_before == 0L, m_before := NA_real_]
  X[, `:=`(n_status = shift(status, -1L), n_wear = shift(wear, -1L), n_analysed = shift(analysed, -1L)), by = person_id]
  X[, c("cs", "cn") := NULL][]
}

# Recording history on these day statuses: the length of the most recent completed spell that was not fully recorded, whether there was an earlier spell, and the days since
# the last return to full recording (0 before the first spell). A spell is a run of analysed days that are not fully recorded; it is completed by the first fully
# recorded day after it (an excluded day ends a run without a return).
stm_history <- function(X) {
  X[, nonfull := analysed & !full]
  X[, `:=`(prev_nf = shift(nonfull, fill = FALSE), run = rleid(nonfull)), by = person_id]
  X[, run_len := fifelse(nonfull, .N, NA_integer_), by = .(person_id, run)]
  X[, ret := full & prev_nf]
  X[, spell_at_ret := fifelse(ret, shift(run_len), NA_integer_), by = person_id]
  X[, `:=`(prev_spell = nafill(spell_at_ret, "locf"), ret_day = nafill(fifelse(ret, as.integer(date), NA_integer_), "locf")), by = person_id]
  X[, `:=`(any_prior = !is.na(prev_spell), dsr = fifelse(!is.na(ret_day), as.integer(date) - ret_day, 0L))]
  X[, c("nonfull", "prev_nf", "run", "run_len", "ret", "spell_at_ret", "ret_day") := NULL][]
}

# ---- risk set ---------------------------------------------------------------------------------------------------------------------------------------

# Outcome category of the status of day d + 1 (NA = excluded day): 111 -> 0, 000 -> 1, 100 -> 2, 110 -> 3, every other pattern -> 4.
stm_outcome <- function(status) {
  fcase(is.na(status), NA_integer_, status == "111", 0L, status == "000", 1L, status == "100", 2L, status == "110", 3L, default = 4L)
}

# The risk set with outcome, exposure and covariates; flow = what each condition removed. persons: the person table; pcfg: the pipeline config.
stm_risk_set <- function(X, persons, pcfg, scfg) {
  n_full <- X[full %in% TRUE, .N]; cand <- X[qual %in% TRUE]; n_w <- nrow(cand)
  cand <- cand[n_before >= scfg$min_earlier_days]; n_b <- nrow(cand)
  cand <- cand[!is.na(n_analysed)]; n_l <- nrow(cand)                     # d is not the last day of the window
  cand <- cand[n_analysed == TRUE]; n_x <- nrow(cand)                     # d + 1 is not an excluded day
  cand[, outcome := stm_outcome(n_status)]
  flow <- data.table(step = c("days d with all three streams (fully recorded)", sprintf(">= %g h of heart rate on d (and steps)", scfg$wear_hours),
                              sprintf(">= %d earlier baseline days", scfg$min_earlier_days), "d is not the last day of the window", "d + 1 is not an excluded day"),
                     n_rows = c(n_full, n_w, n_b, n_l, n_x), n_dropped = c(NA, n_full - n_w, n_w - n_b, n_b - n_l, n_l - n_x))
  R <- stm_covariates(cand, persons, pcfg, scfg)
  list(R = R, flow = flow, scales = attr(R, "scales"))
}

# Exposure components, covariates and the spline basis of the risk-set rows; the scaling constants are returned as an attribute.
stm_covariates <- function(R, persons, pcfg, scfg) {
  cov <- person_covariates(persons, pcfg)[, .(person_id, cohort, sex2, age10, cci, cci_missing)]
  R <- cov[R, on = "person_id"]
  fam <- device_family(R$dev_version, pcfg); fam[is.na(fam)] <- "unknown"
  tab <- table(fam[fam != "unknown"]); ref <- if (length(tab)) names(which.max(tab)) else "unknown"
  R[, device := factor(fam, levels = unique(c(ref, names(sort(tab, decreasing = TRUE)), "unknown")))]
  R[, dow_next := factor(STM_WEEKDAYS[(data.table::wday(date + 1L) + 5L) %% 7L + 1L], levels = STM_WEEKDAYS)]
  R[, `:=`(prev_spell_l = fifelse(any_prior, log2(1 + prev_spell), 0), any_prior_spell = as.integer(any_prior), dsr_l = log2(1 + dsr))]
  sc <- list(sd_w = sqrt(mean((R$x - R$m_before)^2)), mu_b = mean(R$m_before), sd_b = stats::sd(R$m_before))
  R[, `:=`(w = (x - m_before) / sc$sd_w, b = (m_before - sc$mu_b) / sc$sd_b)]
  bs <- unclass(splines::ns(R$study_day, df = 3L)); R[, c("ns1", "ns2", "ns3") := list(bs[, 1L], bs[, 2L], bs[, 3L])]
  R[, `:=`(pw20 = (wear >= 1200 & n_wear >= 1200) %in% TRUE, pw22 = (wear >= 1320 & n_wear >= 1320) %in% TRUE)]
  setattr(R, "scales", data.table(sd_w = sc$sd_w, mu_b = sc$mu_b, sd_b = sc$sd_b, n_rows = nrow(R), n_people = uniqueN(R$person_id)))
  R[]
}

# Transitions to each category, overall and by cohort, against the event-count thresholds (decided before any association is looked at).
stm_event_counts <- function(R, scfg) {
  one <- function(d, label) { n <- d[, .N, by = outcome]; ev <- function(k) { v <- n[outcome == k, N]; if (length(v)) v else 0L }
    data.table(group = label, n_rows = nrow(d), n_people = uniqueN(d$person_id), category = 0:4, label = STM_LABELS[as.character(0:4)], events = vapply(0:4, ev, integer(1)))[
      , share_pct := 100 * events / n_rows] }
  t <- rbind(one(R, "all"), rbindlist(lapply(levels(R$cohort), function(co) one(R[cohort == co], co))))
  t[, `:=`(effect_claim_allowed = category > 0L & events >= scfg$min_events_effect, equivalence_claim_allowed = category > 0L & events >= scfg$min_events_equiv)][]
}

# What category 4 contains: the original patterns (activity, heart rate, sleep), overall and by cohort.
stm_cat4_composition <- function(R) {
  d <- R[outcome == 4L]; if (!nrow(d)) return(data.table())
  rbind(d[, .(group = "all", n = .N), by = .(pattern = n_status)], d[, .(n = .N), by = .(group = as.character(cohort), pattern = n_status)])[
    , share_pct := 100 * n / sum(n), by = group][order(group, -n)]
}

# ---- formula and fits ---------------------------------------------------------------------------------------------------------------------------

# Fixed terms: w, b and the covariates that vary in the data.
stm_terms <- function(dat) {
  vary <- function(a) a %in% names(dat) && uniqueN(dat[[a]][!is.na(dat[[a]])]) > 1L
  c("w", "b", Filter(vary, STM_COVARIATES))
}

# The rows of one cause-specific model: outcome 0 against category k (k = 1-4), or against sleep-missing nights split by heart-rate wear (kind "sub": cutoff 20 / 22, worn TRUE / FALSE).
stm_cause_data <- function(R, spec) {
  if (spec$kind == "cat") return(droplevels(R[outcome %in% c(0L, spec$k)][, y := as.integer(outcome == spec$k)]))
  pw <- R[[paste0("pw", spec$cut)]]
  droplevels(R[outcome == 0L | (outcome == 3L & pw == spec$worn)][, y := as.integer(outcome == 3L)])
}

# Model A: one cause-specific mixed model (random person intercept, the pipeline's engine settings) and everything derived from it.
stm_fit_cause <- function(dat, spec, pcfg, scfg, seed) {
  t0 <- Sys.time(); rhs <- stm_terms(dat)
  info <- data.table(model = spec$id, n_rows = nrow(dat), n_people = uniqueN(dat$person_id), n_events = sum(dat$y), status = "ok", singular = NA, sigma = NA_real_, note = "")
  if (info$n_events < 2L) return(list(info = info[, `:=`(status = "not_estimable", note = "no events")]))
  fit <- tryCatch(suppressMessages(fit_binary_mixed(paste0("y ~ ", paste(rhs, collapse = " + "), " + (1 | person_id)"), dat, pcfg)), error = function(e) e)
  if (inherits(fit, "error")) return(list(info = info[, `:=`(status = "failed", note = conditionMessage(fit))]))
  info[, `:=`(status = if (fit_converged(fit)) "ok" else "not_converged", singular = lme4::isSingular(fit), sigma = sqrt(as.numeric(lme4::VarCorr(fit)$person_id[1L, 1L])))]
  parts <- list(beta = lme4::fixef(fit), V = as.matrix(stats::vcov(fit)), sigma = info$sigma, rhs = rhs)
  coefs <- stm_coef_table(parts, spec)
  c(list(info = info[, elapsed_sec := as.numeric(difftime(Sys.time(), t0, units = "secs"))], coef = coefs), stm_absolute_risk(parts, dat, spec, scfg, pcfg, seed))
}

# Coefficients with odds ratios and Wald 95% and 90% intervals (the 90% interval is the one the equivalence test uses).
stm_coef_table <- function(parts, spec) {
  se <- sqrt(diag(parts$V)); b <- parts$beta
  data.table(model = spec$id, term = names(b), estimate = b, se = se, z = b / se, p = 2 * stats::pnorm(-abs(b / se)), or = exp(b), or_low = exp(b - 1.96 * se), or_high = exp(b + 1.96 * se),
             or90_low = exp(b - stats::qnorm(0.95) * se), or90_high = exp(b + stats::qnorm(0.95) * se))
}

# Summary of a probability matrix P (rows: w_grid, columns: point estimate then draws) into the curve and the risk difference (+1 SD against -1 SD).
stm_curve_rows <- function(P, model, source, scfg) {
  q <- function(M, p) apply(M, 1L, stats::quantile, p, names = FALSE)
  rows <- function(M) data.table(prob_pct = M[, 1L], low = q(M[, -1L, drop = FALSE], 0.025), high = q(M[, -1L, drop = FALSE], 0.975))
  curve <- cbind(data.table(model = model, source = source, w = scfg$w_grid), rows(P))
  i1 <- which(scfg$w_grid == 1); im <- which(scfg$w_grid == -1)
  rd <- cbind(data.table(model = model, source = source), rows(P[i1, , drop = FALSE] - P[im, , drop = FALSE]))
  setnames(rd, c("prob_pct", "low", "high"), c("risk_diff_pct", "rd_low", "rd_high"))
  list(curve = curve, rd = rd)
}

# Population-average probability of the category at w values, averaged over a seeded sample of rows (covariates as observed) and over the random-intercept distribution
# (quantile grid), with intervals from draws of the fixed effects (the random-effect SD held fixed).
stm_absolute_risk <- function(parts, dat, spec, scfg, pcfg, seed) {
  stopifnot(all(c(-1, 1) %in% scfg$w_grid))
  sm <- as.data.frame(dat[with_seed(seed, sort(sample.int(nrow(dat), min(scfg$n_profile_rows, nrow(dat)))))])
  Bm <- cbind(parts$beta, t(with_seed(seed + 1L, wr_mvn(scfg$n_beta_draws, parts$beta, parts$V))))
  u <- stats::qnorm(stats::ppoints(scfg$n_re_grid)) * parts$sigma
  at <- function(wv) { s <- sm; s$w <- wv
    lin <- stats::model.matrix(stats::reformulate(parts$rhs), s)[, names(parts$beta), drop = FALSE] %*% Bm
    Reduce(`+`, lapply(u, function(uj) colMeans(stats::plogis(lin + uj)))) / length(u) }
  P <- 100 * t(vapply(scfg$w_grid, at, numeric(ncol(Bm))))
  s <- stm_curve_rows(P, spec$id, "mixed", scfg)
  list(risk_curve = s$curve, risk_diff = s$rd)
}

# ---- tests ---------------------------------------------------------------------------------------------------------------------------------------------

# H1 (two-sided Wald test of w, category 1; no adjustment: the single primary test) and the equivalence tests (two one-sided tests against the margin; 90% interval within the
# margin) with the Holm adjustment across the equivalence tests of categories 2-4 (the unadjusted p stays), the claims allowed by the event counts and the verdict of the negative
# control (category 4): pass = 90% interval within the margin (and enough events), fail = whole 90% interval outside the margin, otherwise inconclusive.
# coef: stacked coefficient tables of the cause-specific fits (term "w"); info: their info rows.
stm_tests <- function(coef, info, scfg) {
  w <- merge(coef[term == "w"], info[, .(model, n_events)], by = "model")
  lo <- log(scfg$equiv_margin[1L]); hi <- log(scfg$equiv_margin[2L])
  w[, `:=`(p_wald = p, p_lower = 1 - stats::pnorm((estimate - lo) / se), p_upper = stats::pnorm((estimate - hi) / se))]
  w[, `:=`(p_tost = pmax(p_lower, p_upper), equivalent = or90_low >= scfg$equiv_margin[1L] & or90_high <= scfg$equiv_margin[2L],
           effect_claim_allowed = n_events >= scfg$min_events_effect, equivalence_claim_allowed = n_events >= scfg$min_events_equiv)]
  w[, c("p_lower", "p_upper") := NULL]
  eq <- w$model %in% c("c2", "c3", "c4"); w[, p_tost_holm := NA_real_]; w[eq, p_tost_holm := stats::p.adjust(p_tost, "holm")]
  w[, equivalent_holm := p_tost_holm < 0.05]
  w[, negative_control := fcase(model != "c4", NA_character_, equivalent & equivalence_claim_allowed, "pass", or90_low > scfg$equiv_margin[2L] | or90_high < scfg$equiv_margin[1L], "fail", default = "inconclusive")][]
}

# ---- model B: multinomial logit with person-clustered standard errors -------------------------------------------------------------------------

# Multinomial logit (reference outcome 0, categories 1-4) by Newton-Raphson on all rows. Returns the coefficients (p x 4), the inverse information (the bread), the per-person scores
# and the person-clustered covariance of the stacked coefficient vector (category blocks; HC0 with the g / (g - 1) factor, as the pooled logistic fits).
stm_multinom <- function(dat, rhs, max_iter = 30L, tol = 1e-8) {
  X <- stats::model.matrix(stats::reformulate(rhs), dat); K <- 4L
  qq <- qr(crossprod(X)); if (qq$rank < ncol(X)) X <- X[, sort(qq$pivot[seq_len(qq$rank)]), drop = FALSE]
  p <- ncol(X); Y <- vapply(seq_len(K), function(k) as.numeric(dat$outcome == k), numeric(nrow(X)))
  probs <- function(B) { E <- pmin(X %*% B, 700); mx <- pmax(0, apply(E, 1L, max)); ex <- exp(E - mx); den <- exp(-mx) + rowSums(ex); list(P = ex / den, ll = sum(Y * E) - sum(mx + log(den))) }
  info_mat <- function(P) { I <- matrix(0, K * p, K * p)
    for (j in seq_len(K)) for (k in j:K) {
      blk <- crossprod(X, X * (if (j == k) P[, j] * (1 - P[, j]) else -P[, j] * P[, k]))
      I[(j - 1L) * p + seq_len(p), (k - 1L) * p + seq_len(p)] <- blk; if (j != k) I[(k - 1L) * p + seq_len(p), (j - 1L) * p + seq_len(p)] <- t(blk)
    }
    I }
  B <- matrix(0, p, K, dimnames = list(colnames(X), seq_len(K))); cur <- probs(B); conv <- FALSE; it <- 0L
  repeat {
    it <- it + 1L; G <- crossprod(X, Y - cur$P)
    step <- matrix(solve(info_mat(cur$P), as.vector(G)), p, K); s <- 1
    repeat { Bn <- B + s * step; new <- probs(Bn); if (is.finite(new$ll) && new$ll >= cur$ll - 1e-10) break; s <- s / 2; if (s < 1e-6) break }
    done <- abs(new$ll - cur$ll) < tol * (abs(cur$ll) + 1) || max(abs(s * step)) < 1e-7
    B <- Bn; cur <- new; if (done) { conv <- TRUE; break }; if (it >= max_iter) break
  }
  P <- cur$P; bread <- solve(info_mat(P))
  S <- do.call(cbind, lapply(seq_len(K), function(k) rowsum(X * (Y[, k] - P[, k]), dat$person_id))); g <- nrow(S)
  V <- (g / (g - 1)) * bread %*% crossprod(S) %*% bread
  nm <- as.vector(outer(colnames(X), seq_len(K), function(a, k) paste0(k, ":", a))); dimnames(V) <- list(nm, nm)
  list(beta = B, V = V, loglik = cur$ll, iterations = it, converged = conv, n_rows = nrow(X), n_people = g, X_names = colnames(X), rhs = rhs)
}

# Contrasts of the w coefficient between two categories with the joint cluster covariance; Holm-adjusted p and Holm-width (Bonferroni) intervals over the contrasts.
# pairs: list of c(first, second) category numbers; ratio_of_or = ratio of the two odds ratios per +1 SD (population-averaged scale).
stm_contrasts <- function(fit, pairs, term = "w") {
  idx <- function(k) match(paste0(k, ":", term), rownames(fit$V))
  out <- rbindlist(lapply(pairs, function(pr) {
    i <- idx(pr[1L]); j <- idx(pr[2L]); d <- fit$beta[term, pr[1L]] - fit$beta[term, pr[2L]]; se <- sqrt(fit$V[i, i] + fit$V[j, j] - 2 * fit$V[i, j])
    data.table(contrast = sprintf("category %d vs category %d", pr[1L], pr[2L]), difference = d, se = se, ratio_of_or = exp(d), ratio_low = exp(d - 1.96 * se), ratio_high = exp(d + 1.96 * se),
               z = d / se, p = 2 * stats::pnorm(-abs(d / se)))
  }))
  m <- nrow(out); out[, `:=`(p_holm = stats::p.adjust(p, "holm"), ratio_low_holm = exp(difference - stats::qnorm(1 - 0.025 / m) * se), ratio_high_holm = exp(difference + stats::qnorm(1 - 0.025 / m) * se))][]
}

# Predicted probability of every category (0-4) across w from the multinomial fit: covariates as observed in a seeded sample of rows, intervals from draws of the coefficients
# with their cluster covariance.
stm_multi_risk <- function(fit, dat, scfg, seed) {
  stopifnot(all(c(-1, 1) %in% scfg$w_grid)); K <- 4L; p <- length(fit$X_names)
  sm <- as.data.frame(dat[with_seed(seed, sort(sample.int(nrow(dat), min(scfg$n_profile_rows, nrow(dat)))))])
  Bv <- cbind(as.vector(fit$beta), t(with_seed(seed + 1L, wr_mvn(scfg$n_beta_draws, as.vector(fit$beta), fit$V)))); nd <- ncol(Bv)
  Big <- matrix(0, p, K * nd); for (k in seq_len(K)) Big[, (seq_len(nd) - 1L) * K + k] <- Bv[(k - 1L) * p + seq_len(p), ]
  rhs <- fit$rhs
  at <- function(wv) { s <- sm; s$w <- wv; Xs <- stats::model.matrix(stats::reformulate(rhs), s)[, fit$X_names, drop = FALSE]; ex <- exp(pmin(Xs %*% Big, 700))
    ix <- lapply(seq_len(K), function(k) (seq_len(nd) - 1L) * K + k); den <- 1 + Reduce(`+`, lapply(ix, function(i) ex[, i]))
    rbind(colMeans(1 / den), do.call(rbind, lapply(ix, function(i) colMeans(ex[, i] / den)))) }       # 5 x nd
  A <- lapply(scfg$w_grid, at)
  parts <- lapply(0:4, function(k) stm_curve_rows(100 * t(vapply(A, function(m) m[k + 1L, ], numeric(nd))), paste0("c", k), "multinomial", scfg))
  list(risk_curve = rbindlist(lapply(parts, `[[`, "curve")), risk_diff = rbindlist(lapply(parts, `[[`, "rd")))
}

# Model B: fit, coefficient table (all terms, cluster SE), H2 contrasts, probabilities.
stm_fit_multinomial <- function(dat, pcfg, scfg, seed) {
  t0 <- Sys.time(); rhs <- stm_terms(dat)
  fit <- stm_multinom(dat, rhs, scfg$multinom_max_iter)
  se <- sqrt(diag(fit$V)); b <- as.vector(fit$beta)
  coef <- data.table(model = paste0("c", sub(":.*", "", rownames(fit$V))), term = sub("^[0-9]:", "", rownames(fit$V)), estimate = b, se = se, z = b / se, p = 2 * stats::pnorm(-abs(b / se)),
                     or = exp(b), or_low = exp(b - 1.96 * se), or_high = exp(b + 1.96 * se), or90_low = exp(b - stats::qnorm(0.95) * se), or90_high = exp(b + stats::qnorm(0.95) * se))
  info <- data.table(model = "multinomial", n_rows = fit$n_rows, n_people = fit$n_people, n_events = sum(dat$outcome > 0L), status = if (fit$converged) "ok" else "not_converged",
                     iterations = fit$iterations, loglik = fit$loglik, elapsed_sec = as.numeric(difftime(Sys.time(), t0, units = "secs")))
  c(list(info = info, coef = coef, contrasts = stm_contrasts(fit, list(c(1L, 2L), c(1L, 3L)))), stm_multi_risk(fit, dat, scfg, seed))
}

# ---- sensitivity, checks --------------------------------------------------------------------------------------------------------------------------------

# Rows of a sensitivity variant: first 14 days of follow-up excluded / only each person's first transition (rows up to and including the first outcome other than 0).
# (The wear thresholds of day d are variants of the DAY table: stm_days(D, stm_config(list(wear_hours = 16 / 22))).)
stm_variant_rows <- function(R, variant, scfg) {
  switch(variant, first14 = R[study_day > scfg$first_days_excluded],
         first_transition = R[, .SD[cumsum(shift(outcome != 0L, fill = FALSE)) == 0L], by = person_id], R)
}

# Conditional logistic regression stratified by person for category 1 (time-varying terms only).
stm_clogit <- function(R, scfg) {
  d <- R[outcome %in% c(0L, 1L)][, y := as.integer(outcome == 1L)]
  d <- d[person_id %in% d[y == 1L, unique(person_id)]]
  tv <- Filter(function(a) a %in% names(d) && uniqueN(d[[a]]) > 1L, c("w", "dow_next", "prev_spell_l", "any_prior_spell", "dsr_l", "ns1", "ns2", "ns3"))
  fit <- survival::clogit(stats::as.formula(paste("y ~", paste(tv, collapse = " + "), "+ strata(person_id)")), data = d, method = "efron")
  cf <- summary(fit)$coefficients["w", ]
  data.table(model = "clogit_c1", n_rows = nrow(d), n_people = uniqueN(d$person_id), n_events = sum(d$y), estimate = cf[["coef"]], se = cf[["se(coef)"]], or = exp(cf[["coef"]]),
             or_low = exp(cf[["coef"]] - 1.96 * cf[["se(coef)"]]), or_high = exp(cf[["coef"]] + 1.96 * cf[["se(coef)"]]))
}

# Agreement of lme4 (nAGQ = 0) with glmmTMB for category 1 on a seeded random subsample of people.
stm_verify_engine <- function(R, pcfg, scfg) {
  ids <- with_seed(scfg$seed, sample(unique(R$person_id), min(scfg$verify_people, uniqueN(R$person_id))))
  dat <- stm_cause_data(R[person_id %in% ids], list(kind = "cat", k = 1L, id = "c1")); f <- paste0("y ~ ", paste(stm_terms(dat), collapse = " + "), " + (1 | person_id)")
  est <- function(engine) { c2 <- pcfg; c2$model$engine <- engine; fit <- fit_binary_mixed(f, dat, c2); tidy_fixed(fit)[term == "w", estimate] }
  e <- tryCatch(c(lme4 = est("lme4"), glmmTMB = est("glmmTMB")), error = function(e) c(lme4 = NA_real_, glmmTMB = NA_real_))
  data.table(n_people = length(ids), n_events = sum(dat$y), w_lme4 = e[["lme4"]], w_glmmTMB = e[["glmmTMB"]], difference = e[["lme4"]] - e[["glmmTMB"]])
}

# Model-free: the share of risk-set rows followed by each category, by decile of w, with person-clustered 95% intervals.
stm_deciles <- function(R, scfg) {
  R <- copy(R)[, decile := cut(w, stats::quantile(w, seq(0, 1, length.out = scfg$deciles + 1L)), labels = FALSE, include.lowest = TRUE)]
  rbindlist(lapply(1:4, function(k) R[, { ci <- cluster_mean_ci(as.numeric(outcome == k), person_id)
    .(category = k, label = STM_LABELS[[as.character(k)]], mean_w = mean(w), n_rows = .N, share_pct = 100 * ci$mean, low = pmax(0, 100 * ci$lo), high = 100 * ci$hi) }, by = decile]))
}

# ---- one run: units with checkpoints --------------------------------------------------------------------------------------------------------------------

# The cause-specific models of a risk set: categories 1-4 and (main run only) the secondary split of category 3 at 20 and 22 hours.
stm_specs <- function(secondary = TRUE) {
  s <- lapply(1:4, function(k) list(id = paste0("c", k), kind = "cat", k = k))
  if (secondary) for (ct in c(20, 22)) for (wn in c(TRUE, FALSE))
    s[[length(s) + 1L]] <- list(id = sprintf("c3_%s_%d", if (wn) "probably_worn" else "probably_not_worn", ct), kind = "sub", k = 3L, cut = ct, worn = wn)
  s
}

# Fit (or resume) every model of a risk set; returns the unit results. variant: label of the run (part of every checkpoint key); dkey: fingerprint of the data.
# multinomial = FALSE skips model B (for the sensitivity variants that only need model A).
stm_run_models <- function(R, variant, pcfg, scfg, ckpt_dir, dkey, secondary = TRUE, multinomial = TRUE) {
  unit <- function(id, compute) run_unit(ckpt_dir, list(kind = "streammodel", stream = variant, id = id), ckpt_key(pcfg, "streammodel", dkey, list(id = id, variant = variant, scfg = scfg, n = nrow(R))), compute, pcfg)
  out <- lapply(stm_specs(secondary), function(sp) unit(sp$id, function() {
    dat <- stm_cause_data(R, sp); wr_check_memory(nrow(dat), 40L, pcfg)
    isolated(pcfg, stm_fit_cause(dat, sp, pcfg, scfg, scfg$seed + match(sp$id, vapply(stm_specs(), `[[`, "", "id"))))
  })$value)
  names(out) <- vapply(stm_specs(secondary), `[[`, "", "id")
  if (multinomial) out$multinomial <- unit("multinomial", function() {
    dat <- R[, c("outcome", "person_id", stm_terms(R)), with = FALSE]; wr_check_memory(nrow(dat), 120L, pcfg)
    isolated(pcfg, stm_fit_multinomial(dat, pcfg, scfg, scfg$seed + 100L))
  })$value
  out
}

# Collect the unit results into tables (variant column added).
stm_collect <- function(out, scfg, variant) {
  mix <- Filter(function(v) !is.null(v$info) && !identical(v$info$model, "multinomial"), out); mn <- out$multinomial
  bind <- function(nm, from = mix) rbindlist(lapply(from, `[[`, nm), fill = TRUE)[, variant := variant]
  info <- rbind(bind("info"), if (!is.null(mn)) bind("info", list(mn)), fill = TRUE)
  coef <- bind("coef")
  tests_models <- c("c1", "c2", "c3", "c4", "c3_probably_worn_20", "c3_probably_not_worn_20", "c3_probably_worn_22", "c3_probably_not_worn_22")
  list(info = info, coefficients = coef,
       tests = if (nrow(coef)) stm_tests(coef[model %in% tests_models], info, scfg)[, variant := variant] else data.table(),
       multinomial_coefficients = if (!is.null(mn)) bind("coef", list(mn)) else data.table(),
       contrasts = if (!is.null(mn)) bind("contrasts", list(mn)) else data.table(),
       risk_curve = rbind(bind("risk_curve"), if (!is.null(mn)) bind("risk_curve", list(mn)), fill = TRUE),
       risk_diff = rbind(bind("risk_diff"), if (!is.null(mn)) bind("risk_diff", list(mn)), fill = TRUE))
}
