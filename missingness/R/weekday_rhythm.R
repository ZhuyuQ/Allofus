# ==============================================================================
# weekday_rhythm.R
# ------------------------------------------------------------------------------
# Does the watch come off (gap start) or go back on (gap end, return) on particular weekdays? Is that a shared calendar routine or a
# charging cycle (days since return)? Does it differ between people (personal rhythm) and between cohorts?
#
# Terms (used exactly in comments, tables and captions)
#   observed day, gap, gap length    as in gaps.R
#   gap start             the first missing day of a gap: a gap starts on day d + 1 when d is observed and d + 1 is missing
#   gap end (return)      the first observed day after a gap
#   days so far in the gap  for a day t at or inside a gap that started on s: t - s
#   days since return     for an observed day d: d - r + 1 with r the last return (the return day counts as 1); before a person's first gap:
#                         d - first record + 1. Levels 1, ..., 21, 22+
#   weekday               Monday (0) ... Sunday (6) of the day that carries the outcome; the SLEEP stream is aligned to the night it STARTED
#                         (stored date - 1) in wr_align_date(), the only place that shift is made
#   weekly rhythm         a 7-day repeating pattern in the chance of a gap start or gap end
#   shared rhythm         the weekly pattern common to everyone (fixed effects)
#   personal rhythm       how each person's weekly pattern differs from the shared one (person random effects, first harmonic)
#
# Models are binomial mixed models (lme4::glmer, the pipeline's engine settings), one checkpointed unit per fit; see wr_specs().
# Rows: gap-start rows (every observed day except a person's last; no baseline and no record on d - 1 are required) and gap-end rows (one row per
# day a gap could end). Dates are used as stored (local time is not known).
# ==============================================================================

WR_WEEKDAYS <- c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")
WR_RE_TEXT  <- c(correlated = "(1 + c1 + s1 | person_id)",
                 independent_harmonics = "(1 | person_id) + (0 + c1 | person_id) + (0 + s1 | person_id)",
                 intercept = "(1 | person_id)")
WR_ADJUSTERS <- list(start = c("cohort", "age10", "sex2", "ns1", "ns2", "ns3", "device", "cci", "cci_missing", "prev_gap_len_l", "any_prior_gap"),
                     end   = c("cohort", "age10", "sex2", "ns1", "ns2", "ns3", "device", "cci", "cci_missing"))

# ---- calendar -------------------------------------------------------------------------------------------------------------

# The date whose weekday and holiday status a row uses. Sleep dates are wake-up dates, so the night started the day before.
wr_align_date <- function(date, stream, cfg) {
  if (identical(stream, "sleep") && identical(cfg$weekday_rhythm$sleep_alignment, "night_start")) date - 1L else date
}

# Weekday index of the aligned date: 0 = Monday ... 6 = Sunday.
wr_weekday <- function(date, stream, cfg) as.integer((data.table::wday(wr_align_date(date, stream, cfg)) + 5L) %% 7L)

# The n-th (1-based; -1 = last) given weekday (0 = Monday) of a month.
wr_nth_weekday <- function(year, month, weekday, n) {
  first <- as.Date(sprintf("%d-%02d-01", year, month))
  if (n > 0L) return(first + ((weekday - (data.table::wday(first) + 5L) %% 7L) %% 7L) + 7L * (n - 1L))
  last <- seq(first, by = "month", length.out = 2L)[2L] - 1L
  last - (((data.table::wday(last) + 5L) %% 7L - weekday) %% 7L)
}

# US federal holidays of the given years: the actual date and, for fixed-date holidays on a weekend, the observed date (Saturday -> Friday,
# Sunday -> Monday). Built by rule; no package and no network. Juneteenth from 2021.
wr_us_holidays <- function(years) {
  fixed <- data.table(holiday = c("New Year's Day", "Juneteenth", "Independence Day", "Veterans Day", "Christmas Day"),
                      month = c(1L, 6L, 7L, 11L, 12L), day = c(1L, 19L, 4L, 11L, 25L))
  out <- rbindlist(lapply(years, function(y) {
    f <- copy(fixed)[holiday != "Juneteenth" | y >= 2021L]
    f[, date := as.Date(sprintf("%d-%02d-%02d", y, month, day))]
    f <- f[, .(holiday, date, observed = date + fcase(data.table::wday(date) == 7L, -1L, data.table::wday(date) == 1L, 1L, default = 0L))]
    moving <- data.table(holiday = c("Martin Luther King Jr. Day", "Presidents Day", "Memorial Day", "Labor Day", "Columbus Day", "Thanksgiving Day"),
                         month = c(1L, 2L, 5L, 9L, 10L, 11L), weekday = c(0L, 0L, 0L, 0L, 0L, 3L), nth = c(3L, 3L, -1L, 1L, 2L, 4L))
    m <- moving[, .(holiday, date = as.Date(mapply(function(mo, wd, n) wr_nth_weekday(y, mo, wd, n), month, weekday, nth), origin = "1970-01-01"))]
    rbind(f, m[, observed := date])
  }))
  setorder(out, date)
  out[]
}

# TRUE where a date lies within +/- `window` days of any holiday date (vectorised).
wr_near_holiday <- function(dates, holiday_dates, window) {
  h <- sort(unique(as.integer(holiday_dates))); x <- as.integer(dates)
  i <- findInterval(x, h)
  d_lo <- ifelse(i >= 1L, x - h[pmax(i, 1L)], Inf)
  d_hi <- ifelse(i < length(h), h[pmin(i + 1L, length(h))] - x, Inf)
  pmin(d_lo, d_hi) <= window
}

# All dates of a holiday lookup: actual and observed dates of the years around the data.
wr_holiday_dates <- function(dates) {
  yrs <- seq(min(data.table::year(dates)) - 1L, max(data.table::year(dates)) + 1L)
  h <- wr_us_holidays(yrs)
  list(table = h, dates = unique(c(h$date, h$observed)))
}

# ---- rows -----------------------------------------------------------------------------------------------------------------

# Person covariates joined onto a row table (the same covariates as the primary model).
wr_join_covariates <- function(rows, persons, cfg) {
  cov <- person_covariates(persons, cfg)
  setkey(cov, person_id)
  cov[rows, on = "person_id"]
}

# Weekday, first-harmonic columns and the days-since-return / days-so-far factors of a row table.
wr_add_design <- function(rows, stream, cfg) {
  wd <- wr_weekday(rows$out_date, stream, cfg)
  w <- factor(WR_WEEKDAYS[wd + 1L], levels = WR_WEEKDAYS)
  stats::contrasts(w) <- stats::contr.sum(7L)
  set(rows, j = c("wd", "c1", "s1", "weekday"), value = list(wd, cos(2 * pi * wd / 7), sin(2 * pi * wd / 7), w))
  rows
}

# Natural spline of study day (3 df) as plain columns ns1-ns3 (computed on the rows of the fit it is used in).
wr_add_spline <- function(rows) {
  b <- unclass(splines::ns(rows$study_day, df = 3L))
  set(rows, j = c("ns1", "ns2", "ns3"), value = list(b[, 1L], b[, 2L], b[, 3L]))
  rows
}

# Gap-start rows: one per observed day d except each person's last; outcome start = 1 if a gap starts on d + 1.
# attr(, "flow") counts what was left out.
wr_gap_start_rows <- function(days, persons, stream, cfg) {
  rs <- build_risk_set(days, cfg)
  rows <- rs[, .(person_id, date, out_date = date + 1L, start = as.integer(gap_len_next >= 1L), gap_length = gap_len_next, study_day,
                 any_prior_gap, prev_gap_len_l = log2(1 + prev_gap_len),
                 dsr = fifelse(any_prior_gap, days_since_return + 1L, study_day))]
  set(rows, j = "dsr_f", value = factor(pmin(rows$dsr, cfg$weekday_rhythm$dsr_max + 1L), levels = seq_len(cfg$weekday_rhythm$dsr_max + 1L),
                                       labels = c(as.character(seq_len(cfg$weekday_rhythm$dsr_max)), paste0(cfg$weekday_rhythm$dsr_max + 1L, "+"))))
  rows <- wr_add_design(wr_join_covariates(rows, persons, cfg), stream, cfg)
  setattr(rows, "flow", data.table(rows_type = "gap start", step = c("observed days with a next observed day or gap", "left out by construction: each person's last observed day"),
                                   n_rows = c(nrow(rows), uniqueN(days$person_id)), n_people = c(uniqueN(rows$person_id), NA_integer_)))
  rows[]
}

# Gap-end rows: for a gap of length L starting on s, the days t = s + 1, ..., s + L (days so far in the gap = t - s); outcome end = 1 only on the return
# day t = s + L. Rows with more than days_so_far_max days so far (14 always end, because the longest gap is 14 days) are dropped and counted.
wr_gap_end_rows <- function(days, persons, stream, cfg) {
  gaps <- get_gap_runs(days)
  idx <- rep.int(seq_len(nrow(gaps)), gaps$gap_length); dsf <- sequence(gaps$gap_length)
  rows <- data.table(person_id = gaps$person_id[idx], dsf = dsf, gap_length = gaps$gap_length[idx])
  rows[, `:=`(out_date = gaps$gap_start[idx] + dsf, end = as.integer(dsf == gap_length))]
  n_all <- nrow(rows); rows <- rows[dsf <= cfg$weekday_rhythm$days_so_far_max]
  first <- days[, .(first_date = min(date)), by = person_id]
  rows[first, study_day := as.integer(out_date - i.first_date) + 1L, on = "person_id"]
  set(rows, j = "dsf_f", value = factor(rows$dsf, levels = seq_len(cfg$weekday_rhythm$days_so_far_max)))
  rows <- wr_add_design(wr_join_covariates(rows, persons, cfg), stream, cfg)
  setattr(rows, "flow", data.table(rows_type = "gap end", step = c("days in gaps (days so far 1 to gap length)", "dropped: more than days_so_far_max days so far (the gap must end)"),
                                   n_rows = c(n_all, n_all - nrow(rows)), n_people = c(uniqueN(gaps$person_id), NA_integer_)))
  rows[]
}

# ---- specification of the fits ---------------------------------------------------------------------------------------------

# One row per fit. rows: gap-start or gap-end rows; dsr: days since return (start) / days so far in the gap (end) as a factor; weekday: sum-to-zero factor;
# cohort_x: weekday:cohort; re: "auto" (chosen once per stream and row type by S1 / E1) or "intercept"; filter: the rows of the fit; activity: "predictors" adds
# the stream's primary w_ and b_ terms; curve / personal: derived outputs wanted from this fit.
wr_specs <- function() {
  s <- function(fit_id, rows, dsr = FALSE, weekday = FALSE, cohort_x = FALSE, re = "auto", filter = "none", activity = "none", curve = FALSE, used_for = "")
    data.table(fit_id, rows, dsr, weekday, cohort_x, re, filter, activity, curve, personal = fit_id == "S1", used_for)
  rbindlist(list(
    s("S0", "start", used_for = "null of H1 (starts)"),
    s("S1", "start", weekday = TRUE, used_for = "H1 (starts), H3, H4 null, Figures 1 and 3"),
    s("S2", "start", dsr = TRUE, curve = TRUE, used_for = "null of H2"),
    s("S3", "start", dsr = TRUE, weekday = TRUE, curve = TRUE, used_for = "H2, Figure 2"),
    s("S4", "start", weekday = TRUE, cohort_x = TRUE, used_for = "H4 (starts), Figure 4"),
    s("S5", "start", weekday = TRUE, re = "intercept", used_for = "null of H3"),
    s("E0", "end", dsr = TRUE, used_for = "null of H1 (ends)"),
    s("E1", "end", dsr = TRUE, weekday = TRUE, used_for = "H1 (ends), H4 null, Figure 1"),
    s("E4", "end", dsr = TRUE, weekday = TRUE, cohort_x = TRUE, used_for = "H4 (ends), Figure 4"),
    s("S1_holidays", "start", weekday = TRUE, filter = "holidays", used_for = "sensitivity: holidays"),
    s("S3_holidays", "start", dsr = TRUE, weekday = TRUE, filter = "holidays", used_for = "sensitivity: holidays"),
    s("E1_holidays", "end", dsr = TRUE, weekday = TRUE, filter = "holidays", used_for = "sensitivity: holidays"),
    s("S1_before_first_gap", "start", weekday = TRUE, filter = "first_gap", used_for = "sensitivity: before the first gap"),
    s("S3_before_first_gap", "start", dsr = TRUE, weekday = TRUE, filter = "first_gap", used_for = "sensitivity: before the first gap"),
    s("S1_activity_same_rows", "start", weekday = TRUE, filter = "primary_rows", used_for = "sensitivity: the primary model's rows, no predictors"),
    s("S1_activity", "start", weekday = TRUE, filter = "primary_rows", activity = "predictors", used_for = "sensitivity: the primary model's rows, with the predictors of d - 1")))
}

# The pairs of fits behind each likelihood-ratio test: id, family (Holm family), outcome, big, small.
wr_test_pairs <- function() {
  data.table(id = c("H1s", "H1e", "H2", "H2r", "H3", "H4s", "H4e"),
             family = c("H1", "H1", "H2", NA, "H3", "H4", "H4"),
             outcome = c("gap start", "gap end", "gap start", "gap start", "gap start", "gap start", "gap end"),
             big = c("S1", "E1", "S3", "S3", "S1", "S4", "E4"), small = c("S0", "E0", "S2", "S1", "S5", "S1", "E1"))
}

# Fixed-effect terms of a fit: the adjusters that exist and vary in the data, then days since return / days so far, weekday, weekday:cohort, predictors.
wr_terms <- function(spec, dat, predictors = character()) {
  adj <- WR_ADJUSTERS[[spec$rows]]
  adj <- adj[adj %in% names(dat) & vapply(adj, function(a) a %in% names(dat) && uniqueN(dat[[a]][!is.na(dat[[a]])]) > 1L, logical(1))]
  c(adj, if (spec$dsr) c(start = "dsr_f", end = "dsf_f")[[spec$rows]], if (spec$weekday) "weekday",
    if (spec$cohort_x && "cohort" %in% adj) "weekday:cohort", predictors)
}

# The formula of a fit: list(text, rhs) (random-effects structure from WR_RE_TEXT).
wr_formula <- function(spec, dat, re_structure, predictors = character()) {
  rhs <- wr_terms(spec, dat, predictors)
  re <- WR_RE_TEXT[[if (identical(spec$re, "intercept")) "intercept" else re_structure]]
  list(text = paste0(if (spec$rows == "start") "start" else "end", " ~ ", paste(rhs, collapse = " + "), " + ", re), rhs = rhs)
}

# ---- one fit ----------------------------------------------------------------------------------------------------------------

# Fit and collect the warnings; returns list(fit, notes).
wr_fit_core <- function(dat, formula_text, cfg) {
  notes <- character()
  fit <- tryCatch(withCallingHandlers(suppressMessages(fit_binary_mixed(formula_text, dat, cfg)),
                                      warning = function(w) { notes <<- c(notes, conditionMessage(w)); invokeRestart("muffleWarning") }),
                  error = function(e) e)
  list(fit = fit, notes = notes)
}

# Converged (and, if `strict`, not singular)?
wr_fit_ok <- function(fit, strict = FALSE) {
  !inherits(fit, "error") && fit_converged(fit) && !(strict && lme4::isSingular(fit))
}

# Random-effects covariance as a 3 x 3 matrix over (Intercept, c1, s1); absent terms are 0.
wr_re_cov <- function(fit) {
  nm <- c("(Intercept)", "c1", "s1")
  S <- matrix(0, 3L, 3L, dimnames = list(nm, nm))
  for (g in lme4::VarCorr(fit)) { m <- as.matrix(g); S[rownames(m), rownames(m)] <- m }
  S
}

# Draws from N(mean, V) by an eigen decomposition (robust to a slightly indefinite V).
wr_mvn <- function(n, mean, V) {
  ev <- eigen(V, symmetric = TRUE)
  A <- ev$vectors %*% diag(sqrt(pmax(ev$values, 0)), length(ev$values))
  sweep(matrix(stats::rnorm(n * length(mean)), n) %*% t(A), 2L, mean, "+")
}

# Fixed-effect table of a fit.
wr_coef_table <- function(fit) {
  cf <- stats::coef(summary(fit))
  data.table(term = rownames(cf), estimate = cf[, 1L], se = cf[, 2L], z = cf[, 3L], p = cf[, 4L])
}

# What the derived quantities need from a fit: coefficients, their covariance, the random-effects covariance and the fixed terms.
wr_fit_parts <- function(fit, rhs) {
  list(beta = lme4::fixef(fit), V = as.matrix(stats::vcov(fit)), Sigma = wr_re_cov(fit), rhs = rhs)
}

# ---- population-average predictions ------------------------------------------------------------------------------------------

# Predicted probability, averaged over the rows of `sm` and over the random effects, for every variant in `variants` (a list of data frames).
# Returns a matrix variants x (1 + draws): the first column at the estimated coefficients, the others at draws from their normal distribution
# (the random-effects covariance is held fixed). Vectorised over rows and coefficient draws; one loop over the shared random-effects draws.
wr_predict_matrix <- function(parts, variants, n_beta, n_re, seed) {
  draws <- with_seed(seed, list(b = wr_mvn(n_beta, parts$beta, parts$V), u = wr_mvn(n_re, c(0, 0, 0), parts$Sigma)))
  Bm <- cbind(parts$beta, t(draws$b))
  one <- function(sm) {
    X <- stats::model.matrix(stats::reformulate(parts$rhs), sm)[, names(parts$beta), drop = FALSE]
    lin <- X %*% Bm
    R <- cbind(1, sm$c1, sm$s1) %*% t(draws$u)
    acc <- numeric(ncol(Bm))
    for (j in seq_len(ncol(R))) acc <- acc + colMeans(stats::plogis(lin + R[, j]))
    acc / ncol(R)
  }
  t(vapply(variants, one, numeric(ncol(Bm))))
}

# A seeded sample of rows as a data frame.
wr_sample_rows <- function(dat, n, seed) {
  as.data.frame(dat[with_seed(seed, sort(sample.int(nrow(dat), min(n, nrow(dat)))))])
}

# The sample with the weekday set to each of the 7 days.
wr_weekday_variants <- function(sm) {
  lapply(0:6, function(k) {
    v <- sm
    w <- factor(rep(WR_WEEKDAYS[k + 1L], nrow(sm)), levels = WR_WEEKDAYS); stats::contrasts(w) <- stats::contr.sum(7L)
    v$weekday <- w; v$c1 <- cos(2 * pi * k / 7); v$s1 <- sin(2 * pi * k / 7)
    v
  })
}

# Profile and summary from the matrix of predicted probabilities (weekdays x (1 + draws)).
# ratio_to_average_day = p_k / mean(p); peak_to_trough = max p / min p: ratios of predicted probabilities, not odds ratios.
wr_profile_from_matrix <- function(P) {
  q <- function(x) stats::quantile(x, c(0.025, 0.975), names = FALSE)
  ratio <- sweep(P, 2L, colMeans(P), "/")
  ptr <- apply(P, 2L, max) / apply(P, 2L, min)
  prof <- data.table(weekday = 0:6, pct = 100 * P[, 1L], pct_low = 100 * apply(P[, -1L, drop = FALSE], 1L, function(x) q(x)[1L]),
                     pct_high = 100 * apply(P[, -1L, drop = FALSE], 1L, function(x) q(x)[2L]), ratio_to_average_day = ratio[, 1L],
                     ratio_low = apply(ratio[, -1L, drop = FALSE], 1L, function(x) q(x)[1L]), ratio_high = apply(ratio[, -1L, drop = FALSE], 1L, function(x) q(x)[2L]))
  summ <- data.table(peak_to_trough = ptr[1L], ptr_low = q(ptr[-1L])[1L], ptr_high = q(ptr[-1L])[2L], peak_weekday = which.max(P[, 1L]) - 1L,
                     trough_weekday = which.min(P[, 1L]) - 1L, average_pct = 100 * mean(P[, 1L]))
  list(profile = prof, summary = summ)
}

# Weekday profile (and summary) of the rows `dat` for one fit.
wr_profile_for <- function(parts, dat, cfg, seed) {
  sm <- wr_sample_rows(dat, cfg$weekday_rhythm$n_profile_rows, seed)
  wr_profile_from_matrix(wr_predict_matrix(parts, wr_weekday_variants(sm), cfg$weekday_rhythm$n_beta_draws, cfg$weekday_rhythm$n_re_draws, seed))
}

# Weekday profile of a fit: overall, or one per cohort for the weekday:cohort fits.
wr_weekday_outputs <- function(parts, dat, spec, cfg, seed) {
  if (!spec$cohort_x) { r <- wr_profile_for(parts, dat, cfg, seed); return(list(profile = r$profile, summary = r$summary)) }
  parts_c <- lapply(levels(dat$cohort), function(co) { r <- wr_profile_for(parts, dat[cohort == co], cfg, seed + match(co, levels(dat$cohort)))
                                                         list(profile = cbind(cohort = co, r$profile), summary = cbind(cohort = co, r$summary)) })
  list(profile = rbindlist(lapply(parts_c, `[[`, "profile")), summary = rbindlist(lapply(parts_c, `[[`, "summary")))
}

# Days since return curve: predicted % by level, averaged over the observed weekday mix (and the random effects).
wr_dsr_curve <- function(parts, dat, cfg, seed) {
  sm <- wr_sample_rows(dat, cfg$weekday_rhythm$n_profile_rows, seed)
  lv <- levels(dat$dsr_f)
  variants <- lapply(lv, function(l) { v <- sm; v$dsr_f <- factor(rep(l, nrow(sm)), levels = lv); v })
  P <- 100 * wr_predict_matrix(parts, variants, cfg$weekday_rhythm$n_curve_draws, cfg$weekday_rhythm$n_re_draws, seed)
  data.table(level = lv, pct = P[, 1L], pct_low = apply(P[, -1L, drop = FALSE], 1L, stats::quantile, 0.025, names = FALSE),
             pct_high = apply(P[, -1L, drop = FALSE], 1L, stats::quantile, 0.975, names = FALSE))
}

# ---- personal rhythm --------------------------------------------------------------------------------------------------------

# Fixed weekday effects (log-odds vs the average day, sum to zero) projected on the first harmonic: (a, b) with e_k ~ a cos + b sin.
wr_first_harmonic <- function(beta) {
  e <- c(beta[paste0("weekday", 1:6)], -sum(beta[paste0("weekday", 1:6)]))
  k <- 0:6
  c(a = (2 / 7) * sum(e * cos(2 * pi * k / 7)), b = (2 / 7) * sum(e * sin(2 * pi * k / 7)))
}

# Peak weekday (0-7, a position on the weekly circle) of a first-harmonic vector (a, b).
wr_peak_day <- function(a, b) (atan2(b, a) * 7 / (2 * pi)) %% 7

# Personal rhythm of a S1 fit. Each person's first-harmonic rhythm is (a + u_c, b + u_s); amplitude A, peak-to-trough exp(2 A) (odds scale).
# Reported from seeded draws of (u_c, u_s) from the estimated covariance (model based); BLUPs are only binned for Figure 3 (shrunk toward the shared rhythm).
wr_personal_outputs <- function(fit, parts, cfg, seed) {
  hb <- wr_first_harmonic(parts$beta); Sg <- parts$Sigma[c("c1", "s1"), c("c1", "s1")]
  u <- with_seed(seed + 7L, wr_mvn(cfg$weekday_rhythm$n_person_draws, c(0, 0), Sg))
  x <- hb[["a"]] + u[, 1L]; y <- hb[["b"]] + u[, 2L]
  amp <- sqrt(x^2 + y^2); peak <- wr_peak_day(x, y); shared <- wr_peak_day(hb[["a"]], hb[["b"]])
  dist <- abs(((peak - shared) + 3.5) %% 7 - 3.5)
  ptr <- exp(2 * amp)
  list(personal = data.table(sd_intercept = sqrt(parts$Sigma[1L, 1L]), sd_c1 = sqrt(Sg[1L, 1L]), sd_s1 = sqrt(Sg[2L, 2L]),
                             cor_c1_s1 = if (all(diag(Sg) > 0)) Sg[1L, 2L] / sqrt(Sg[1L, 1L] * Sg[2L, 2L]) else NA_real_,
                             shared_a = hb[["a"]], shared_b = hb[["b"]], shared_amplitude = sqrt(sum(hb^2)), shared_peak_day = shared,
                             person_ptr_median = stats::median(ptr), person_ptr_p10 = stats::quantile(ptr, 0.1, names = FALSE),
                             person_ptr_p90 = stats::quantile(ptr, 0.9, names = FALSE), share_near_shared_peak = mean(dist <= 1)),
       blup_bins = wr_blup_bins(fit, hb, cfg), ellipse = wr_ellipse(hb, Sg))
}

# BLUPs of (u_c, u_s) -> (a + u_c, b + u_s) in square bins; bins with fewer than min_people people are suppressed (n = NA).
wr_blup_bins <- function(fit, hb, cfg) {
  re <- lme4::ranef(fit)
  get <- function(nm) { for (g in re) if (nm %in% names(g)) return(g[[nm]]); rep(0, nrow(re[[1L]])) }
  x <- hb[["a"]] + get("c1"); y <- hb[["b"]] + get("s1")
  w <- max(abs(c(x, y))) / 12 + 1e-9
  b <- data.table(x = round(x / w) * w, y = round(y / w) * w)[, .(n = .N), by = .(x, y)]
  b[, `:=`(n_total = length(x), width = w)]
  b[n < cfg$weekday_rhythm$min_people, n := NA_integer_]
  b[]
}

# 95% region of the model-based personal rhythm (a + u_c, b + u_s).
wr_ellipse <- function(hb, Sg) {
  th <- seq(0, 2 * pi, length.out = 101L)
  R <- eigen(Sg, symmetric = TRUE)
  A <- R$vectors %*% diag(sqrt(pmax(R$values, 0)), 2L) * sqrt(stats::qchisq(0.95, 2))
  data.table(x = hb[["a"]] + A[1L, 1L] * cos(th) + A[1L, 2L] * sin(th), y = hb[["b"]] + A[2L, 1L] * cos(th) + A[2L, 2L] * sin(th))
}

# ---- one unit: filter the rows, fit, extract ----------------------------------------------------------------------------------

# Rows of a fit after its filter. `primary_rows` is a function returning the primary model's rows (only needed for those fits).
wr_apply_filter <- function(rows, spec, stream, cfg, primary_rows = NULL) {
  dat <- switch(spec$filter,
    none = copy(rows),
    holidays = rows[!wr_near_holiday(wr_align_date(out_date, stream, cfg), wr_holiday_dates(rows$out_date)$dates, cfg$weekday_rhythm$holiday_window)],
    first_gap = rows[any_prior_gap == FALSE],
    primary_rows = primary_rows(),
    stop("unknown row filter: ", spec$filter, call. = FALSE))
  preds <- attr(dat, "predictors")
  needed <- intersect(c(WR_ADJUSTERS[[spec$rows]], "weekday", "person_id"), names(dat))
  dat <- wr_add_spline(dat[stats::complete.cases(dat[, ..needed])])
  if (!is.null(preds)) setattr(dat, "predictors", preds)
  dat
}

# Check that a fit of this size fits in the free memory (same rule and message style as the other models).
wr_check_memory <- function(n_rows, n_cols, cfg) {
  if (!isTRUE(cfg$resources$memory_guard)) return(invisible(TRUE))
  need <- estimate_fit_memory_gb(n_rows, n_cols, cfg$resources$bytes_per_cell); head <- memory_headroom_gb(cfg)
  if (need > head && need <= cfg$resources$memory_budget_gb) head <- wait_for_headroom(cfg, need)
  if (need > head) stop(sprintf("this fit needs about %.0f GB but only %.0f GB can be used (free memory minus the %.0f GB reserve). Nothing was subsampled; restart the kernel and run only the setup and load cells, or set cfg$resources$memory_guard <- FALSE to try anyway.",
                                need, head, cfg$resources$memory_reserve_gb), call. = FALSE)
  invisible(TRUE)
}

# The empty result of a fit that was not run or failed.
wr_info_row <- function(spec, stream, dat, status, note = "", re = NA_character_) {
  data.table(stream = stream, fit_id = spec$fit_id, rows = spec$rows, n_rows = nrow(dat), n_people = uniqueN(dat$person_id),
             n_events = sum(dat[[if (spec$rows == "start") "start" else "end"]]), status = status, converged = NA, singular = NA,
             re_structure = re, loglik = NA_real_, df = NA_integer_, elapsed_sec = NA_real_, note = note)
}

# Enough events to fit (the pipeline's thresholds)?
wr_enough_events <- function(dat, spec, cfg) {
  y <- dat[[if (spec$rows == "start") "start" else "end"]]
  sum(y) >= cfg$model$min_events && uniqueN(dat$person_id[y == 1L]) >= cfg$model$min_persons_with_events
}

# The fit with the given random-effects structure; decide = TRUE (S1 / E1): try the correlated structure and, if it does not converge or is
# singular, independent harmonics. The choice is then used for every other fit of the stream and row type.
wr_fit_attempts <- function(dat, spec, cfg, re_structure, decide, predictors) {
  attempt <- function(re) { f <- wr_formula(spec, dat, re, predictors); c(wr_fit_core(dat, f$text, cfg), list(f = f, re = re)) }
  a <- attempt(if (decide) "correlated" else re_structure)
  if (decide && !wr_fit_ok(a$fit, strict = TRUE)) a <- attempt("independent_harmonics")
  a
}

# The info row of a fitted model (status, convergence, singularity, random-effects structure, log-likelihood).
wr_fit_info <- function(a, dat, spec, stream) {
  ok <- wr_fit_ok(a$fit)
  info <- wr_info_row(spec, stream, dat, if (ok) "ok" else "not_converged", paste(unique(a$notes), collapse = " | "), a$re)
  info[, `:=`(converged = ok, singular = lme4::isSingular(a$fit), loglik = as.numeric(stats::logLik(a$fit)), df = attr(stats::logLik(a$fit), "df"))]
  info[]
}

# Coefficients plus the derived outputs this fit is wanted for (weekday profile, days-since-return curve, personal rhythm); the wishes are data (wr_specs()).
wr_fit_outputs <- function(a, dat, spec, cfg, seed) {
  parts <- wr_fit_parts(a$fit, a$f$rhs)
  wanted <- list(list(spec$weekday, function() wr_weekday_outputs(parts, dat, spec, cfg, seed)),
                 list(spec$curve, function() list(curve = wr_dsr_curve(parts, dat, cfg, seed))),
                 list(spec$personal && spec$re != "intercept", function() wr_personal_outputs(a$fit, parts, cfg, seed)))
  c(list(coef = wr_coef_table(a$fit)), unlist(lapply(Filter(function(w) w[[1L]], wanted), function(w) w[[2L]]()), recursive = FALSE))
}

# Fit one specification and extract everything derived from it (the value of a unit).
wr_fit_value <- function(dat, spec, stream, cfg, re_structure, seed, decide = FALSE, predictors = character()) {
  t0 <- Sys.time()
  if (!wr_enough_events(dat, spec, cfg)) return(list(info = wr_info_row(spec, stream, dat, "not_estimable", "too few events")))
  a <- wr_fit_attempts(dat, spec, cfg, re_structure, decide, predictors)
  if (inherits(a$fit, "error")) return(list(info = wr_info_row(spec, stream, dat, "failed", conditionMessage(a$fit), a$re)))
  out <- c(list(info = wr_fit_info(a, dat, spec, stream)), wr_fit_outputs(a, dat, spec, cfg, seed))
  out$info[, elapsed_sec := as.numeric(difftime(Sys.time(), t0, units = "secs"))]
  out
}

# Seed of a fit: fixed by the configured seed, the stream and the fit id (not by the order the units run in).
wr_seed <- function(cfg, stream, fit_id) cfg$seed + 100000L + 1000L * match(stream, c("activity", "sleep", "hr")) + match(fit_id, wr_specs()$fit_id)

# The rows of the primary model (record on d - 1, complete predictors) in the format of the gap-start rows, and the predictor terms.
wr_primary_rows <- function(bundle, prim, stream, cfg) {
  spec <- as.list(prim); spec$predictors <- list(spec$predictors[[1L]])
  d <- build_onset_data(bundle, spec, cfg, list(keep_date = TRUE))
  d[, `:=`(start = y_ge1, out_date = date + 1L)]
  wr_add_design(d, stream, cfg)
  setattr(d, "predictors", as.vector(rbind(paste0("w_", spec$predictors[[1L]]), paste0("b_", spec$predictors[[1L]]))))
  d
}

# ---- the stage of one stream -------------------------------------------------------------------------------------------------

# The rows a fit uses and its fitted value; runs in a worker. `primary` caches the primary model's rows (only the activity fits need them).
wr_fit_compute <- function(sp, s, start_rows, end_rows, get_primary, re, decide, cfg) {
  rows <- if (sp$rows == "start") start_rows$get() else end_rows$get()
  dat <- wr_apply_filter(rows, sp, s, cfg, if (sp$filter == "primary_rows") get_primary else NULL)
  preds <- if (sp$activity == "predictors") attr(dat, "predictors") else character()
  wr_check_memory(nrow(dat), 30L + 21L * sp$dsr + 18L * sp$cohort_x + 3L, cfg)
  isolated(cfg, wr_fit_value(dat, sp, s, cfg, re, wr_seed(cfg, s, sp$fit_id), decide, preds))
}

# A function that runs (or resumes) one unit of the weekday stage; its status is recorded in `status` (an environment).
wr_unit_runner <- function(s, cfg, ckpt_dir, dkey, status) {
  function(id, compute, extra = NULL) {
    u <- run_unit(ckpt_dir, list(kind = "weekday", stream = s, id = id),
                  ckpt_key(cfg, "weekday", dkey, c(list(id = id, weekday_rhythm = cfg$weekday_rhythm), extra)), compute, cfg)
    assign(id, u$status, envir = status); u
  }
}

# Runs (or resumes) every fit of one stream as a checkpointed unit; the descriptions are the unit "DESC".
#   start_rows, end_rows, bundle: lazy objects; prim: the primary row of the onset grid; dkey: fingerprint of the stream's data
# S1 and E1 run first and decide the random-effects structure of their row type; every other fit follows it. Returns the status of each unit.
run_weekday_stage <- function(s, start_rows, end_rows, bundle, prim, cfg, ckpt_dir, dkey) {
  status <- new.env(); unit <- wr_unit_runner(s, cfg, ckpt_dir, dkey, status); specs <- wr_specs(); re <- list(); primary <- NULL
  get_primary <- function() { if (is.null(primary)) primary <<- wr_primary_rows(bundle$get(), prim, s, cfg); primary }
  fit_unit <- function(sp, decide) unit(sp$fit_id, function() wr_fit_compute(sp, s, start_rows, end_rows, get_primary, re[[sp$rows]], decide, cfg),
                                        list(spec = as.list(sp), re = re[[sp$rows]]))
  for (decider in c("S1", "E1")) re[specs[fit_id == decider]$rows] <- list(fit_unit(specs[fit_id == decider], TRUE)$value$info$re_structure)
  for (sp in split(specs[!fit_id %in% c("S1", "E1")], by = "fit_id")) {
    if (is.null(re[[sp$rows]])) assign(sp$fit_id, "skipped", envir = status) else fit_unit(sp, FALSE)
  }
  unit("DESC", function() wr_describe(start_rows$get(), end_rows$get(), s, cfg))
  unlist(as.list(status))
}

# ---- descriptions (no model) -------------------------------------------------------------------------------------------------

# Rate of an event with a person-clustered interval, one row.
wr_crude_cell <- function(y, id) {
  ci <- cluster_mean_ci(as.numeric(y), id)
  data.table(n_days = length(y), n_events = sum(y), n_people = ci$n_persons, rate_pct = 100 * ci$mean, rate_low = 100 * ci$lo, rate_high = 100 * ci$hi)
}

# Crude weekday rates by cohort (and all cohorts) for one outcome.
wr_crude <- function(rows, ycol) {
  by_co <- rows[, wr_crude_cell(get(ycol), person_id), by = .(cohort = as.character(cohort), weekday = wd)]
  all <- rows[, wr_crude_cell(get(ycol), person_id), by = .(weekday = wd)][, cohort := "All"]
  rbind(by_co, all, use.names = TRUE)[, outcome := if (ycol == "start") "gap start" else "gap end"][]
}

# Share of gaps of each length (1-7, 8+) by start weekday; each weekday adds to 100%.
wr_gap_length_table <- function(rows) {
  g <- rows[start == 1L, .(weekday = wd, len = factor(pmin(gap_length, 8L), levels = 1:8, labels = c(1:7, "8+")), person_id)]
  t <- g[, .(n_gaps = .N, n_people = uniqueN(person_id)), by = .(weekday, len)]
  t[, share_pct := 100 * n_gaps / sum(n_gaps), by = weekday]
  t[]
}

# Days at risk, gap starts and crude rate by weekday x days since return; ratio to the overall rate.
wr_dsr_cells <- function(rows, by_extra = character()) {
  overall <- rows[, 100 * mean(start)]
  t <- rows[, .(n_days = .N, n_events = sum(start), n_people = uniqueN(person_id)), by = c(by_extra, "wd", "dsr_f")]
  t[, `:=`(rate_pct = 100 * n_events / n_days, ratio_to_overall = (100 * n_events / n_days) / overall)]
  setnames(t, c("wd", "dsr_f"), c("weekday", "days_since_return"))
  t[]
}

# Crude gap-start rate by days since return for each device family with at least min_people people.
wr_dsr_by_device <- function(rows, cfg) {
  if (!"device" %in% names(rows)) return(data.table())
  keep <- rows[, .(p = uniqueN(person_id)), by = device][p >= cfg$weekday_rhythm$min_people, device]
  rows[device %in% keep, .(n_days = .N, n_events = sum(start), n_people = uniqueN(person_id), rate_pct = 100 * mean(start)), by = .(device = as.character(device), days_since_return = dsr_f)]
}

# Weekday of each person's first record (people counts).
wr_first_record_weekday <- function(rows) {
  rows[study_day == 1L, .(n_people = uniqueN(person_id)), by = .(cohort = as.character(cohort), weekday = wd)]
}

# Cells with 1 to min_people - 1 people are not shown: the count becomes NA with people_label "<20" and the values are NA; cells with fewer than
# min_cell_days days at risk are flagged `low_exposure` (kept in the table, drawn grey).
wr_suppress <- function(tab, value_cols, cfg, days_col = NULL) {
  if (!nrow(tab)) return(tab)
  tab <- copy(tab); lim <- cfg$weekday_rhythm$min_people
  small <- !is.na(tab$n_people) & tab$n_people > 0L & tab$n_people < lim
  set(tab, j = "people_label", value = ifelse(small, paste0("<", lim), as.character(tab$n_people)))
  for (v in intersect(value_cols, names(tab))) set(tab, which(small), v, NA)
  set(tab, which(small), "n_people", NA_integer_)
  set(tab, j = "suppressed", value = small)
  if (!is.null(days_col)) set(tab, j = "low_exposure", value = tab[[days_col]] < cfg$weekday_rhythm$min_cell_days)
  tab[]
}

# All descriptions of one stream (the value of the unit "DESC").
wr_describe <- function(start_rows, end_rows, stream, cfg) {
  vc <- c("n_days", "n_events", "rate_pct", "rate_low", "rate_high", "ratio_to_overall", "n_gaps", "share_pct")
  hol <- wr_holiday_dates(c(start_rows$out_date, end_rows$out_date))$table
  mk <- function(x) { if (!nrow(x)) return(x); x <- copy(x); x[, stream := stream]; setcolorder(x, "stream"); x[] }
  list(crude = mk(wr_suppress(rbind(wr_crude(start_rows, "start"), wr_crude(end_rows, "end")), vc, cfg, "n_days")),
       gap_length = mk(wr_suppress(wr_gap_length_table(start_rows), vc, cfg)),
       dsr_cells = mk(wr_suppress(wr_dsr_cells(start_rows), vc, cfg, "n_days")),
       dsr_by_device = mk(wr_suppress(wr_dsr_by_device(start_rows, cfg), vc, cfg, "n_days")),
       first_record_weekday = mk(wr_suppress(wr_first_record_weekday(start_rows), vc, cfg)),
       flow = mk(rbind(attr(start_rows, "flow"), attr(end_rows, "flow"), fill = TRUE)),
       holidays = hol)
}

# ---- results from the checkpoints ----------------------------------------------------------------------------------------------

# Note of a test row.
wr_test_note <- function(ok, same, id) {
  fcase(!ok, "a fit is missing or failed", !same, "random-effects structures differ: not tested",
        id == "H3", "variances on their boundary: the chi-square p value is conservative (too large)", default = "")
}

# Do both fits of a pair exist and have a log-likelihood?
wr_pair_ok <- function(big, small) nrow(big) == 1L && nrow(small) == 1L && !is.na(big$loglik) && !is.na(small$loglik)

# Likelihood-ratio statistic and its df (from the fits), NA when a fit is missing.
wr_lr_stats <- function(big, small) {
  if (!wr_pair_ok(big, small)) return(list(chisq = NA_real_, df = NA_integer_))
  list(chisq = max(0, 2 * (big$loglik - small$loglik)), df = big$df - small$df)
}

# One likelihood-ratio test of a stream: the big and the small fit of the pair p.
wr_test_row <- function(fits, s, p) {
  big <- fits[stream == s & fit_id == p$big]; small <- fits[stream == s & fit_id == p$small]
  ok <- wr_pair_ok(big, small); lr <- wr_lr_stats(big, small)
  same <- ok && (p$id == "H3" || identical(big$re_structure, small$re_structure))
  data.table(id = p$id, family = p$family, stream = s, outcome = p$outcome, fit_big = p$big, fit_small = p$small, chisq = lr$chisq, df = lr$df,
             p = if (same) stats::pchisq(lr$chisq, lr$df, lower.tail = FALSE) else NA_real_, re_structure = if (ok) big$re_structure else NA_character_,
             note = wr_test_note(ok, same, p$id))
}

# Likelihood-ratio tests of the specified pairs, Holm within each family (H2r has no family and is not adjusted). Both fits of a test must use the
# same random-effects structure (except H3, whose point is the difference); otherwise the test is NA with a note.
wr_tests <- function(fits) {
  pairs <- wr_test_pairs()
  rows <- rbindlist(lapply(unique(fits$stream), function(s) rbindlist(lapply(seq_len(nrow(pairs)), function(i) wr_test_row(fits, s, pairs[i])))))
  rows[, p_holm := stats::p.adjust(p, "holm"), by = family]
  rows[is.na(family), p_holm := NA_real_]
  rows[]
}

# Share of the weekday pattern explained by days since return: 1 - log(ptr of S3) / log(ptr of S1), per stream (point estimate).
wr_shrinkage <- function(summ) {
  s1 <- summ[fit_id == "S1", .(stream, ptr_S1 = peak_to_trough)]; s3 <- summ[fit_id == "S3", .(stream, ptr_S3 = peak_to_trough)]
  m <- merge(s1, s3, by = "stream")
  m[, shrinkage := 1 - log(ptr_S3) / log(ptr_S1)][]
}

# Rows, people and events of every fit, and the rows left out against all the rows of its type, with the reason (the filter of the fit).
wr_fit_flow <- function(fits, flow) {
  reasons <- c(none = "dropped by construction (more than days_so_far_max days so far in a gap) or a missing adjuster", holidays = "outcome day within the holiday window of a US federal holiday",
               first_gap = "after the person's first gap start", primary_rows = "no record on d - 1, or incomplete predictors of the primary model")
  spec <- wr_specs()[, .(fit_id, rows_type = fifelse(rows == "start", "gap start", "gap end"), reason = reasons[filter])]
  base <- flow[grepl("^(observed days|days in gaps)", step), .(stream, rows_type, base_rows = n_rows)]
  m <- merge(merge(fits[, .(stream, fit_id, n_rows, n_people, n_events)], spec, by = "fit_id"), base, by = c("stream", "rows_type"))
  m[, .(stream, rows_type, step = paste0("fit ", fit_id), n_rows, n_people, n_events, n_left_out = base_rows - n_rows, reason)]
}

# Collect every unit of the weekday stage into the output tables (names without the prefix weekday_rhythm_).
collect_weekday_results <- function(ckpt_dir, cfg, data_keys) {
  cu <- collect_units(ckpt_dir, "weekday", cfg, data_keys)
  is_desc <- cu$units$id == "DESC"
  fitv <- cu$values[!is_desc]; fit_streams <- cu$units$stream[!is_desc]; fit_ids <- cu$units$id[!is_desc]
  desc <- cu$values[is_desc]
  stack <- function(name) rbindlist(Map(function(v, st, id) if (!is.null(v[[name]])) cbind(stream = st, fit_id = id, v[[name]]), fitv, fit_streams, fit_ids), fill = TRUE)
  fits <- rbindlist(lapply(fitv, `[[`, "info"), fill = TRUE)
  d <- function(name) rbindlist(lapply(desc, `[[`, name), fill = TRUE)
  out <- list(fits = fits, coefficients = stack("coef"), profile = stack("profile"), summary = stack("summary"), dsr_curve = stack("curve"), personal = stack("personal"),
              crude = d("crude"), gap_length = d("gap_length"), dsr_cells = d("dsr_cells"), dsr_by_device = d("dsr_by_device"),
              first_record_weekday = d("first_record_weekday"), flow = d("flow"), holidays = unique(d("holidays")[, .(holiday, date, observed)]),
              blup_bins = stack("blup_bins"), ellipse = stack("ellipse"), failures = cu$failures)
  if (nrow(fits)) { out$tests <- wr_tests(fits); out$flow <- rbind(out$flow, wr_fit_flow(fits, out$flow), fill = TRUE); if (nrow(out$summary)) out$shrinkage <- wr_shrinkage(out$summary) }
  out
}


# ---- output names, figures list and the plain-language summary ---------------------------------------------------------------------

# The tables of the stage under their output names (weekday_rhythm_<name>); the binned BLUPs and the ellipse only feed Figure 3.
wr_name_tables <- function(wr) {
  keep <- setdiff(names(wr), c("blup_bins", "ellipse"))
  stats::setNames(wr[keep], paste0("weekday_rhythm_", keep))
}

# One line per stream and outcome with the main effect sizes (first), then the tests.
write_weekday_summary <- function(wr, cfg, outdir) {
  nm <- c(activity = "activity", sleep = "sleep", hr = "heart rate"); fmt <- function(x, d = 2) formatC(x, digits = d, format = "f")
  L <- c("# Weekly rhythm of gap starts and gap ends", "",
         "Effect sizes first. With tens of thousands of gaps every test is likely to be significant: the p values mostly show that a pattern is not exactly zero.",
         "Ratios are ratios of predicted probabilities (population average), not odds ratios. Associations, not causes.", "", "## Shared weekly rhythm (weekday only; S1 gap starts, E1 gap ends)")
  s <- wr$summary[fit_id %in% c("S1", "E1")]
  for (i in seq_len(nrow(s))) L <- c(L, sprintf("- %s, gap %s: the busiest weekday is %s and the quietest %s; peak-to-trough ratio %s (95%% interval %s-%s); average %s%%.",
                                                nm[[s$stream[i]]], if (s$fit_id[i] == "S1") "start" else "end", WR_WEEKDAYS[s$peak_weekday[i] + 1L], WR_WEEKDAYS[s$trough_weekday[i] + 1L],
                                                fmt(s$peak_to_trough[i]), fmt(s$ptr_low[i]), fmt(s$ptr_high[i]), fmt(s$average_pct[i])))
  L <- c(L, "", "## Calendar routine or charging cycle (gap starts)")
  for (i in seq_len(nrow(wr$shrinkage))) L <- c(L, sprintf("- %s: days since return explain %s%% of the weekday pattern (peak-to-trough %s without, %s with days since return).", nm[[wr$shrinkage$stream[i]]],
                                                           fmt(100 * wr$shrinkage$shrinkage[i], 0), fmt(wr$shrinkage$ptr_S1[i]), fmt(wr$shrinkage$ptr_S3[i])))
  L <- c(L, "- Reading rule: if the weekday pattern holds after days since return it is a calendar routine; if it shrinks to about nothing it was a cycle; if both matter, both are present.")
  if (nz(wr$personal)) {
    L <- c(L, "", "## Personal rhythm (gap starts)")
    for (i in seq_len(nrow(wr$personal))) L <- c(L, sprintf("- %s: person peak-to-trough ratio, median %s (10th-90th percentile %s-%s); %s%% of people peak within one day of the shared peak (model based, not shrunk).",
                                                            nm[[wr$personal$stream[i]]], fmt(wr$personal$person_ptr_median[i]), fmt(wr$personal$person_ptr_p10[i]), fmt(wr$personal$person_ptr_p90[i]),
                                                            fmt(100 * wr$personal$share_near_shared_peak[i], 0)))
  }
  t <- wr$tests[!is.na(p)]
  L <- c(L, "", "## Tests (likelihood-ratio; Holm within the family H1-H4)", sprintf("- %s %s: chisq %s on %s df, Holm p %s.", t$id, paste0("(", sub(" ", "", t$outcome), ", ", t$stream, ")"), fmt(t$chisq, 1), t$df, wr_ptext(t$p_holm)))
  bad <- wr$fits[status != "ok"]
  if (nrow(bad)) L <- c(L, "", "## Fits that did not run cleanly", sprintf("- %s %s: %s %s", bad$stream, bad$fit_id, bad$status, bad$note))
  writeLines(L, file.path(outdir, "results_summary.md"))
  invisible(L)
}
