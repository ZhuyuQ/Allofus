# ==============================================================================
# temporal.R
# ------------------------------------------------------------------------------
# Temporal structure of missingness at the DAY level:
#   1. onset vs continuation (do missing days come in runs?)
#   2. missingness by study week and by calendar year
#   3. day-of-week clustering
#
# All counts come from the gap table and arithmetic on each person's window; the
# full calendar grid is never built. Uncertainty is a person-level (cluster)
# bootstrap, because days of one person are not independent.
#
# Note on edge effects: the first and last day of every window are observed by
# construction (window (a)), so missingness is slightly under-estimated in the
# very first and last study week.
# ==============================================================================

# ---- 1. onset vs continuation ----------------------------------------------------------

# Per person counts needed for the transition probabilities.
#   at_risk   : observed days that have a next day              (n_observed - 1)
#   onsets    : of those, days followed by a missing day          (= n_gaps)
#   miss_days : missing days                                       (= n_missing)
#   continue  : missing days followed by another missing day       (= n_missing - n_gaps)
transition_person_counts <- function(pg) {
  pg[, .(person_id, cohort, stream, at_risk = n_observed - 1L, onsets = n_gaps,
         miss_days = n_missing, continue_days = n_missing - n_gaps, window_days)]
}

# P(missing tomorrow | observed today)   = onset probability
# P(missing tomorrow | missing today)    = continuation probability
# persistence ratio = continuation / overall missing share. Under "missing days are
# scattered at random" it is 1; above 1 missing days come in runs.
transition_table <- function(counts, cfg) {
  one <- function(sub) {
    stat <- function(s) {
      p_on <- sum(s$onsets) / sum(s$at_risk)
      p_co <- if (sum(s$miss_days) > 0) sum(s$continue_days) / sum(s$miss_days) else NA_real_
      p_mi <- sum(s$miss_days) / sum(s$window_days)
      c(p_onset = p_on, p_continue = p_co, p_missing = p_mi, persistence_ratio = p_co / p_mi)
    }
    est <- stat(sub); n <- nrow(sub)
    draws <- if (n >= 2L) vapply(seq_len(cfg$boot$B), function(b) stat(sub[sample.int(n, n, TRUE)]), numeric(4))
             else matrix(NA_real_, 4, 1)
    ci <- apply(draws, 1L, stats::quantile, probs = c(0.025, 0.975), na.rm = TRUE, names = FALSE)
    out <- data.table(n_persons = n)
    for (i in seq_along(est)) {
      nm <- names(est)[i]
      out[, (nm) := est[i]]; out[, (paste0(nm, "_low")) := ci[1L, i]]; out[, (paste0(nm, "_high")) := ci[2L, i]]
    }
    out
  }
  with_seed(cfg$seed, counts[, one(.SD), by = .(stream, cohort)])
}

# ---- 2. study week / calendar year ----------------------------------------------------------

# Missing days by study week (week 1 = first 7 days of the person's window).
# Returns person x week counts: num = missing days, den = window days in that week.
study_week_counts <- function(pg, gaps) {
  wk <- pg[, .(week = seq_len(ceiling(window_days / 7)), window_days, first_date, cohort),
           by = person_id]
  wk[, den := pmin(window_days, 7L * week) - 7L * (week - 1L)]
  gd <- expand_gap_days(gaps)
  gd[pg, on = "person_id", first_date := i.first_date]
  gd[, week := (as.integer(date - first_date)) %/% 7L + 1L]
  miss <- gd[, .(num = .N), by = .(person_id, week)]
  wk <- merge(wk, miss, by = c("person_id", "week"), all.x = TRUE)
  wk[is.na(num), num := 0L]
  wk[, .(person_id, cohort, week, num, den)]
}

# Missing days by calendar year: window days per person-year come from the overlap of
# the window with each calendar year.
calendar_year_counts <- function(pg, gaps) {
  yrs <- seq(year(min(pg$first_date)), year(max(pg$last_date)))
  win <- rbindlist(lapply(yrs, function(y) {
    lo <- as.Date(sprintf("%d-01-01", y)); hi <- as.Date(sprintf("%d-12-31", y))
    pg[, .(person_id, cohort, year = y,
           den = pmax(0L, as.integer(pmin(last_date, hi) - pmax(first_date, lo)) + 1L))]
  }))
  win <- win[den > 0L]
  gd <- expand_gap_days(gaps)[, .(num = .N), by = .(person_id, year = year(date))]
  out <- merge(win, gd, by = c("person_id", "year"), all.x = TRUE)
  out[is.na(num), num := 0L]
  out[]
}

# Pooled missing share per group with a cluster-bootstrap interval; only groups with
# at least `min_persons` contributing people are returned.
rate_table <- function(counts, group_cols, cfg, min_persons = 30L) {
  g <- counts[, .(n = uniqueN(person_id)), by = group_cols][n >= min_persons]
  d <- merge(counts, g[, ..group_cols], by = group_cols)
  cluster_boot_ratio(d, group_cols, num = "num", den = "den", B = cfg$boot$B, seed = cfg$seed)
}

# ---- 3. day of week ------------------------------------------------------------------------------------

DOW_LABELS <- c("Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat")

# Number of each weekday inside a window [first_date, first_date + window_days - 1].
# Day i of the window has weekday (f + i) mod 7; weekday w occurs at i = r, r+7, ...
# with r = (w - f) mod 7, so the count is floor((W - 1 - r) / 7) + 1 (never below 0).
window_weekday_counts <- function(first_date, window_days) {
  f <- wday(first_date) - 1L                      # 0 = Sunday
  rbindlist(lapply(0:6, function(w) {
    r <- (w - f) %% 7L
    data.table(row = seq_along(first_date),
               dow = factor(DOW_LABELS[w + 1L], levels = DOW_LABELS),
               den = pmax(0L, (window_days - 1L - r) %/% 7L + 1L))
  }))
}

# Person x weekday counts: den = days of that weekday in the window, num = missing ones.
dow_counts <- function(pg, gaps) {
  wc <- window_weekday_counts(pg$first_date, pg$window_days)
  base <- wc[, .(person_id = pg$person_id[row], cohort = pg$cohort[row], dow, den)]
  gd <- expand_gap_days(gaps)
  gd[, dow := factor(DOW_LABELS[wday(date)], levels = DOW_LABELS)]
  miss <- gd[, .(num = .N), by = .(person_id, dow)]
  out <- merge(base, miss, by = c("person_id", "dow"), all.x = TRUE)
  out[is.na(num), num := 0L]
  out[]
}

# Does the missing share depend on the weekday, and does that differ by cohort?
#   m0: cohort            m1: + weekday            m2: + weekday x cohort
# Likelihood-ratio tests: m1 vs m0 (weekday effect), m2 vs m1 (differs by cohort).
# A person random intercept allows for people who are simply more gap-prone.
fit_dow_model <- function(counts, persons, cfg, adjusted = FALSE) {
  d <- merge(counts[, .(person_id, dow, num, den)],
             persons[, .(person_id, cohort, age10 = (age_first - 60) / 10,
                         sex2 = factor(as.character(sex), levels = c("Female", "Male")))],
             by = "person_id")
  d[, obs := den - num]
  adj <- if (adjusted) c("age10", "sex2") else character()
  if (adjusted) d <- d[complete.cases(d[, ..adj])]
  fit <- function(terms) {
    f <- stats::as.formula(paste("cbind(num, obs) ~", paste(c(terms, adj, "(1 | person_id)"), collapse = " + ")))
    lme4::glmer(f, data = d, family = stats::binomial(), nAGQ = cfg$model$nAGQ,
                control = lme4::glmerControl(optimizer = "bobyqa"))
  }
  m0 <- fit("cohort"); m1 <- fit(c("dow", "cohort")); m2 <- fit(c("dow * cohort"))
  t1 <- lrt(m0, m1); t2 <- lrt(m1, m2)
  tests <- data.table(test = c("weekday effect", "weekday x cohort"), chisq = c(t1$chisq, t2$chisq),
                      df = c(t1$df, t2$df), p = c(t1$p, t2$p), adjusted = adjusted,
                      n_persons = uniqueN(d$person_id))
  cf <- summary(m1)$coefficients
  cf <- cf[grepl("^dow", rownames(cf)), , drop = FALSE]
  coefs <- data.table(term = sub("^dow", "", rownames(cf)), estimate = cf[, 1], se = cf[, 2], p = cf[, 4],
                      or_vs_sunday = exp(cf[, 1]), or_low = exp(cf[, 1] - 1.96 * cf[, 2]),
                      or_high = exp(cf[, 1] + 1.96 * cf[, 2]), adjusted = adjusted)
  list(tests = tests, coefs = coefs)
}
