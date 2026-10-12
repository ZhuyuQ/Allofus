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

# ---- 4. attained age (continuous) --------------------------------------------------------------------------------

# The persons' exact birth dates are not kept, only the age in whole years at the first record
# (`age_first`). The attained age at day e (0 = first record) is therefore taken as
#     age_first + 0.5 + e / 365.25 ,
# i.e. the middle of the unknown birthday year; each person's age is uncertain by about +/- 0.5 year, which smooths
# the curve slightly but does not bias it.
# Missing days and window days are then counted by whole year of attained age. A day belongs to year `a` when
# a <= attained age < a + 1; the same thresholds are used for both counts, so the window days of a person add up
# exactly to window_days and the missing days to n_missing.
age_year_of_day <- function(age0, e) {
  a <- floor(age0 + e / 365.25)
  a <- a + (ceiling((a + 1 - age0) * 365.25) <= e) - (ceiling((a - age0) * 365.25) > e)   # exact at the boundaries
  a
}

# Person x attained-age counts: den = window days spent at that age, num = missing ones.
# half_width = 0: one row per whole year of age. half_width = h > 0: a MOVING WINDOW - the row for age a pools the days
# at ages a - h ... a + h (h = 2: a 5-year window), which keeps enough people per age to be reported and gives a smooth
# curve over a continuous age axis. Each person's days then appear in 2h + 1 neighbouring ages (the bootstrap still
# resamples whole people, so the bands stay valid).
age_counts <- function(pg, gaps, persons, half_width = 0L) {
  p <- merge(pg[, .(person_id, cohort, first_date, window_days)], persons[, .(person_id, age_first)], by = "person_id")
  p <- p[!is.na(age_first)][, age0 := age_first + 0.5]
  lo <- age_year_of_day(p$age0, 0L); hi <- age_year_of_day(p$age0, p$window_days - 1L)
  grid <- p[rep(seq_len(nrow(p)), hi - lo + 1L), .(person_id, cohort, age0, window_days)]
  grid[, age := rep(lo, hi - lo + 1L) + sequence(hi - lo + 1L) - 1L]
  e_lo <- pmax(0, ceiling((grid$age - grid$age0) * 365.25))
  e_hi <- pmin(grid$window_days - 1, ceiling((grid$age + 1 - grid$age0) * 365.25) - 1)
  grid[, den := pmax(0, e_hi - e_lo + 1)]

  gd <- expand_gap_days(gaps)
  gd[p, on = "person_id", `:=`(age0 = i.age0, first_date = i.first_date)]
  gd <- gd[!is.na(age0)]
  gd[, age := age_year_of_day(age0, as.integer(date - first_date))]
  miss <- gd[, .(num = .N), by = .(person_id, age)]
  out <- merge(grid[den > 0, .(person_id, cohort, age, den = as.integer(den))], miss, by = c("person_id", "age"), all.x = TRUE)
  out[is.na(num), num := 0L]
  out <- out[, .(person_id, cohort, age, num, den)]
  if (half_width > 0L) {
    out <- rbindlist(lapply(-half_width:half_width, function(off) copy(out)[, age := age - off]))
    out <- out[, .(num = sum(num), den = sum(den)), by = .(person_id, cohort, age)]
  }
  out[]
}

# ---- the person x unit counts and the statistics behind the "missing by ..." figures --------------------------------------

# Person x (study week | calendar year | weekday | age) counts of missing days (num) and window days (den), all streams.
# `age` has one row per whole year of age (for the tests), `age_smooth` the 5-year moving window (for the figure).
missing_rate_counts <- function(pg, gaps, persons, streams, cfg) {
  per_stream <- function(f) rbindlist(lapply(streams, function(s) cbind(stream = s, f(pg[stream == s], gaps[[s]]))))
  list(week = per_stream(study_week_counts), year = per_stream(calendar_year_counts), dow = per_stream(dow_counts),
       age = per_stream(function(p, g) age_counts(p, g, persons, 0L)),
       age_smooth = per_stream(function(p, g) age_counts(p, g, persons, cfg$strata$age_half_width)))
}

# Wald chi-square test of the linear hypothesis R b = 0 with covariance V.
wald_chisq <- function(b, V, R) {
  Rb <- R %*% b
  stat <- tryCatch(drop(t(Rb) %*% solve(R %*% V %*% t(R)) %*% Rb), error = function(e) NA_real_)
  list(chisq = stat, df = nrow(R), p = if (is.na(stat)) NA_real_ else stats::pchisq(stat, nrow(R), lower.tail = FALSE))
}

# Pooled logistic regression of the aggregated counts, cbind(missing, observed) ~ formula, with PERSON-CLUSTERED sandwich
# covariance (sandwich::vcovCL). The estimates are those of the binomial model; the clustering makes the standard errors valid
# for days that are correlated within a person (and for extra-binomial variation). Returns NULL when the model cannot be fitted.
cluster_logit <- function(d, rhs) {
  d <- copy(d)[, obs := den - num]
  fit <- tryCatch(suppressWarnings(stats::glm(stats::as.formula(paste("cbind(num, obs) ~", rhs)), data = d, family = stats::binomial())),
                  error = function(e) NULL)
  if (is.null(fit) || !isTRUE(fit$converged)) return(NULL)
  V <- tryCatch(sandwich::vcovCL(fit, cluster = d$person_id), error = function(e) NULL)
  if (is.null(V)) return(NULL)
  b <- stats::coef(fit); ok <- !is.na(b)
  list(b = b[ok], V = V[ok, ok, drop = FALSE])
}

TREND_VARIABLES <- list(study_week   = list(table = "week", col = "week", unit = 4,  label = "per 4 study weeks"),
                        calendar_year = list(table = "year", col = "year", unit = 1,  label = "per calendar year"),
                        age           = list(table = "age",  col = "age",  unit = 10, label = "per 10 years of age"))

# Is the missing share rising or falling with study week, calendar year or age, and does that differ between cohorts?
# For every stream and variable x (centred, in units `unit`):
#   * per cohort: logit(P(day missing)) = a + b x  ->  odds ratio exp(b) per unit, 95% CI, z test (robust SE);
#     Holm-adjusted over the cohorts of one stream x variable;
#   * all cohorts: the same with cohort as a covariate (a common slope);
#   * heterogeneity: Wald chi-square test that the slopes differ between cohorts (cohort x x interaction, k - 1 df).
# A LINEAR trend on the logit scale is a summary; the figure shows the actual shape (the model does not test curvature).
missing_trend_statistics <- function(counts, min_persons = 30L) {
  trend <- list(); het <- list()
  for (v in names(TREND_VARIABLES)) {
    tv <- TREND_VARIABLES[[v]]
    for (s in unique(counts[[tv$table]]$stream)) {
      d <- counts[[tv$table]][stream == s]; d[, xs := (get(tv$col) - mean(get(tv$col))) / tv$unit]
      d[, cohort := droplevels(factor(cohort))]
      for (co in levels(d$cohort)) {
        dc <- d[cohort == co]
        row <- data.table(stream = s, variable = v, unit = tv$label, cohort = co, n_persons = uniqueN(dc$person_id), n_rows = nrow(dc))
        m <- if (row$n_persons >= min_persons && uniqueN(dc$xs) >= 2L) cluster_logit(dc, "xs") else NULL
        trend[[length(trend) + 1L]] <- if (is.null(m) || !"xs" %in% names(m$b)) cbind(row, status = "not_estimable") else {
          est <- m$b[["xs"]]; se <- sqrt(m$V["xs", "xs"])
          cbind(row, estimate = est, se_robust = se, z = est / se, p = 2 * stats::pnorm(-abs(est / se)),
                odds_ratio = exp(est), or_low = exp(est - 1.96 * se), or_high = exp(est + 1.96 * se), status = "ok")
        }
      }
      m_all <- if (nlevels(d$cohort) >= 2L) cluster_logit(d, "cohort + xs") else cluster_logit(d, "xs")
      if (!is.null(m_all) && "xs" %in% names(m_all$b)) {
        est <- m_all$b[["xs"]]; se <- sqrt(m_all$V["xs", "xs"])
        trend[[length(trend) + 1L]] <- data.table(stream = s, variable = v, unit = tv$label, cohort = "All cohorts (adjusted for cohort)",
          n_persons = uniqueN(d$person_id), n_rows = nrow(d), estimate = est, se_robust = se, z = est / se,
          p = 2 * stats::pnorm(-abs(est / se)), odds_ratio = exp(est), or_low = exp(est - 1.96 * se), or_high = exp(est + 1.96 * se), status = "ok")
      }
      m_int <- if (nlevels(d$cohort) >= 2L) cluster_logit(d, "cohort * xs") else NULL
      het[[length(het) + 1L]] <- if (is.null(m_int)) data.table(stream = s, variable = v, test = "slope differs between cohorts",
                                                                chisq = NA_real_, df = NA_integer_, p = NA_real_, status = "not_estimable") else {
        idx <- grep(":xs$|^xs:", names(m_int$b)); R <- diag(length(m_int$b))[idx, , drop = FALSE]
        w <- wald_chisq(m_int$b, m_int$V, R)
        data.table(stream = s, variable = v, test = "slope differs between cohorts", chisq = w$chisq, df = w$df, p = w$p, status = "ok")
      }
    }
  }
  trend <- rbindlist(trend, fill = TRUE); het <- rbindlist(het, fill = TRUE)
  if (nrow(trend) > 0L) trend[!grepl("^All cohorts", cohort), p_holm := stats::p.adjust(p, "holm"), by = .(stream, variable)]
  if (nrow(het) > 0L) het[, p_holm := stats::p.adjust(p, "holm"), by = variable]
  list(missing_trend_tests = trend, missing_trend_heterogeneity = het)
}

# Day-of-week clustering: is the missing share different on some weekdays?
#   * per cohort: Wald chi-square test of the 6 weekday terms (robust SE), Holm over the cohorts of one stream;
#   * all cohorts (adjusted for cohort): the same with cohort as a covariate;
#   * weekday x cohort interaction: do the weekday patterns differ between cohorts (6 x (k - 1) df).
# Effect sizes: odds ratio of a day being missing for each weekday versus Sunday, per cohort, with 95% CI.
# (Complements the mixed model in weekday_tests / weekday_odds_ratios, which needs the temporal stage.)
missing_weekday_statistics <- function(counts, min_persons = 30L) {
  tests <- list(); ors <- list()
  for (s in unique(counts$dow$stream)) {
    d <- counts$dow[stream == s]; d[, `:=`(dow = factor(dow, levels = DOW_LABELS), cohort = droplevels(factor(cohort)))]
    one <- function(dd, label, rhs) {
      m <- if (uniqueN(dd$person_id) >= min_persons) cluster_logit(dd, rhs) else NULL
      if (is.null(m)) return(NULL)
      idx <- grep("^dow", names(m$b)); if (!length(idx)) return(NULL)
      R <- diag(length(m$b))[idx, , drop = FALSE]; w <- wald_chisq(m$b, m$V, R)
      list(test = data.table(stream = s, cohort = label, test = "weekday effect", n_persons = uniqueN(dd$person_id),
                             chisq = w$chisq, df = w$df, p = w$p),
           or = data.table(stream = s, cohort = label, weekday = sub("^dow", "", names(m$b)[idx]), estimate = m$b[idx],
                           se_robust = sqrt(diag(m$V)[idx]))[, `:=`(odds_ratio_vs_sunday = exp(estimate), or_low = exp(estimate - 1.96 * se_robust),
                                                                    or_high = exp(estimate + 1.96 * se_robust), p = 2 * stats::pnorm(-abs(estimate / se_robust)))])
    }
    for (co in levels(d$cohort)) { r <- one(d[cohort == co], co, "dow"); if (!is.null(r)) { tests[[length(tests) + 1L]] <- r$test; ors[[length(ors) + 1L]] <- r$or } }
    r <- one(d, "All cohorts (adjusted for cohort)", if (nlevels(d$cohort) >= 2L) "cohort + dow" else "dow")
    if (!is.null(r)) { tests[[length(tests) + 1L]] <- r$test; ors[[length(ors) + 1L]] <- r$or }
    if (nlevels(d$cohort) >= 2L) {
      m <- cluster_logit(d, "cohort * dow")
      if (!is.null(m)) {
        idx <- grep("^cohort.*:dow|^dow.*:cohort", names(m$b)); w <- wald_chisq(m$b, m$V, diag(length(m$b))[idx, , drop = FALSE])
        tests[[length(tests) + 1L]] <- data.table(stream = s, cohort = "weekday pattern differs between cohorts", test = "weekday x cohort",
                                                  n_persons = uniqueN(d$person_id), chisq = w$chisq, df = w$df, p = w$p)
      }
    }
  }
  tests <- rbindlist(tests, fill = TRUE); ors <- rbindlist(ors, fill = TRUE)
  if (nrow(tests) > 0L) tests[test == "weekday effect" & !grepl("^All cohorts", cohort), p_holm := stats::p.adjust(p, "holm"), by = stream]
  list(missing_weekday_tests = tests, missing_weekday_odds_ratios = ors)
}

# Everything behind the four "missing by ..." figures: the pooled missing share per cohort with a cluster-bootstrap
# interval (the lines and bands) and the statistics tables (trend, weekday, heterogeneity).
missing_rate_statistics <- function(counts, cfg) {
  c(missing_trend_statistics(counts), missing_weekday_statistics(counts))
}

# Pooled missing share per cohort with a cluster-bootstrap interval, by study week, calendar year,
# weekday and attained age, for every stream, plus the statistics tables. Only counts are used, so this is cheap.
missing_rate_tables <- function(pg, gaps, persons, streams, cfg) {
  cn <- missing_rate_counts(pg, gaps, persons, streams, cfg)
  c(list(study_week_rates    = rate_table(cn$week, c("stream", "cohort", "week"), cfg),
         calendar_year_rates = rate_table(cn$year, c("stream", "cohort", "year"), cfg),
         weekday_rates       = rate_table(cn$dow, c("stream", "cohort", "dow"), cfg, min_persons = 1L),
         age_rates           = rate_table(cn$age_smooth, c("stream", "cohort", "age"), cfg)),
    missing_rate_statistics(cn, cfg))
}
