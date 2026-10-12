# ==============================================================================
# summarise.R
# ------------------------------------------------------------------------------
# Person-level missingness -> cohort tables and cohort comparisons.
# The unit of comparison is the PERSON (one row each, from summarise_person_gaps()),
# so people with long follow-up do not dominate and tests are not pseudo-replicated.
# ==============================================================================

# Table-2-style summary per cohort: "median (Q1-Q3)" of follow-up, record days,
# missing % and longest gap, plus the share of people with any gap / a gap >= 7 days.
cohort_missingness_table <- function(pg) {
  pg[, .(
    n_persons          = .N,
    followup_days      = format_med_iqr(window_days),
    record_days        = format_med_iqr(n_observed),
    missing_pct        = format_med_iqr(pct_missing),
    max_gap_days       = format_med_iqr(gap_max),
    pct_any_gap        = round(100 * mean(n_gaps > 0), 1),
    pct_gap_ge7        = round(100 * mean(gap_max >= 7), 1),
    gaps_per_person    = format_med_iqr(n_gaps)
  ), by = .(stream, cohort)][order(stream, cohort)]
}

# Cliff's delta (effect size, -1..1) from the Mann-Whitney statistic:
# P(x > y) - P(x < y). Positive: x tends to be larger.
cliffs_delta <- function(x, y) {
  x <- x[!is.na(x)]; y <- y[!is.na(y)]
  if (length(x) == 0L || length(y) == 0L) return(NA_real_)
  w <- suppressWarnings(stats::wilcox.test(x, y, exact = FALSE))$statistic
  unname(2 * w / (length(x) * length(y)) - 1)
}

# Kruskal-Wallis across cohorts + all pairwise Wilcoxon tests (Holm-adjusted) with
# Cliff's delta, for one person-level variable within one stream.
compare_cohorts <- function(pg, var, stream_name = NULL) {
  d <- if (is.null(stream_name)) pg else pg[stream == stream_name]
  d <- d[!is.na(get(var))]
  kw <- stats::kruskal.test(d[[var]] ~ d$cohort)
  overall <- data.table(stream = stream_name %||% "all", variable = var,
                        chisq = unname(kw$statistic), df = unname(kw$parameter), p = kw$p.value)
  lv <- levels(droplevels(d$cohort))
  pairs <- if (length(lv) >= 2L) t(utils::combn(lv, 2L)) else matrix(character(), 0, 2)
  pw <- rbindlist(lapply(seq_len(nrow(pairs)), function(i) {
    x <- d[cohort == pairs[i, 1L]][[var]]; y <- d[cohort == pairs[i, 2L]][[var]]
    data.table(stream = stream_name %||% "all", variable = var, group1 = pairs[i, 1L],
               group2 = pairs[i, 2L], n1 = length(x), n2 = length(y),
               median1 = stats::median(x), median2 = stats::median(y),
               p_raw = suppressWarnings(stats::wilcox.test(x, y, exact = FALSE)$p.value),
               cliffs_delta = cliffs_delta(x, y))
  }))
  if (nrow(pw) > 0L) pw[, p_holm := stats::p.adjust(p_raw, "holm")]
  list(overall = overall, pairwise = pw)
}

# Compare the share of missing days between cohorts with a quasi-binomial model:
#   cbind(missing days, observed days) ~ cohort [+ age + sex]
# so a person's weight reflects how many days they were followed. Odds ratios are
# for the odds that a day inside the window is missing, versus `ref`.
fit_missing_fraction <- function(pg, persons, adjusted = FALSE, ref = "Healthy") {
  d <- merge(pg[, .(person_id, n_missing, n_observed)],
             persons[, .(person_id, cohort, age10 = (age_first - 60) / 10,
                         sex2 = factor(as.character(sex), levels = c("Female", "Male")))],
             by = "person_id")
  d[, cohort := stats::relevel(factor(cohort), ref = ref)]
  rhs <- c("cohort", if (adjusted) c("age10", "sex2"))
  if (adjusted) d <- d[complete.cases(d)]
  if (adjusted && uniqueN(d$sex2) < 2L) rhs <- setdiff(rhs, "sex2")
  fit <- stats::glm(stats::as.formula(paste("cbind(n_missing, n_observed) ~", paste(rhs, collapse = " + "))),
                    data = d, family = stats::quasibinomial())
  cf <- summary(fit)$coefficients
  out <- data.table(term = rownames(cf), estimate = cf[, 1], se = cf[, 2], p = cf[, 4])
  out[, `:=`(or = exp(estimate), or_low = exp(estimate - 1.96 * se), or_high = exp(estimate + 1.96 * se),
             adjusted = adjusted, n_persons = nrow(d), dispersion = summary(fit)$dispersion)]
  out[]
}

# ---- person-level % missing by race / sex / CCI group / age group, within each cohort ---------------------

# Person-level % missing with the strata attached (one row per person x stream).
#   pg      person-level gap summary with cohort and stream (see summarise_person_gaps)
#   persons the person table (sex, age_group, race_plot, cci_group)
missing_pct_long <- function(pg, persons) {
  strata <- persons[, .(person_id, sex, age_group, race_plot, cci_group)]
  out <- merge(pg[, .(person_id, cohort, stream, metric = stream, value = pct_missing)], strata, by = "person_id")
  out[, cohort := factor(cohort, levels = levels(persons$cohort))]
  out[]
}

# Epsilon-squared effect size of a Kruskal-Wallis test: H / (n - 1), range 0..1 (Tomczak & Tomczak 2014).
kw_epsilon_sq <- function(H, n) unname(H / (n - 1))

# One test for the values of `value` across the groups `g` (groups with fewer than `min_n` people are left out):
#   >= 3 groups: Kruskal-Wallis (statistic = H, df = groups - 1, effect = epsilon-squared)
#   2 groups   : Wilcoxon rank-sum, normal approximation (statistic = W, effect = Cliff's delta, first group vs second)
#   < 2 groups : no test
rank_test_groups <- function(value, g, min_n = 5L) {
  g <- droplevels(factor(g)); ok <- !is.na(value) & !is.na(g)
  value <- value[ok]; g <- droplevels(g[ok])
  keep <- names(which(table(g) >= min_n)); g <- factor(g, levels = keep)
  use <- !is.na(g); value <- value[use]; g <- droplevels(g[use])
  k <- nlevels(g); n <- length(value)
  none <- list(test = "none", n_total = n, n_groups = k, statistic = NA_real_, df = NA_real_, p = NA_real_,
               effect = NA_real_, effect_name = NA_character_)
  if (k < 2L) return(none)
  if (k == 2L) {
    x <- value[g == levels(g)[1L]]; y <- value[g == levels(g)[2L]]
    w <- suppressWarnings(stats::wilcox.test(x, y, exact = FALSE))
    return(list(test = "Wilcoxon rank-sum", n_total = n, n_groups = 2L, statistic = unname(w$statistic), df = NA_real_,
                p = w$p.value, effect = cliffs_delta(x, y), effect_name = "Cliff's delta"))
  }
  kw <- stats::kruskal.test(value ~ g)
  list(test = "Kruskal-Wallis", n_total = n, n_groups = k, statistic = unname(kw$statistic), df = unname(kw$parameter),
       p = kw$p.value, effect = kw_epsilon_sq(kw$statistic, n), effect_name = "epsilon-squared")
}

# Tests for the person-level % missing by stratum (`stratum` = sex | age_group | race_plot | cci_group):
#   family "groups_within_cohort": do the levels of the stratum differ in % missing, within each cohort?
#                                  (the question of the figure; Holm-adjusted over the cohorts of one stream x stratum)
#   family "cohorts_within_level": do the cohorts differ within each level of the stratum? (the comparison of
#                                  the earlier cohort figures; Holm-adjusted over the levels)
# Every test is a rank test on one value per person, so people with long follow-up do not dominate.
missing_pct_stratum_tests <- function(long, stratum, min_n = 5L) {
  rows <- list()
  for (s in unique(long$stream)) {
    d <- long[stream == s & !is.na(get(stratum))]
    d[, lev := factor(get(stratum))]
    for (co in levels(d$cohort)) {
      sub <- d[cohort == co]
      if (nrow(sub) == 0L) next
      rows[[length(rows) + 1L]] <- c(list(stream = s, stratum_var = stratum, family = "groups_within_cohort", group = co),
                                     rank_test_groups(sub$value, sub$lev, min_n))
    }
    for (lv in levels(d$lev)) {
      sub <- d[lev == lv]
      rows[[length(rows) + 1L]] <- c(list(stream = s, stratum_var = stratum, family = "cohorts_within_level", group = lv),
                                     rank_test_groups(sub$value, droplevels(sub$cohort), min_n))
    }
  }
  out <- rbindlist(rows, fill = TRUE)
  if (nrow(out) > 0L) out[, p_holm := stats::p.adjust(p, "holm"), by = .(stream, stratum_var, family)]
  out[]
}

# Odds ratios for the share of missing days by stratum level, WITHIN each cohort and stream:
# quasi-binomial GLM  cbind(missing days, observed days) ~ level  (reference = the first level with >= min_n people), one row per
# person so each person is weighted by the days followed. Reported: odds ratio of a day being missing, level vs reference, with
# 95% CI and the t-test p value of the quasi-binomial model, the dispersion, and the omnibus F test of the stratum
# (all levels at once) with its Holm-adjusted p value over the cohorts of one stream x stratum.
# These are unadjusted for the other strata (age, sex, race and CCI are correlated); the cohort comparison adjusted for age and sex is in
# missing_fraction_models.
missing_pct_stratum_models <- function(pg, persons, min_n = 5L) {
  strata <- c("sex", "age_group", "race_plot", "cci_group")
  d0 <- merge(pg[, .(person_id, stream, cohort, n_missing, n_observed)], persons[, c("person_id", strata), with = FALSE], by = "person_id")
  d0[, cohort := factor(cohort, levels = levels(persons$cohort))]
  rows <- list()
  for (st in strata) for (s in unique(d0$stream)) for (co in levels(d0$cohort)) {
    d <- d0[stream == s & cohort == co & !is.na(get(st))]; d[, lev := factor(get(st))]
    keep <- names(which(table(d$lev) >= min_n)); d <- d[lev %in% keep][, lev := droplevels(lev)]
    if (nlevels(d$lev) < 2L) next
    fit <- tryCatch(stats::glm(cbind(n_missing, n_observed) ~ lev, data = d, family = stats::quasibinomial()), error = function(e) NULL)
    if (is.null(fit)) next
    cf <- summary(fit)$coefficients; cf <- cf[-1L, , drop = FALSE]
    an <- tryCatch(stats::anova(fit, test = "F"), error = function(e) NULL)
    rows[[length(rows) + 1L]] <- data.table(
      stratum_var = st, stream = s, cohort = co, level = sub("^lev", "", rownames(cf)), reference = levels(d$lev)[1L],
      n_persons_level = as.integer(table(d$lev)[sub("^lev", "", rownames(cf))]), n_persons_reference = sum(d$lev == levels(d$lev)[1L]),
      estimate = cf[, 1L], se = cf[, 2L], odds_ratio = exp(cf[, 1L]), or_low = exp(cf[, 1L] - 1.96 * cf[, 2L]), or_high = exp(cf[, 1L] + 1.96 * cf[, 2L]),
      p = cf[, 4L], dispersion = summary(fit)$dispersion,
      omnibus_F = if (!is.null(an)) an[["F"]][2L] else NA_real_, omnibus_df = nlevels(d$lev) - 1L, omnibus_p = if (!is.null(an)) an[["Pr(>F)"]][2L] else NA_real_)
  }
  out <- rbindlist(rows, fill = TRUE)
  if (nrow(out) > 0L) {
    om <- unique(out[, .(stratum_var, stream, cohort, omnibus_p)])[, omnibus_p_holm := stats::p.adjust(omnibus_p, "holm"), by = .(stratum_var, stream)]
    out <- merge(out, om[, .(stratum_var, stream, cohort, omnibus_p_holm)], by = c("stratum_var", "stream", "cohort"), sort = FALSE)
  }
  out[]
}

# All four strata: per-cell summary of the person-level % missing and the tests.
missing_pct_by_stratum <- function(pg, persons, min_n = 5L) {
  long <- missing_pct_long(pg, persons)
  strata <- c("sex", "age_group", "race_plot", "cci_group")
  list(long = long,
       summary = rbindlist(lapply(strata, function(st) cbind(stratum_var = st, setnames(
         metric_summary_table(long, st), st, "level"))), fill = TRUE),
       tests = rbindlist(lapply(strata, function(st) missing_pct_stratum_tests(long, st, min_n)), fill = TRUE),
       models = missing_pct_stratum_models(pg, persons, min_n))
}
