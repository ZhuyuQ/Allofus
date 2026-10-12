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
