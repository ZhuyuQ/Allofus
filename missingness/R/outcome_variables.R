# ==============================================================================
# outcome_variables.R
# ------------------------------------------------------------------------------
# How does the PERSON-LEVEL outcome "% of days missing" (one value per person and stream, summarise_person_gaps()) relate to each
# variable of the models, one variable at a time? Descriptive and associational; the gap-onset / gap-level models are fitted elsewhere.
#
# Variables (person level): cohort, sex, watch version (device family, "Missing" = no watch version, see device.R), age at the first
# record, and the between-person level of every model predictor of the stream (the person's mean over observed days: steps, sedentary
# minutes, MVPA, time in bed, sleep efficiency, heart-rate minutes; log1p-transformed ones are shown as geometric means).
# Day-level variables of the models (weekday, study time, calendar year) are shown at the day level in the temporal figures
# (missing_by_weekday / study_week / year); the gap-history terms (previous gap, any earlier gap, days since return) are built from the
# outcome itself and are not compared with it here.
#
# Tests (every one is on one value per person, so people with long follow-up do not dominate):
#   categorical : Kruskal-Wallis (>= 3 levels, effect epsilon-squared) or Wilcoxon rank-sum (2 levels, Cliff's delta) with a percentile
#                 bootstrap 95% interval of the effect size; each level against the reference level: Wilcoxon with the Hodges-Lehmann
#                 difference in percentage points and its 95% interval, Cliff's delta, Holm over the comparisons of one variable.
#   continuous  : Spearman rank correlation with the Bonett-Wright 95% interval.
#   all variables of one stream: Holm adjustment of the overall p values (one family per stream).
#   adjusted    : quasi-binomial model  cbind(missing days, observed days) ~ variable + the other person-level covariates
#                 (cohort, sex, age): odds ratio that a day is missing, per level against the reference (categorical) or per +1 SD
#                 (continuous), with the F test of the variable; the unadjusted fit is reported next to it.
# Levels with fewer than cfg$outcome_vars$min_n people are drawn but left out of the tests.
# ==============================================================================

OV_CATEGORICAL <- c(cohort = "Cohort", sex = "Sex", watch = "Watch version")

# Display names / units of the continuous variables (model predictors): "gm" = geometric mean (the predictor is log1p-transformed).
OV_PRED_LABELS <- c(act_steps = "Steps per day", act_sedentary = "Sedentary minutes per day", act_mvpa = "MVPA minutes per day",
                    sleep_inbed = "Time in bed (minutes)", sleep_eff = "Sleep efficiency (%)",
                    hr_cover = "Heart-rate wear minutes per day", hr_elev = "Elevated heart-rate minutes per day")

# Person means of the model predictors on the model's own (transformed) scale, over the person's observed days.
#   inputs: stream_inputs() (day tables by stream); preds: predictor names. Returns person_id + one column per predictor
#   (the mean of the transformed value) and <pred>__display (shown value: geometric mean for log1p predictors, else the mean).
person_predictor_means <- function(inputs, preds, cfg) {
  out <- NULL
  for (p in preds) {
    sp <- cfg$predictors[[p]]; tbl <- inputs[[sp$stream]]
    if (is.null(tbl) || !sp$col %in% names(tbl)) next
    x <- tbl[[sp$col]]
    x <- if (identical(sp$transform, "log1p")) log1p(pmax(x, 0)) else x
    m <- data.table(person_id = tbl$person_id, x = x)[!is.na(x), .(m = mean(x)), by = person_id]
    setnames(m, "m", p)
    m[, (paste0(p, "__display")) := if (identical(sp$transform, "log1p")) expm1(get(p)) else get(p)]
    out <- if (is.null(out)) m else merge(out, m, by = "person_id", all = TRUE)
  }
  if (is.null(out)) data.table(person_id = numeric()) else out[]
}

# One row per person: outcome, categorical variables and the continuous variables of this stream.
outcome_variable_data <- function(pg, persons, device_person, means, cfg) {
  d <- merge(pg[, .(person_id, n_missing, n_observed, window_days, pct_missing)],
             data.table::as.data.table(persons)[, .(person_id, cohort, sex, age_first)], by = "person_id")
  d[, cohort := factor(as.character(cohort), levels = cfg$cohorts)]
  if (!is.null(device_person)) {
    d <- merge(d, data.table::as.data.table(device_person)[, .(person_id, watch = device_family)], by = "person_id", all.x = TRUE)
    d[is.na(watch), watch := cfg$device$missing_label]               # nobody is dropped: no device record = missing
    known <- d[watch != cfg$device$missing_label, .N, by = watch][order(-N), watch]
    known <- c(setdiff(known, cfg$device$other_label), intersect(known, cfg$device$other_label))
    d[, watch := factor(watch, levels = c(cfg$device$missing_label, known))]
  }
  if (nrow(means)) d <- merge(d, means, by = "person_id", all.x = TRUE)
  d[]
}

# Bonett-Wright interval of a Spearman correlation.
spearman_ci <- function(rho, n, level = 0.95) {
  if (!is.finite(rho) || n < 5L) return(c(NA_real_, NA_real_))
  z <- atanh(rho); se <- sqrt((1 + rho^2 / 2) / (n - 3)); q <- stats::qnorm(1 - (1 - level) / 2)
  tanh(c(z - q * se, z + q * se))
}

# Percentile bootstrap interval of the effect size of a rank test on `value` across the groups `g` (people resampled).
boot_rank_effect <- function(value, g, B, level = 0.95) {
  g <- droplevels(factor(g)); n <- length(value); k <- nlevels(g)
  if (k < 2L || B < 2L) return(c(NA_real_, NA_real_))
  eff <- function(v, gg) {
    if (k == 2L) cliffs_delta(v[gg == levels(g)[1L]], v[gg == levels(g)[2L]])
    else { kw <- suppressWarnings(stats::kruskal.test(v, gg)); unname(kw$statistic / (length(v) - 1L)) }
  }
  b <- vapply(seq_len(B), function(i) { ix <- sample.int(n, n, replace = TRUE); gg <- g[ix]; if (nlevels(droplevels(gg)) < k) NA_real_ else eff(value[ix], gg) }, numeric(1))
  unname(stats::quantile(b, c((1 - level) / 2, 1 - (1 - level) / 2), na.rm = TRUE))
}

# Quasi-binomial fit of cbind(missing, observed) ~ rhs on d (one row per person); NULL if it cannot be fitted.
qb_fit <- function(d, rhs) tryCatch(suppressWarnings(stats::glm(stats::as.formula(paste("cbind(n_missing, n_observed) ~", rhs)), data = d, family = stats::quasibinomial())), error = function(e) NULL)

# Odds-ratio table of the coefficient rows of a quasi-binomial fit whose names start with `prefix`.
qb_rows <- function(fit, prefix, d, adjustment, extra = list()) {
  cf <- summary(fit)$coefficients; cf <- cf[startsWith(rownames(cf), prefix), , drop = FALSE]
  if (!nrow(cf)) return(NULL)
  dfres <- fit$df.residual; q <- stats::qt(0.975, dfres)
  data.table(term = sub(paste0("^", prefix), "", rownames(cf)), adjustment = adjustment, odds_ratio = exp(cf[, 1L]),
             or_low = exp(cf[, 1L] - q * cf[, 2L]), or_high = exp(cf[, 1L] + q * cf[, 2L]), p = cf[, 4L],
             dispersion = summary(fit)$dispersion, n_persons = nrow(d))
}

# The analysis of ONE stream. Returns list(data =, tests =, contrasts =, summary =, models =, bins =).
outcome_variable_analysis <- function(d, stream, cfg, cont_vars = NULL) {
  oc <- cfg$outcome_vars; min_n <- oc$min_n
  d <- copy(d); d[, age10 := age_first / 10]
  cats <- intersect(names(OV_CATEGORICAL), names(d))
  tests <- list(); contrasts <- list(); summ <- list(); models <- list()

  # reference levels
  ref_of <- list(cohort = cfg$reference_cohort %||% cfg$cohorts[1L], sex = "Female",
                 watch = if ("watch" %in% cats) { kn <- d[watch != cfg$device$missing_label, .N, by = watch][order(-N), as.character(watch)]; if (length(kn)) kn[1L] else NA_character_ })
  adj_of <- list(cohort = c("sex", "age10"), sex = c("cohort", "age10"), watch = c("cohort", "sex", "age10"))

  for (v in cats) {
    x <- d[!is.na(get(v)) & !is.na(pct_missing)]
    x[, lev := droplevels(factor(get(v)))]
    ct <- x[, .N, by = lev]
    summ[[v]] <- x[, .(n = .N, mean = mean(pct_missing), median = stats::median(pct_missing), q1 = stats::quantile(pct_missing, 0.25), q3 = stats::quantile(pct_missing, 0.75)), by = lev][
      , `:=`(stream = stream, variable = v, level = as.character(lev))][, lev := NULL][]
    use <- as.character(ct[N >= min_n, lev]); xx <- x[as.character(lev) %in% use][, lev := droplevels(lev)]
    rt <- rank_test_groups(xx$pct_missing, xx$lev, min_n)
    ci <- if (rt$test == "none") c(NA_real_, NA_real_) else boot_rank_effect(xx$pct_missing, xx$lev, oc$boot_B)
    tests[[v]] <- data.table(stream = stream, variable = v, variable_label = OV_CATEGORICAL[[v]], type = "categorical", n = rt$n_total, test = rt$test,
                             statistic = rt$statistic, df = rt$df, p = rt$p, effect_name = rt$effect_name, effect = rt$effect, effect_low = ci[1L], effect_high = ci[2L],
                             note = sprintf("levels with < %d people left out of the test: %s", min_n, paste(setdiff(as.character(ct$lev), use), collapse = ", ")))
    tests[[v]][note == sprintf("levels with < %d people left out of the test: ", min_n), note := ""]
    # each level against the reference
    ref <- ref_of[[v]]
    if (!is.na(ref) && ref %in% use) {
      r <- xx[as.character(lev) == ref, pct_missing]
      cr <- rbindlist(lapply(setdiff(use, ref), function(lv) {
        a <- xx[as.character(lev) == lv, pct_missing]
        w <- suppressWarnings(stats::wilcox.test(a, r, conf.int = TRUE, exact = FALSE))
        data.table(stream = stream, variable = v, level = lv, reference = ref, n_level = length(a), n_reference = length(r),
                   median_level = stats::median(a), median_reference = stats::median(r), hl_difference_pct_points = unname(w$estimate),
                   ci_low = w$conf.int[1L], ci_high = w$conf.int[2L], cliffs_delta = cliffs_delta(a, r), p = w$p.value)
      }))
      if (nrow(cr)) cr[, p_holm := stats::p.adjust(p, "holm")]
      contrasts[[v]] <- cr
      # quasi-binomial models, unadjusted and adjusted
      g <- xx[, c("n_missing", "n_observed", "cohort", "sex", "age10"), with = FALSE]; g[, lev := stats::relevel(xx$lev, ref = ref)]
      adj <- adj_of[[v]]
      g_adj <- copy(g); for (a in intersect(adj, c("cohort", "sex"))) { tt <- table(g_adj[[a]]); keep <- names(tt)[tt >= min_n]; g_adj <- g_adj[as.character(get(a)) %in% keep] }
      g_adj <- g_adj[stats::complete.cases(g_adj)]
      for (a in intersect(adj, c("cohort", "sex"))) g_adj[, (a) := droplevels(factor(get(a)))]
      if (v == "cohort") g_adj[, cohort := NULL]
      f0 <- qb_fit(g, "lev"); rhs_adj <- paste(c("lev", setdiff(adj, if (v == "cohort") "cohort")), collapse = " + ")
      f1 <- if (nlevels(droplevels(g_adj$lev)) >= 2L) qb_fit(g_adj, rhs_adj) else NULL
      for (pair in list(list(f0, g, "unadjusted"), list(f1, g_adj, paste("adjusted for", paste(setdiff(adj, if (v == "cohort") "cohort"), collapse = ", "))))) {
        fit <- pair[[1L]]; if (is.null(fit)) next
        rows <- qb_rows(fit, "lev", pair[[2L]], pair[[3L]])
        if (is.null(rows)) next
        fit_null <- tryCatch(stats::update(fit, . ~ . - lev), error = function(e) NULL)
        an <- if (!is.null(fit_null)) tryCatch(stats::anova(fit_null, fit, test = "F"), error = function(e) NULL)
        rows[, `:=`(stream = stream, variable = v, reference = ref, omnibus_F = if (!is.null(an)) an$F[2L] else NA_real_, omnibus_p = if (!is.null(an)) an$`Pr(>F)`[2L] else NA_real_)]
        models[[length(models) + 1L]] <- rows
      }
    }
  }

  # continuous variables
  bins <- list(); cv_all <- c(list(age_first = list(label = "Age at first record (years)", display = "age_first", log_x = FALSE)), cont_vars)
  for (v in names(cv_all)) {
    info <- cv_all[[v]]
    if (!v %in% names(d) || !info$display %in% names(d)) next
    x <- d[!is.na(get(v)) & !is.na(pct_missing)]
    n <- nrow(x); if (n < 20L || uniqueN(x[[v]]) < 3L) next
    ct <- suppressWarnings(stats::cor.test(x[[v]], x$pct_missing, method = "spearman", exact = FALSE))
    ci <- spearman_ci(unname(ct$estimate), n)
    tests[[v]] <- data.table(stream = stream, variable = v, variable_label = info$label, type = "continuous", n = n, test = "Spearman rank correlation",
                             statistic = unname(ct$statistic), df = NA_real_, p = ct$p.value, effect_name = "Spearman rho", effect = unname(ct$estimate),
                             effect_low = ci[1L], effect_high = ci[2L], note = "")
    # quantile bins for the figure
    br <- unique(stats::quantile(x[[v]], seq(0, 1, length.out = oc$bins + 1L), na.rm = TRUE))
    if (length(br) >= 3L) {
      x[, bin := cut(get(v), br, include.lowest = TRUE, labels = FALSE)]
      bb <- x[, .(n = .N, x_display = as.numeric(stats::median(get(info$display))), mean_pct = mean(pct_missing), sd = stats::sd(pct_missing)), by = bin][order(bin)]
      bb[, `:=`(se = sd / sqrt(n), stream = stream, variable = v, variable_label = info$label, log_x = isTRUE(info$log_x))]
      bb[, `:=`(low = mean_pct - 1.96 * se, high = mean_pct + 1.96 * se)]
      bins[[v]] <- bb[, .(stream, variable, variable_label, log_x, bin, n, x_display, mean_pct, low, high)]
    }
    # per +1 SD, unadjusted and adjusted for the other person-level covariates
    g <- x[, c("n_missing", "n_observed", "cohort", "sex", "age10"), with = FALSE]; g[, z := as.numeric(scale(x[[v]]))]
    g[, sex := droplevels(factor(sex))]; g <- g[sex %in% names(which(table(sex) >= min_n))][, sex := droplevels(sex)]
    adjrhs <- paste(c("z", "cohort", if (uniqueN(g$sex) > 1L) "sex", if (v != "age_first") "age10"), collapse = " + ")
    for (pair in list(list("z", "unadjusted"), list(adjrhs, paste("adjusted for", paste(setdiff(strsplit(adjrhs, " \\+ ")[[1L]], "z"), collapse = ", "))))) {
      fit <- qb_fit(g[stats::complete.cases(g)], pair[[1L]]); if (is.null(fit)) next
      rows <- qb_rows(fit, "z", g[stats::complete.cases(g)], pair[[2L]]); if (is.null(rows)) next
      rows[, `:=`(stream = stream, variable = v, reference = "per +1 SD", term = "per +1 SD", omnibus_F = NA_real_, omnibus_p = NA_real_)]
      models[[length(models) + 1L]] <- rows
    }
  }

  tests <- rbindlist(tests, fill = TRUE)
  if (nrow(tests)) tests[, p_holm := stats::p.adjust(p, "holm")]                          # one family: all variables of this stream
  models <- rbindlist(models, fill = TRUE)
  if (nrow(models)) setcolorder(models, c("stream", "variable", "term", "reference", "adjustment"))
  list(data = d, tests = tests[], contrasts = rbindlist(contrasts, fill = TRUE), summary = rbindlist(summ, fill = TRUE),
       models = models[], bins = rbindlist(bins, fill = TRUE))
}

# Continuous variables of a stream with their labels (the model predictors of the stream; see the header).
outcome_continuous_spec <- function(stream, cfg) {
  preds <- unique(c(cfg$own_predictors[[stream]], cfg$activity_predictors))
  out <- lapply(preds, function(p) list(label = paste0(OV_PRED_LABELS[[p]] %||% nice_label(p), if (identical(cfg$predictors[[p]]$transform, "log1p")) ", geometric mean" else ", mean"),
                                        display = paste0(p, "__display"), log_x = identical(cfg$predictors[[p]]$transform, "log1p")))
  setNames(out, preds)
}
