# ==============================================================================
# wear_consistency.R
# ------------------------------------------------------------------------------
# Stratified analysis at the end of the workflow: do the findings differ between people who wore the device
# consistently and people who wore it inconsistently?
#
# MEASURE (reproducible from the day tables only)
#   wear consistency of a person in a stream = OBSERVED-DAY COVERAGE
#       coverage = observed days / days between the person's first and last record in that stream       (0 < coverage <= 1)
#   It is the complement of the "% of days missing" of the missingness summaries and uses the same observation window
#   (first to last record), so it needs no information beyond the data the pipeline already reads. People were selected
#   upstream on >= 70% observed days (cfg$selection), so coverage lies in [0.70, 1] and "inconsistent" means "more gaps
#   than the typical participant", not "barely wore it".
# CLASSIFICATION RULE (fixed in cfg$wear BEFORE any stratified model is fitted; written to wear_consistency_classification)
#   method "median" (default): Consistent = coverage >= the stream's median coverage over the analysed people, else
#   Inconsistent (two strata of similar size; the cut-off is a property of the sample, not tuned to any result);
#   method "fixed": Consistent = coverage >= cfg$wear$threshold. The stratum is defined in the stream that is modelled
#   (the stream whose gaps are the outcome).
# WHAT THE STRATIFICATION CAN AND CANNOT SHOW
#   Coverage is computed from the same observed days that define the gap outcome. A person with many gaps has a lower
#   coverage and contributes more events, so the strata differ by construction in how often gaps occur; the stratified
#   models compare the WITHIN-person association (prior-day steps -> gap tomorrow) between people with different wear
#   patterns. A difference is effect modification by wear pattern, not a causal effect of consistent wear, and the strata
#   are not exchangeable (see the composition table). Strata are compared with a formal test, never by which stratum
#   has a significant estimate.
# ==============================================================================

WEAR_LEVELS <- c("Consistent", "Inconsistent")

# ---- classification --------------------------------------------------------------------------------

# Person-level coverage from the person-level gap summary (summarise_person_gaps): person_id, n_observed, window_days, ...
wear_coverage_table <- function(pg) {
  out <- pg[, .(person_id, window_days, n_observed, coverage = n_observed / window_days,
                n_gaps, gap_max, longest_streak)]
  out[]
}

# Add the stratum to a coverage table. attr(, "criteria") records the rule and the realized cut-off.
classify_wear_consistency <- function(cov_tab, cfg) {
  method <- cfg$wear$method
  stopifnot(method %in% c("median", "fixed"))
  thr <- if (method == "fixed") cfg$wear$threshold else stats::median(cov_tab$coverage)
  out <- copy(cov_tab)
  out[, stratum := factor(fifelse(coverage >= thr, "Consistent", "Inconsistent"), levels = WEAR_LEVELS)]
  setattr(out, "criteria", data.table(
    method = method, threshold = thr, n_people = nrow(out),
    n_consistent = sum(out$stratum == "Consistent"), n_inconsistent = sum(out$stratum == "Inconsistent"),
    rule = sprintf("Consistent: observed-day coverage >= %.4f (%s); Inconsistent: below", thr,
                   if (method == "median") "median of the analysed people" else "fixed threshold")))
  out[]
}

# ---- descriptive statistics -------------------------------------------------------------------------------

smd_continuous <- function(x, g) {
  a <- x[g == WEAR_LEVELS[1L]]; b <- x[g == WEAR_LEVELS[2L]]
  s <- sqrt((stats::var(a, na.rm = TRUE) + stats::var(b, na.rm = TRUE)) / 2)
  if (!is.finite(s) || s == 0) return(NA_real_)
  (mean(a, na.rm = TRUE) - mean(b, na.rm = TRUE)) / s
}
smd_binary <- function(x, g) {
  p1 <- mean(x[g == WEAR_LEVELS[1L]], na.rm = TRUE); p2 <- mean(x[g == WEAR_LEVELS[2L]], na.rm = TRUE)
  s <- sqrt((p1 * (1 - p1) + p2 * (1 - p2)) / 2)
  if (!is.finite(s) || s == 0) return(NA_real_)
  (p1 - p2) / s
}

# Descriptive tables of the wear-consistency strata, one stream at a time.
#   cov_tab : classify_wear_consistency() output of ONE stream          persons : the person table
# Returns list(by_person, summary, composition, classification).
wear_consistency_describe <- function(cov_tab, persons, stream, cfg) {
  d <- merge(cov_tab, persons[, .(person_id, cohort, age_first, sex, cci_score)], by = "person_id", all.x = TRUE)
  d[, female := as.integer(sex == "Female")]
  one <- function(x, label) {
    data.table(stream = stream, group = label, n_people = nrow(x),
               coverage_median = stats::median(x$coverage), coverage_q25 = stats::quantile(x$coverage, 0.25, names = FALSE),
               coverage_q75 = stats::quantile(x$coverage, 0.75, names = FALSE),
               follow_up_days_median = stats::median(x$window_days), observed_days_median = stats::median(x$n_observed),
               gaps_per_person_median = stats::median(x$n_gaps), gaps_total = sum(x$n_gaps),
               gaps_per_100_observed_days = 100 * sum(x$n_gaps) / sum(x$n_observed),
               pct_with_any_gap = 100 * mean(x$n_gaps > 0), longest_gap_median = stats::median(x$gap_max),
               longest_streak_median = stats::median(x$longest_streak),
               age_mean = mean(x$age_first, na.rm = TRUE), age_sd = stats::sd(x$age_first, na.rm = TRUE),
               pct_female = 100 * mean(x$female, na.rm = TRUE), cci_median = stats::median(x$cci_score, na.rm = TRUE))
  }
  summ <- rbindlist(c(list(one(d, "All")), lapply(WEAR_LEVELS, function(lv) one(d[stratum == lv], lv))))
  coh <- dcast(d[, .N, by = .(group = as.character(stratum), cohort)], group ~ cohort, value.var = "N", fill = 0L)
  coh <- rbind(coh, cbind(group = "All", as.data.table(as.list(colSums(coh[, -1L, with = FALSE])))))
  coh_cols <- setdiff(names(coh), "group"); tot <- rowSums(coh[, coh_cols, with = FALSE])
  coh_pct <- copy(coh)
  coh_pct[, (coh_cols) := lapply(.SD, function(v) 100 * v / tot), .SDcols = coh_cols]
  setnames(coh_pct, coh_cols, paste0("pct_", coh_cols))
  summ <- merge(summ, coh_pct, by = "group", all.x = TRUE, sort = FALSE)
  # standardized mean differences (Consistent - Inconsistent); |SMD| > 0.1 is read as a relevant imbalance
  comp <- rbindlist(c(list(
    data.table(variable = "Age at first record (years)", smd = smd_continuous(d$age_first, d$stratum)),
    data.table(variable = "Female (%)", smd = smd_binary(d$female, d$stratum)),
    data.table(variable = "CCI score", smd = smd_continuous(d$cci_score, d$stratum)),
    data.table(variable = "Follow-up length (days)", smd = smd_continuous(d$window_days, d$stratum))),
    lapply(levels(d$cohort), function(co) data.table(variable = paste("Cohort:", co), smd = smd_binary(as.integer(d$cohort == co), d$stratum)))))
  comp[, `:=`(stream = stream, abs_smd = abs(smd), imbalanced = abs(smd) > 0.1)]
  setcolorder(comp, c("stream", "variable"))
  list(by_person = d[, .(stream = stream, person_id, cohort, stratum, coverage, window_days, n_observed, n_gaps, gap_max, longest_streak)],
       summary = summ, composition = comp, classification = cbind(stream = stream, attr(cov_tab, "criteria")))
}

# ---- stratified models, as independent checkpointed units -------------------------------------------------------

# The units that fit the applicable models within each stratum of ONE stream (same specifications as the primary
# analysis, only the people differ):
#   wear_onset : the onset spec(s) of cfg$wear$onset_roles (default: the primary joint adjusted >= 1-day model), per stratum
#   wear_extra : the likelihood-ratio test of (primary term x stratum), with the stratum main effect in both models
#   wear_gapB  : model B (adjusted, first full day back) per stratum
#   wear_gapC  : model C (adjusted) per stratum
# map : data.table(person_id, stratum). bundle / gap_data : lazy objects (built only when some unit has to run).
run_wear_stage <- function(s, bundle, gap_data, grid, map, cfg, ckpt_dir, dkey) {
  status <- character()
  crit <- list(method = cfg$wear$method, threshold = signif(attr(map, "criteria")$threshold, 10),
               people = ckpt_hash(list(consistent = sort(map$person_id[map$stratum == "Consistent"]),
                                       inconsistent = sort(map$person_id[map$stratum == "Inconsistent"]))))
  n_by <- table(map$stratum)
  usable <- names(n_by)[n_by >= cfg$wear$min_persons_stratum]
  for (lv in setdiff(WEAR_LEVELS, usable)) status[paste0("skipped: ", lv)] <- "too few people"
  if (length(usable) == 0L) return(status)
  unit <- function(kind, id, compute, spec_extra = NULL) {
    u <- run_unit(ckpt_dir, list(kind = kind, stream = s, id = id),
                  ckpt_key(cfg, kind, dkey, c(spec_extra, list(criteria = crit, schema = 1L)), cfg$model$sample_persons), compute, cfg)
    status[paste(kind, id, sep = ": ")] <<- u$status
  }

  g <- grid[stream == s & role %in% cfg$wear$onset_roles]
  for (i in seq_len(nrow(g))) {
    spec <- g[i]
    for (lv in usable) {
      unit("wear_onset", paste(spec$spec_id, lv, sep = " | "), function() {
        r <- fit_onset_spec(bundle$get(), spec, cfg, variant = list(strata = map, stratum = lv)); r$fit <- NULL; r
      }, c(spec_key(spec), list(stratum = lv)))
    }
    if (length(usable) == 2L) {
      unit("wear_extra", paste("interaction", spec$spec_id, sep = " | "), function()
        cbind(stream = s, spec_id = spec$spec_id, test_group_interaction(bundle$get(), spec, cfg, by = "stratum", variant = list(strata = map))),
        c(spec_key(spec), list(test = "stratum_interaction")))
    }
  }

  if (!is.null(gap_data)) {
    for (lv in usable) {
      ids_lv <- map$person_id[map$stratum == lv]
      unit("wear_gapB", paste("B_back2_adj", lv, sep = " | "), function() {
        gdat <- gap_data$get()[person_id %in% ids_lv]
        isolated(cfg, fit_gap_dose_response(gdat, cfg, "back2_dev", adjusted = TRUE, adjust_pre = FALSE, id = "B_back2_adj"))
      }, list(model = "B_back2_adj", stratum = lv))
      unit("wear_gapC", paste("C_adj", lv, sep = " | "), function() {
        gdat <- gap_data$get()[person_id %in% ids_lv]
        isolated(cfg, fit_gap_length_ordinal(gdat, cfg, adjusted = TRUE, id = "C_adj"))
      }, list(model = "C_adj", stratum = lv))
    }
  }
  status
}

# Compare the stratum-specific estimates of the SAME term. The strata contain different people, so the two estimates are
# independent and  z = (b_c - b_i) / sqrt(se_c^2 + se_i^2)  tests "same effect in both strata" (a Wald test of the difference,
# equivalent to the interaction term in a model with separate error structures). Reference stratum: Consistent.
# For log-odds models the exponentiated difference is the RATIO OF ODDS RATIOS (Inconsistent / Consistent).
compare_strata <- function(b1, se1, b2, se2, log_scale = TRUE) {
  d <- b2 - b1; se <- sqrt(se1^2 + se2^2); z <- d / se
  p <- 2 * stats::pnorm(abs(z), lower.tail = FALSE)
  out <- data.table(difference = d, difference_low = d - 1.96 * se, difference_high = d + 1.96 * se, difference_se = se, z = z, p_wald = p)
  if (log_scale) out[, `:=`(ratio_of_or = exp(difference), ratio_low = exp(difference_low), ratio_high = exp(difference_high))]
  out[]
}

# Everything finished, all streams: coefficient tables, model information, diagnostics, the formal comparison and the
# sparse-data flags.
collect_wear_results <- function(ckpt_dir, cfg, data_keys) {
  stratum_of <- function(ids) sub("^.* \\| ", "", ids)
  grab <- function(cu, part, add_stratum = TRUE) rbindlist(Map(function(v, st, id) {
    x <- v[[part]]; if (is.null(x) || !is.data.frame(x) || nrow(x) == 0L) return(NULL)
    if (!"stream" %in% names(x)) x <- cbind(stream = st, x); if (add_stratum && !"wear_stratum" %in% names(x)) x[, wear_stratum := stratum_of(id)]
    x[]
  }, cu$values, cu$units$stream, cu$units$id), fill = TRUE)
  on <- collect_units(ckpt_dir, "wear_onset", cfg, data_keys); ex <- collect_units(ckpt_dir, "wear_extra", cfg, data_keys)
  gb <- collect_units(ckpt_dir, "wear_gapB", cfg, data_keys); gc <- collect_units(ckpt_dir, "wear_gapC", cfg, data_keys)
  out <- list(
    wear_strata_onset_coefficients = grab(on, "coef"), wear_strata_onset_models = grab(on, "info"),
    wear_strata_onset_collinearity = grab(on, "gvif"), wear_strata_onset_term_tests = grab(on, "wald"),
    wear_strata_interaction = rbindlist(ex$values, fill = TRUE),
    wear_strata_gapB = grab(gb, "coef"), wear_strata_gapB_all_terms = grab(gb, "coef_all"),
    wear_strata_gapC = grab(gc, "coef"), wear_strata_gapC_all_terms = grab(gc, "coef_all"),
    wear_strata_gap_collinearity = rbindlist(list(grab(gb, "gvif"), grab(gc, "gvif")), fill = TRUE),
    wear_failures = rbindlist(list(on$failures, ex$failures, gb$failures, gc$failures), fill = TRUE))
  out$wear_strata_comparison <- wear_strata_comparison(out, cfg)
  out
}

# One row per (model, stream, term of interest): both strata side by side, the difference between them and the sparse-data
# flags. "Term of interest": the primary within-person term of the onset model (OR per +1 SD higher steps), the gap-length
# term of model B (change in activity per doubling of the gap, SD units) and the prior-day steps term of model C (OR of a
# longer gap per +1 SD higher steps).
wear_strata_comparison <- function(res, cfg) {
  rows <- list()
  add <- function(model, tab, which_term, effect_is_or) {
    if (is.null(tab) || nrow(tab) == 0L) return(invisible(NULL))
    t1 <- tab[term == which_term & !is.na(estimate)]
    if (nrow(t1) == 0L) return(invisible(NULL))
    for (s in unique(t1$stream)) {
      a <- t1[stream == s & wear_stratum == WEAR_LEVELS[1L]]; b <- t1[stream == s & wear_stratum == WEAR_LEVELS[2L]]
      if (nrow(a) != 1L || nrow(b) != 1L) next
      cmp <- compare_strata(a$estimate, a$se, b$estimate, b$se, log_scale = effect_is_or)
      one <- function(x) if (effect_is_or) c(exp(x$estimate), exp(x$conf_low), exp(x$conf_high)) else c(x$estimate, x$conf_low, x$conf_high)
      ea <- one(a); eb <- one(b)
      rows[[length(rows) + 1L]] <<- cbind(data.table(
        model = model, stream = s, term = which_term, effect_measure = if (effect_is_or) "odds ratio per +1 SD" else "SD change per unit",
        effect_consistent = ea[1L], low_consistent = ea[2L], high_consistent = ea[3L],
        effect_inconsistent = eb[1L], low_inconsistent = eb[2L], high_inconsistent = eb[3L]), cmp)
    }
  }
  add("onset (primary)", if (nz(res$wear_strata_onset_coefficients)) res$wear_strata_onset_coefficients[role == "primary"], cfg$primary_term, TRUE)
  add("model B (adjusted)", if (nz(res$wear_strata_gapB)) res$wear_strata_gapB[term == "log2_gap"], "log2_gap", FALSE)
  add("model C (adjusted)", res$wear_strata_gapC, "w_pre_day_steps", TRUE)
  out <- rbindlist(rows, fill = TRUE)
  if (nrow(out) == 0L) return(out)
  # likelihood-ratio test of the interaction in the pooled onset model (the formal test for the onset question)
  it <- res$wear_strata_interaction
  if (!is.null(it) && nrow(it) > 0L) {
    out[, p_interaction_lrt := NA_real_]
    for (s in unique(it$stream)) out[model == "onset (primary)" & stream == s, p_interaction_lrt := it[stream == s & term == cfg$primary_term, p][1L]]
  }
  out[, p_difference := if ("p_interaction_lrt" %in% names(out)) fifelse(!is.na(p_interaction_lrt), p_interaction_lrt, p_wald) else p_wald]
  out[]
}

# Sparse-data and convergence flags per stratified fit: events (onset) or gaps (B, C) per estimated coefficient, status of the
# fit, and whether the stratum's estimate is usable. events per variable below cfg$wear$epv_min = "sparse".
wear_strata_quality <- function(res, cfg) {
  rows <- list()
  if (nz(res$wear_strata_onset_models)) {
    m <- res$wear_strata_onset_models
    npar <- res$wear_strata_onset_coefficients[term != "(Intercept)", .(n_coefficients = .N), by = .(stream, spec_id, wear_stratum)]
    m <- merge(m, npar, by = c("stream", "spec_id", "wear_stratum"), all.x = TRUE)
    rows[[1L]] <- m[, .(model = paste0("onset: ", spec_id), stream, stratum = wear_stratum, n_people = n_persons, n_events = n_events,
                        n_coefficients, events_per_coefficient = n_events / n_coefficients, status,
                        singular = if ("singular" %in% names(m)) singular else NA)]
  }
  for (nm in c("wear_strata_gapB", "wear_strata_gapC")) {
    m <- res[[nm]]; ta <- res[[paste0(nm, "_all_terms")]]
    if (!nz(m) || !nz(ta)) next
    npar <- ta[, .(n_coefficients = .N), by = .(stream, model_id, wear_stratum)]
    m <- merge(m, npar, by = c("stream", "model_id", "wear_stratum"), all.x = TRUE)
    rows[[length(rows) + 1L]] <- m[, .(model = paste0(if (nm == "wear_strata_gapB") "model B: " else "model C: ", model_id), stream, stratum = wear_stratum,
                                       n_people = n_persons, n_events = n_gaps, n_coefficients, events_per_coefficient = n_gaps / n_coefficients, status,
                                       singular = if ("singular" %in% names(m)) singular else NA)]
  }
  out <- rbindlist(rows, fill = TRUE)
  if (nrow(out) == 0L) return(out)
  out[, sparse := events_per_coefficient < cfg$wear$epv_min]
  out[, usable := status %in% c("ok", "ok_bootstrap_se") & !sparse]
  out[]
}
