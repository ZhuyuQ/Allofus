# ==============================================================================
# covariates.R
# ------------------------------------------------------------------------------
# The person-level covariates of EVERY model (gap onset, models B and C, sensitivity, run-in analysis), built in one place so that all
# models code them identically.
#
#   cohort       factor, reference Healthy (cfg$reference_cohort)
#   age10        age at the FIRST record in years / 10: per 10 years, NOT centred (the intercept refers to age 0 and is never interpreted)
#   sex2         Female (reference) / Male
#   cci          Charlson comorbidity SCORE, per +1 point. A person without a score gets the median and cci_missing = 1, so nobody is dropped
#                and the filler value does not influence the CCI coefficient. The score is strongly tied to the cohorts (it counts diabetes and
#                cancer), which is why its multicollinearity (GVIF) is always reported and three sensitivity analyses exist (see below).
#   cci_wc       sensitivity only: the score minus the mean score of the person's cohort (the part of the score that is NOT cohort membership)
#   cci_grp      sensitivity only: CCI groups 0 / 1-2 / 3-4 / >=5 (+ "Unknown"); reference = the most common group, absent groups drop out
#   device       the person's MAIN Fitbit device family (device.R): "Missing" (no usable device), the wrist-worn families, "Clip-on or phone".
#                Product families, not versions: families run different algorithms and have different capabilities (e.g. heart rate), which
#                changes what is recorded. Wrist families with fewer than cfg$device$min_family_people people are pooled as "Other wrist-worn".
#                Reference = the largest family that is not "Missing" (cfg$device$reference overrides). Time-invariant (main device).
# cfg$model$adjust_cci / adjust_device switch the CCI and device terms off.
# ==============================================================================

# Factor of the device family for the models: small wrist families pooled, reference first, "Missing" last.
device_model_factor <- function(device_family, cfg) {
  miss <- cfg$device$missing_label
  x <- as.character(device_family); x[is.na(x)] <- miss
  tab <- table(x)
  special <- c(miss, cfg$device$nonwrist_label, cfg$device$other_label)
  small <- setdiff(names(tab)[tab < (cfg$device$min_family_people %||% 50L)], special)
  if (length(small)) x[x %in% small] <- cfg$device$other_label
  tab <- table(x)
  known <- setdiff(names(tab), miss)
  ref <- if (length(known)) { r <- cfg$device$reference; if (!is.null(r) && r %in% known) r else names(which.max(tab[known])) } else miss
  lv <- c(ref, setdiff(names(sort(tab[known], decreasing = TRUE)), ref), if (miss %in% names(tab) && !identical(ref, miss)) miss)
  f <- factor(x, levels = lv)
  attr(f, "reference") <- ref
  f
}

# The reference family used in the models (for labels): the same rule as device_model_factor().
device_reference <- function(persons, cfg) {
  if (!"device_family" %in% names(persons)) return(NA_character_)
  attr(device_model_factor(persons$device_family, cfg), "reference")
}

# CCI groups for the sensitivity analysis: reference first (the most common group), "Unknown" last; groups nobody has are dropped.
cci_group_factor <- function(cci_group) {
  g <- as.character(cci_group); g[is.na(g)] <- "Unknown"
  tab <- table(g[g != "Unknown"])
  ref <- if (length(tab)) names(which.max(tab)) else "Unknown"
  lv <- c(ref, setdiff(c("0", "1-2", "3-4", ">=5"), ref), "Unknown")
  factor(g, levels = lv[lv %in% unique(g)])
}

person_covariates <- function(persons, cfg) {
  p <- data.table::as.data.table(persons)
  out <- p[, .(person_id, cohort = cohort_factor(cohort, cfg), sex2 = factor(as.character(sex), levels = c("Female", "Male")), age10 = age_first / 10)]
  if (isTRUE(cfg$model$adjust_cci)) {
    if (!"cci_score" %in% names(p)) stop("the person table has no cci_score (CCI file missing?): set cfg$model$adjust_cci <- FALSE or load the CCI file", call. = FALSE)
    sc <- p$cci_score; med <- stats::median(sc, na.rm = TRUE); if (!is.finite(med)) med <- 0
    cmean <- tapply(sc, as.character(p$cohort), mean, na.rm = TRUE)
    out[, `:=`(cci = ifelse(is.na(sc), med, sc), cci_missing = as.integer(is.na(sc)),
               cci_wc = ifelse(is.na(sc), 0, sc - cmean[as.character(p$cohort)]), cci_grp = cci_group_factor(p$cci_group))]
  }
  if (isTRUE(cfg$model$adjust_device)) {
    if (!"device_family" %in% names(p))
      stop("no device information for the models: run Cell 6c (export) and Cell 8b (run_etl_device) first, or set cfg$model$adjust_device <- FALSE", call. = FALSE)
    out[, device := device_model_factor(p$device_family, cfg)]
  }
  out[]
}
