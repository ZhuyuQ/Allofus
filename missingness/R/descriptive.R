# ==============================================================================
# descriptive.R
# ------------------------------------------------------------------------------
# Summary statistics of the METRIC VALUES (steps, sleep minutes, heart-rate zones ...)
# for each stream, by cohort and by strata (sex, age group, race, CCI group).
# Follows the earlier convention: first a per-person mean/SD/CV, then the
# distribution of person means across people.
# ==============================================================================

# Which metric columns to describe per stream (columns that were dropped by the ETL
# are skipped automatically).
descriptive_metrics <- function() {
  list(
    activity = c("steps", "calories_out", "sedentary_minutes", "lightly_active_minutes",
                 "fairly_active_minutes", "very_active_minutes", "mvpa_minutes"),
    sleep    = c("minute_in_bed", "minute_asleep", "sleep_efficiency", "minute_awake",
                 "minute_deep", "minute_light", "minute_rem", "minute_wake"),
    hr       = c("total_hr_min", "elevated_min", "mins_out_of_range", "mins_fat_burn",
                 "mins_cardio", "mins_peak", "total_cal", "cal_out_of_range", "cal_fat_burn",
                 "cal_cardio", "cal_peak", "hrhi_out_of_range"))
}

# Per-person mean, SD and CV (%) of each metric. Wide table: person_id, n_days,
# <metric>__mean, <metric>__sd, <metric>__cv.
person_metric_stats <- function(days, metrics) {
  m <- intersect(metrics, names(days))
  if (length(m) == 0L) stop("none of the requested metrics are in the day table", call. = FALSE)
  mean_dt <- days[, c(list(n_days = .N), lapply(.SD, function(x) mean(x, na.rm = TRUE))),
                  by = person_id, .SDcols = m]
  sd_dt   <- days[, lapply(.SD, function(x) stats::sd(x, na.rm = TRUE)), by = person_id, .SDcols = m]
  setnames(mean_dt, m, paste0(m, "__mean"))
  setnames(sd_dt,   m, paste0(m, "__sd"))
  out <- merge(mean_dt, sd_dt, by = "person_id")
  for (cl in m) {
    mu <- out[[paste0(cl, "__mean")]]; s <- out[[paste0(cl, "__sd")]]
    set(out, j = paste0(cl, "__cv"), value = fifelse(!is.na(mu) & mu > 0, s / mu * 100, NA_real_))
  }
  for (cl in names(out)) set(out, which(is.nan(out[[cl]])), cl, NA_real_)
  out[]
}

# Long format of the person means (one row per person x metric), with the strata
# columns attached - the input to the plots and summary tables.
person_means_long <- function(pms, persons, stream) {
  mean_cols <- grep("__mean$", names(pms), value = TRUE)
  long <- melt(pms[, c("person_id", mean_cols), with = FALSE], id.vars = "person_id",
               variable.name = "metric", value.name = "value")
  long[, metric := sub("__mean$", "", as.character(metric))]
  long <- merge(long, persons[, .(person_id, cohort, sex, age_group, race_plot, cci_group)],
                by = "person_id")
  long[, stream := stream]
  long[]
}

# N, mean, SD, median, quartiles, P5 and P95 of the person means, by cohort and
# (optionally) by one stratum variable.
metric_summary_table <- function(long, stratum = NULL) {
  by <- c("stream", "metric", "cohort", stratum)
  long[!is.na(value), .(N = .N, Mean = mean(value), SD = stats::sd(value),
                        Median = stats::median(value),
                        Q1 = stats::quantile(value, 0.25, names = FALSE),
                        Q3 = stats::quantile(value, 0.75, names = FALSE),
                        P5 = stats::quantile(value, 0.05, names = FALSE),
                        P95 = stats::quantile(value, 0.95, names = FALSE)), by = by][order(stream, metric, cohort)]
}

# Optional small-cell suppression (off by default, cfg$suppress_n_lt = 0).
suppress_small_cells <- function(tab, cfg, count_col = "N") {
  if (cfg$suppress_n_lt <= 0L) return(tab)
  out <- copy(tab)
  small <- out[[count_col]] < cfg$suppress_n_lt
  num <- setdiff(names(out)[vapply(out, is.numeric, logical(1))], count_col)
  for (cl in num) set(out, which(small), cl, NA_real_)
  set(out, which(small), count_col, NA_integer_)
  out[]
}
