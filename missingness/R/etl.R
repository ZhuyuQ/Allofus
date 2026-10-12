# ==============================================================================
# etl.R
# ------------------------------------------------------------------------------
# Turn each RAW stream into a compact DAY TABLE: exactly one row per
# (person_id, date) holding only the columns the analyses need.
#
#   activity : already one row/day            -> keep metrics, add mvpa_minutes
#   sleep    : main sleep + naps are separate rows -> collapse to one row/day
#   hr       : four zone rows per day (Out of Range, Fat Burn, Cardio, Peak)
#              -> pivot to one wide row/day and derive total/elevated minutes
#
# All three functions are PURE (data.table in -> data.table out); reading and
# writing is in io.R. Day-table semantics: a row EXISTS  <=>  the day is "observed".
# ==============================================================================

# ---- shared step: drop columns that carry no information ---------------------

# Drop columns that are mostly NA (> na_frac) or constant. Key columns are
# protected. The dropped columns and the reason are stored in
# attr(result, "dropped_columns") so nothing disappears silently.
drop_uninformative_columns <- function(dt, protect = c("person_id", "date"),
                                       na_frac = 0.90) {
  if (is.na(na_frac)) {            # rule switched off (e.g. before stacking cohorts)
    data.table::setattr(dt, "dropped_columns", data.table(column = character(), frac_na = numeric(),
                                                          n_unique = integer(), reason = character()))
    return(dt)
  }
  cand <- setdiff(names(dt), protect)
  info <- rbindlist(lapply(cand, function(col) {
    x <- dt[[col]]
    frac_na <- mean(is.na(x))
    n_unique <- if (frac_na < 1) uniqueN(x[!is.na(x)]) else 0L
    reason <- if (frac_na > na_frac) "mostly NA"
              else if (n_unique <= 1L) "constant"
              else NA_character_
    data.table(column = col, frac_na = frac_na, n_unique = n_unique, reason = reason)
  }))
  dropped <- info[!is.na(reason)]
  if (nrow(dropped) > 0L) dt <- dt[, setdiff(names(dt), dropped$column), with = FALSE]
  data.table::setattr(dt, "dropped_columns", dropped)
  dt
}

# ---- activity -----------------------------------------------------------------

# raw columns used: person_id, date, steps, calories_out, sedentary_minutes,
# lightly_/fairly_/very_active_minutes (+ any others, kept if informative).
standardize_activity <- function(raw, na_frac = 0.90) {
  dt <- data.table::as.data.table(raw)
  assert_columns(dt, c("person_id", "date", "steps"), "raw activity")
  dt <- copy(dt)
  dt[, date := as.Date(date)]
  dt <- dt[!is.na(date)]

  # Activity is expected to be one row per day. If a duplicate exists, keep the
  # row with the largest step count (deterministic) and report how many.
  n_before <- nrow(dt)
  setorderv(dt, c("person_id", "date", "steps"), order = c(1L, 1L, -1L), na.last = TRUE)
  dt <- unique(dt, by = c("person_id", "date"))
  if (nrow(dt) < n_before) {
    log_info("activity: collapsed %d duplicate person-day row(s)", n_before - nrow(dt))
  }

  # moderate-to-vigorous minutes = fairly + very active
  if (all(c("fairly_active_minutes", "very_active_minutes") %in% names(dt))) {
    dt[, mvpa_minutes := row_sum_na(fairly_active_minutes, very_active_minutes)]
  }
  drop_uninformative_columns(dt, na_frac = na_frac)
}

# ---- sleep --------------------------------------------------------------------

# Convert the 'true'/'false' text (or logical / 0-1) is_main_sleep to logical.
as_main_sleep_flag <- function(x) {
  if (is.logical(x)) return(x)
  tolower(as.character(x)) %in% c("true", "t", "1")
}

# raw columns used: person_id, sleep_date, is_main_sleep, minute_* columns.
# Output date = sleep_date (the date the sleep record is assigned to).
standardize_sleep <- function(raw, na_frac = 0.90) {
  dt <- data.table::as.data.table(raw)
  assert_columns(dt, c("person_id", "sleep_date", "is_main_sleep",
                       "minute_in_bed", "minute_asleep"), "raw sleep")
  dt <- copy(dt)
  dt[, `:=`(date = as.Date(sleep_date), is_main = as_main_sleep_flag(is_main_sleep))]
  dt <- dt[!is.na(date)]

  minute_cols <- intersect(c("minute_in_bed", "minute_asleep", "minute_awake",
                             "minute_restless", "minute_deep", "minute_light",
                             "minute_rem", "minute_wake", "minute_after_wakeup"),
                           names(dt))

  # Sum main sleep + naps within the day. Efficiency is computed AFTER summing,
  # from the day totals (not as a mean of per-record ratios).
  dt[, is_nap := !is_main]
  flags <- dt[, .(n_records = .N, n_main = sum(is_main), n_nap = sum(is_nap)), by = .(person_id, date)]
  day <- merge(flags, grouped_sum_na(dt, c("person_id", "date"), minute_cols), by = c("person_id", "date"))
  day[, `:=`(has_main_sleep = n_main > 0L, has_nap = n_nap > 0L)]
  day[, c("n_main", "n_nap") := NULL]
  setcolorder(day, c("person_id", "date", "has_main_sleep", "has_nap", "n_records"))

  # Efficiency is only defined when both totals are positive. A day with time in
  # bed but zero minutes asleep is flagged (these are usually incomplete records)
  # and gets NA efficiency instead of a misleading 0.
  day[, asleep_zero_flag := minute_in_bed > 0 & minute_asleep == 0]
  day[, sleep_efficiency := fifelse(minute_in_bed > 0 & minute_asleep > 0,
                                    minute_asleep / minute_in_bed * 100, NA_real_)]
  n_bad <- day[sleep_efficiency > 100 | sleep_efficiency < 0, .N]
  if (n_bad > 0L) {
    log_info("sleep: %d day(s) with efficiency outside [0,100] set to NA", n_bad)
    day[sleep_efficiency > 100 | sleep_efficiency < 0, sleep_efficiency := NA_real_]
  }
  # the day-type flags are needed by the main-sleep-only sensitivity: never drop them
  drop_uninformative_columns(day, protect = c("person_id", "date", "has_main_sleep", "has_nap", "n_records"),
                             na_frac = na_frac)
}

# ---- heart rate ---------------------------------------------------------------

# "Out of Range" -> "out_of_range", "Fat Burn" -> "fat_burn", ...
zone_slug <- function(x) gsub("^_|_$", "", gsub("[^a-z0-9]+", "_", tolower(trimws(x))))

# raw columns: person_id, date, zone_name, min_heart_rate, max_heart_rate,
# minute_in_zone, calorie_count. See the header for the output columns.
#
# Wide columns created (<zone> in out_of_range, fat_burn, cardio, peak):
#   mins_<zone>, cal_<zone>, hrlo_<zone>, hrhi_<zone>   (zone minutes / calories / boundaries)
#   total_hr_min  = sum of the four zone minutes   (~ minutes with heart rate recorded)
#   elevated_min  = fat_burn + cardio + peak minutes
#   total_cal     = sum of the four zone calories
#
# A zone that is absent on a day which has other zone data is set to 0 (the zone
# had no minutes). The number of such fills is logged.
standardize_hr <- function(raw, zones = c("out_of_range", "fat_burn", "cardio", "peak"),
                           elevated = c("fat_burn", "cardio", "peak"),
                           max_min_per_day = 1440, na_frac = 0.90) {
  dt <- data.table::as.data.table(raw)
  assert_columns(dt, c("person_id", "date", "zone_name", "minute_in_zone",
                       "calorie_count"), "raw heart rate")
  dt <- copy(dt)
  dt[, `:=`(date = as.Date(date), zone = zone_slug(zone_name))]
  dt <- dt[!is.na(date)]

  unknown <- setdiff(unique(dt$zone), zones)
  if (length(unknown) > 0L) {
    log_info("hr: ignoring unexpected zone(s): %s", paste(unknown, collapse = ", "))
    dt <- dt[zone %in% zones]
  }
  if (nrow(dt) == 0L) {
    stop("heart-rate data has none of the expected zones (", paste(zones, collapse = ", "),
         "); check zone_name spelling", call. = FALSE)
  }
  for (col in c("min_heart_rate", "max_heart_rate")) {
    if (!col %in% names(dt)) dt[, (col) := NA_real_]
  }

  # one row per person-day-zone. The raw export already has that grain; a zone split across
  # shards (duplicates) is summed, with min/max for the boundary columns.
  if (anyDuplicated(dt, by = c("person_id", "date", "zone")) == 0L) {
    long <- dt[, .(person_id, date, zone, mins = minute_in_zone, cal = calorie_count,
                   hrlo = min_heart_rate, hrhi = max_heart_rate)]
  } else {
    key <- c("person_id", "date", "zone")
    long <- merge(grouped_sum_na(dt, key, c("minute_in_zone", "calorie_count")),
                  dt[, .(hrlo = suppressWarnings(min(min_heart_rate, na.rm = TRUE)),
                         hrhi = suppressWarnings(max(max_heart_rate, na.rm = TRUE))), by = key], by = key)
    long[is.infinite(hrlo), hrlo := NA_real_][is.infinite(hrhi), hrhi := NA_real_]
    setnames(long, c("minute_in_zone", "calorie_count"), c("mins", "cal"))
  }

  wide <- data.table::dcast(long, person_id + date ~ zone,
                            value.var = c("mins", "cal", "hrlo", "hrhi"))
  # make sure every expected column exists even if a zone never appears
  for (v in c("mins", "cal", "hrlo", "hrhi")) {
    for (z in zones) {
      nm <- paste0(v, "_", z)
      if (!nm %in% names(wide)) wide[, (nm) := NA_real_]
    }
  }

  # absent zone on a day that has other zone data -> 0
  fill_zero <- function(cols) {
    has_any <- rowSums(!is.na(wide[, ..cols])) > 0L
    n_fill <- 0L
    for (cl in cols) {
      idx <- which(has_any & is.na(wide[[cl]]))
      n_fill <- n_fill + length(idx)
      if (length(idx) > 0L) set(wide, idx, cl, 0)
    }
    n_fill
  }
  n_m <- fill_zero(paste0("mins_", zones))
  n_c <- fill_zero(paste0("cal_", zones))
  if (n_m + n_c > 0L) log_info("hr: filled %d absent zone minute/calorie value(s) with 0", n_m + n_c)

  mins_cols <- paste0("mins_", zones)
  wide[, total_hr_min := fifelse(rowSums(!is.na(.SD)) == 0L, NA_real_, rowSums(.SD, na.rm = TRUE)),
       .SDcols = mins_cols]
  wide[, elevated_min := fifelse(is.na(total_hr_min), NA_real_,
                                 rowSums(.SD, na.rm = TRUE)),
       .SDcols = paste0("mins_", elevated)]
  wide[, total_cal := fifelse(rowSums(!is.na(.SD)) == 0L, NA_real_, rowSums(.SD, na.rm = TRUE)),
       .SDcols = paste0("cal_", zones)]

  n_over <- wide[total_hr_min > max_min_per_day, .N]
  if (n_over > 0L) log_info("hr: %d day(s) have total zone minutes > %d", n_over, max_min_per_day)

  drop_uninformative_columns(wide, na_frac = na_frac)
}

# ---- dispatcher ---------------------------------------------------------------

# Build the day table for one stream from its raw data.
build_day_table <- function(raw, stream, cfg, apply_drop_rule = TRUE) {
  na_frac <- if (apply_drop_rule) cfg$etl$na_drop_frac else NA_real_
  dt <- switch(stream,
    activity = standardize_activity(raw, na_frac),
    sleep    = standardize_sleep(raw, na_frac),
    hr       = standardize_hr(raw, cfg$etl$hr_zones, cfg$etl$hr_elevated,
                              cfg$etl$hr_max_min_per_day, na_frac),
    stop("unknown stream: ", stream, call. = FALSE))
  assert_unique_keys(dt, c("person_id", "date"), paste(stream, "day table"))
  setkey(dt, person_id, date)
  dt
}

# Check (and report) that the upstream selection rule holds in the day tables.
# Only warns: the data are used as given, the check just documents the truth.
check_selection_rule <- function(days, stream, cfg) {
  s <- cfg$selection
  w <- days[, .(n_days = .N, first = min(date), last = max(date),
                max_gap = if (.N > 1L) max(as.integer(diff(date)) - 1L) else 0L),
            by = person_id]
  w[, pct_obs := 100 * n_days / (as.integer(last - first) + 1L)]
  res <- list(
    n_persons          = nrow(w),
    below_min_days     = sum(w$n_days < s$min_days),
    below_min_pct      = sum(w$pct_obs < s$min_pct_observed),
    above_max_gap      = sum(w$max_gap > s$max_gap_days),
    observed_max_gap   = max(w$max_gap)
  )
  if (res$below_min_days + res$below_min_pct + res$above_max_gap > 0L) {
    warning(sprintf(
      "%s: selection rule not fully met (n<%d days: %d; <%g%% observed: %d; gap>%d days: %d). Results are labelled as given.",
      stream, s$min_days, res$below_min_days, s$min_pct_observed, res$below_min_pct,
      s$max_gap_days, res$above_max_gap), call. = FALSE)
  }
  res
}
