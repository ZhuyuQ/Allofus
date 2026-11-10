# ==============================================================================
# stream_states.R
# ------------------------------------------------------------------------------
# Which streams (activity, heart rate, sleep) are present on each day, and what the combinations mean. Parts A-G, see ../STREAM.md.
# Descriptive only: estimates and their uncertainty, no tests. Shares of days carry person-clustered intervals (cluster_mean_ci); averages of
# person-level values use t intervals over people. Nothing here relates the previous day's steps to the states.
#
# Terms
#   window            a person's first to last record in ANY of the three streams; one row per person-day
#   stream not active the day lies before that stream's first record or after its last (removed and counted, not missing)
#   device on a day   the tracker record of the person on or before the day (as-of join; the latest record wins, ties -> alphabetically last)
#   removed days      stream not active > device without heart rate (the device of that day has no heart-rate sensor) > device unknown and the
#                     person's usual device does not measure heart rate (or there is none); in this order, so the counts add up
#   activity filter   an activity row is set to missing when sedentary minutes cover the whole day and there is no heart-rate record that day
#                     (flagged, not removed, when heart rate was recorded: the watch was worn but recorded no movement)
#   present           activity: a row remains after the filter; heart rate: a record with at least one minute; sleep: a MAIN sleep record dated that day
#                     (naps do not count; Fitbit dates sleep by the wake-up day, so it describes the night before)
#   state             1 all present | 2 all missing | 3 non-wear (activity only) | 4 night missing (activity + heart rate) | 5 other
#   worn day          state 1 or 4
# ==============================================================================

suppressPackageStartupMessages(library(data.table))

ST_STATE_LABELS <- c(`1` = "All present", `2` = "All missing", `3` = "Non-wear", `4` = "Night missing", `5` = "Other")
ST_COMBO_LABELS <- c(`111` = "activity + heart rate + sleep", `110` = "activity + heart rate (sleep missing)", `100` = "activity only", `000` = "all missing",
                     `101` = "activity + sleep (heart rate missing)", `011` = "heart rate + sleep (activity missing)", `010` = "heart rate only", `001` = "sleep only")
ST_COMBO_CATEGORY <- c(`111` = "expected", `110` = "expected", `100` = "expected", `000` = "expected",
                       `101` = "activity + sleep, heart rate missing", `011` = "other unexpected", `010` = "other unexpected", `001` = "other unexpected")
ST_REASONS <- c(stream_not_active = "stream not active", device_without_hr = "device without heart rate", device_unknown_removed = "device unknown, removed")

# ---- settings --------------------------------------------------------------------------------------------------------------------------

# no_hr_versions: device versions WITHOUT a heart-rate sensor (lower case, "fitbit " removed); product knowledge, NOT from the data: check the table
# stream_device_hr_lookup. A version that is not listed is assumed to measure heart rate.
st_config <- function(overrides = list()) {
  cfg <- list(seed = 20260511L,
              no_hr_versions = c("one", "zip", "ultra", "classic", "flex", "flex 2", "force", "alta", "charge", "inspire", "ace", "ace 2", "ace 3", "mobiletrack"),
              keep_unknown_if_no_usual = FALSE,      # device unknown and no usual device at all: removed (the rule needs a usual device that measures heart rate)
              drop_unknown_days = FALSE,             # the check: remove EVERY day with an unknown device
              min_worn_days_habit = 30L, min_days_compare = 3L, min_people = 20L,
              hr_cuts_present = c(10, 20), probably_worn_cuts = c(20, 22),
              unexpected_max_share = 0.01, minutes_exceed_max = 0.001,
              sample_people = 25L, day_plot_max = 365L)
  utils::modifyList(cfg, overrides)
}

# ---- calendar: length of a day ----------------------------------------------------------------------------------------------------------

# The n-th Sunday of a month.
st_nth_sunday <- function(year, month, n) {
  first <- as.Date(sprintf("%d-%02d-01", year, month))
  first + ((8L - data.table::wday(first)) %% 7L) + 7L * (n - 1L)
}

# Shortest and longest possible length of each date in minutes over the US time zones: a day has 1,380 minutes where clocks go forward (2nd Sunday of March),
# 1,500 where they go back (1st Sunday of November) and 1,440 elsewhere, in the zones without daylight saving time as well. The time zone of a person is not known.
st_day_minutes <- function(dates) {
  yr <- unique(data.table::year(dates))
  lk <- data.table(year = yr, spring = st_nth_sunday(yr, 3L, 2L), fall = st_nth_sunday(yr, 11L, 1L))
  y <- match(data.table::year(dates), lk$year)
  list(min = fifelse(dates == lk$spring[y], 1380L, 1440L), max = fifelse(dates == lk$fall[y], 1500L, 1440L))
}

# ---- devices -----------------------------------------------------------------------------------------------------------------------------

st_norm_version <- function(v) sub("^fitbit +", "", tolower(trimws(as.character(v))))

# Does a device version measure heart rate? NA stays NA.
st_version_has_hr <- function(v, cfg) fifelse(is.na(v), NA, !(st_norm_version(v) %in% cfg$no_hr_versions))

# Tracker rows (no scales, with a version) of the people of the analysis from the RAW device export (cfg: the pipeline config for normalise_device_rows).
# Clip-on and phone trackers are kept: they are not in the L1 watch timeline but are devices without heart rate. Returns list(dated, all).
st_trackers <- function(raw, persons, cfg) {
  d <- normalise_device_rows(raw, cfg)[person_id %in% persons$person_id & !is_scale & !is.na(device_version)]
  list(dated = unique(d[!is.na(device_date), .(person_id, device_date, device_version)]), all = d)
}

# The usual device of each person: the version with the most records (ties: the most recent record, then alphabetical).
st_usual_device <- function(trackers, cfg) {
  u <- trackers$all[, .(n = .N, last = if (all(is.na(device_date))) as.Date(NA) else max(device_date, na.rm = TRUE)), by = .(person_id, device_version)]
  setorderv(u, c("person_id", "n", "last", "device_version"), order = c(1L, -1L, -1L, 1L), na.last = TRUE)
  u <- u[, .SD[1L], by = person_id][, .(person_id, usual_version = device_version)]
  u[, usual_hr := st_version_has_hr(usual_version, cfg)][]
}

# The device versions seen, with the heart-rate classification used (to be checked by a person who knows the products).
st_device_lookup <- function(trackers, cfg) {
  t <- trackers$all[, .(n_rows = .N, n_people = uniqueN(person_id)), by = device_version][order(-n_people)]
  t[, `:=`(has_hr = st_version_has_hr(device_version, cfg), basis = fifelse(st_norm_version(device_version) %in% cfg$no_hr_versions, "listed as without heart rate", "assumed to measure heart rate"))]
  t[]
}

# ---- the person-day table ---------------------------------------------------------------------------------------------------------------

# The three day tables of the people (the columns this analysis needs).
st_load_streams <- function(store, persons) {
  keep <- function(t) t[person_id %in% persons$person_id]
  ac <- intersect(c("steps", "sedentary_minutes", "lightly_active_minutes", "fairly_active_minutes", "very_active_minutes", "mvpa_minutes"), store$columns("activity"))
  sc <- intersect(c("has_main_sleep", "has_nap"), store$columns("sleep"))
  if (!"has_main_sleep" %in% sc) stop("the sleep day table has no has_main_sleep column: run the ETL of this version first", call. = FALSE)
  list(act = keep(store$get("activity", ac)), hr = keep(store$get("hr", "total_hr_min")), slp = keep(store$get("sleep", sc)))
}

# First and last record of each stream per person, and the window over the three streams.
st_windows <- function(tabs) {
  rg <- function(t, nm) setnames(t[, .(f = min(date), l = max(date)), by = person_id], c("f", "l"), paste0(c("first_", "last_"), nm))
  w <- Reduce(function(a, b) merge(a, b, by = "person_id", all = TRUE), list(rg(tabs$act, "act"), rg(tabs$slp, "slp"), rg(tabs$hr, "hr")))
  w[, `:=`(first = pmin(first_act, first_slp, first_hr, na.rm = TRUE), last = pmax(last_act, last_slp, last_hr, na.rm = TRUE))]
  w[, n_days := as.integer(last - first) + 1L][]
}

# One row per person-day of every window, with the records of the three streams attached and the 'stream active' flags.
st_grid <- function(tabs, persons) {
  win <- st_windows(tabs)
  D <- data.table(person_id = rep(win$person_id, win$n_days), date = rep(win$first, win$n_days) + sequence(win$n_days) - 1L)
  attach <- function(t, flag) {
    cols <- setdiff(names(t), c("person_id", "date")); set(D, j = flag, value = FALSE)
    D[t, on = c("person_id", "date"), (c(flag, cols)) := c(list(TRUE), mget(paste0("i.", cols)))]
  }
  attach(tabs$act, "act_row"); attach(tabs$hr, "hr_row"); attach(tabs$slp, "slp_row")
  D[win, on = "person_id", `:=`(act_active = date >= i.first_act & date <= i.last_act, slp_active = date >= i.first_slp & date <= i.last_slp,
                                hr_active = date >= i.first_hr & date <= i.last_hr)]
  for (f in c("act_active", "slp_active", "hr_active")) set(D, which(is.na(D[[f]])), f, FALSE)
  setnames(D, "total_hr_min", "hr_min")
  p <- as.data.table(persons)[, .(person_id, cohort, age_group, device_family)]
  D[p, on = "person_id", `:=`(cohort = i.cohort, age_group = i.age_group, device_family = i.device_family)]
  D[]
}

# The device of each day (as-of join) and of each person (usual device), and what they measure.
st_attach_device <- function(D, trackers, usual, cfg) {
  tl <- trackers$dated; setorder(tl, person_id, device_date, device_version); tl <- tl[, .SD[.N], by = .(person_id, device_date)]
  set(D, j = "dev_version", value = tl[D, .(v = device_version), on = .(person_id, device_date = date), roll = TRUE]$v)
  set(D, j = "dev_hr", value = st_version_has_hr(D$dev_version, cfg))
  D[usual, on = "person_id", usual_hr := i.usual_hr]
  D[]
}

# Activity filter and presence flags (independent of the removal rules).
st_flags <- function(D, cfg) {
  dm <- st_day_minutes(D$date)
  oth <- intersect(c("lightly_active_minutes", "fairly_active_minutes", "very_active_minutes"), names(D))
  none <- if (length(oth)) Reduce(`&`, lapply(oth, function(cl) is.na(D[[cl]]) | D[[cl]] == 0)) else TRUE
  set(D, j = "hr_present", value = D$hr_row & !is.na(D$hr_min) & D$hr_min >= 1)
  set(D, j = "covers_day", value = D$act_row & !is.na(D$sedentary_minutes) & D$sedentary_minutes >= dm$min & none)
  set(D, j = "act_filtered", value = D$covers_day & !D$hr_present)
  set(D, j = "motionless_worn", value = D$covers_day & D$hr_present)
  set(D, j = "act_after", value = D$act_row & !D$act_filtered)
  set(D, j = "sleep_main", value = D$slp_row & D$has_main_sleep %in% TRUE)
  set(D, j = c("day_min", "day_max"), value = list(dm$min, dm$max))
  D[]
}

# Which days are removed, in the fixed order (stream not active > device without heart rate > device unknown and the usual device does not measure it).
st_assign_reason <- function(D, cfg) {
  unknown <- is.na(D$dev_version)
  keep_unknown <- !cfg$drop_unknown_days & (D$usual_hr %in% TRUE | (is.na(D$usual_hr) & cfg$keep_unknown_if_no_usual))
  set(D, j = "reason", value = fcase(!(D$act_active & D$slp_active & D$hr_active), "stream_not_active", !unknown & D$dev_hr %in% FALSE, "device_without_hr",
                                     unknown & !keep_unknown, "device_unknown_removed", default = "analysed"))
  set(D, j = "unknown_kept", value = unknown & D$reason == "analysed")
  st_set_states(D)
}

# Combination code ("111", ...) and state (1-5) of every day.
st_set_states <- function(D) {
  a <- D$act_after; h <- D$hr_present; s <- D$sleep_main
  set(D, j = "combo", value = paste0(as.integer(a), as.integer(h), as.integer(s)))
  set(D, j = "combo_before", value = paste0(as.integer(D$act_row), as.integer(h), as.integer(s)))
  set(D, j = "state", value = fcase(a & h & s, 1L, !a & !h & !s, 2L, a & !h & !s, 3L, a & h & !s, 4L, default = 5L))
  set(D, j = "dtype", value = fifelse(data.table::wday(D$date) %in% 2:6, "weekday", "weekend"))
  D[]
}

# Everything per person-day, built once; the removal rule is applied by st_assign_reason() (so the check without unknown-device days reuses it).
st_build_days <- function(store, persons, trackers, cfg, pipeline_cfg) {
  D <- st_grid(st_load_streams(store, persons), persons)
  D <- st_attach_device(D, trackers, st_usual_device(trackers, cfg), cfg)
  D <- st_flags(D, cfg)
  setorder(D, person_id, date)
  D[, cohort := factor(cohort, levels = pipeline_cfg$cohorts)]
  st_assign_reason(D, cfg)
}

# ---- statistics helpers -------------------------------------------------------------------------------------------------------------------

# Share of days (%) with a person-clustered 95% interval, by group. y: the name of a logical column.
st_share_by <- function(dt, y, by = "cohort") {
  col <- y
  dt[, { ci <- cluster_mean_ci(as.numeric(get(col)), person_id)
         .(n_days = .N, n_events = sum(get(col)), n_people = ci$n_persons, share_pct = 100 * ci$mean, low = 100 * ci$lo, high = 100 * ci$hi) }, by = by]
}

# The same for several indicator columns: specs = c(label = column).
st_share_long <- function(dt, specs, by = "cohort") {
  rbindlist(lapply(names(specs), function(lb) cbind(step = lb, st_share_by(dt, specs[[lb]], by))))
}

# Mean of person-level values with a t interval over people.
st_mean_ci <- function(x) {
  x <- x[is.finite(x)]; n <- length(x)
  se <- if (n > 1L) stats::sd(x) / sqrt(n) else NA_real_
  m <- if (n) mean(x) else NA_real_
  data.table(n_people = n, mean = m, low = m - stats::qt(0.975, max(n - 1L, 1L)) * se, high = m + stats::qt(0.975, max(n - 1L, 1L)) * se)
}

# Quantiles of a vector as a one-row table.
st_quant <- function(x, probs = c(0.05, 0.25, 0.5, 0.75, 0.95)) {
  q <- if (length(x)) stats::quantile(x, probs, names = FALSE) else rep(NA_real_, length(probs))
  as.data.table(setNames(as.list(q), paste0("p", probs * 100)))
}

# Cells with 1 to min_people - 1 people are not shown: the count becomes NA (people_label "<20") and the values NA.
st_suppress <- function(tab, value_cols, cfg) {
  if (!nrow(tab)) return(tab)
  tab <- copy(tab); small <- !is.na(tab$n_people) & tab$n_people > 0L & tab$n_people < cfg$min_people
  set(tab, j = "people_label", value = ifelse(small, paste0("<", cfg$min_people), as.character(tab$n_people)))
  for (v in intersect(value_cols, names(tab))) set(tab, which(small), v, NA)
  set(tab, which(small), "n_people", NA_integer_)
  tab[, suppressed := small][]
}

ST_VALUE_COLS <- c("n_days", "n_events", "share_pct", "low", "high", "mean", "mean_state3", "mean_state1", "relative_diff_pct", "p5", "p25", "p50", "p75", "p95", "n_nights")

# ---- Part A: data flow ------------------------------------------------------------------------------------------------------------------

# From all days in the windows to the days analysed, by cohort: days, share of the window days (clustered interval) and people at each step.
st_part_a <- function(D, cfg) {
  W <- D[, .(person_id, cohort, reason, unknown_kept, i_all = TRUE, i_stream = reason == "stream_not_active", i_act = !act_active, i_slp = !slp_active, i_hr = !hr_active,
             i_nohr = reason == "device_without_hr", i_unk = reason == "device_unknown_removed", i_ok = reason == "analysed")]
  W[, i_unk_kept := unknown_kept]
  specs <- c(`all days in the windows` = "i_all", `removed: a stream not active` = "i_stream", `  activity not active (not exclusive)` = "i_act",
             `  sleep not active (not exclusive)` = "i_slp", `  heart rate not active (not exclusive)` = "i_hr", `removed: device without heart rate` = "i_nohr",
             `removed: device unknown (usual device without heart rate, or none)` = "i_unk", `analysed days` = "i_ok", `  of which device unknown, kept (usual device measures heart rate)` = "i_unk_kept")
  flow <- st_share_long(W, specs)
  setnames(flow, c("n_days", "n_events"), c("n_window_days", "n_days_at_step"))            # n_days_at_step = days with the characteristic of the step; n_window_days = all days of the cohort
  flow <- st_suppress(flow, c(ST_VALUE_COLS, "n_window_days", "n_days_at_step"), cfg)
  people <- D[, .(n_people = uniqueN(person_id), n_people_with_analysed_days = uniqueN(person_id[reason == "analysed"]), window_days = .N,
                  analysed_days = sum(reason == "analysed")), by = cohort]
  list(flow = flow, people = people)
}

# ---- Part B: the activity filter ------------------------------------------------------------------------------------------------------------

# On the analysed days: activity rows set to missing, people affected, worn-but-motionless days, and the spread of the share across people.
st_part_b <- function(Dk, cfg) {
  A <- Dk[act_row == TRUE]
  tab <- st_share_long(A, c(`activity rows set to missing by the filter` = "act_filtered", `activity rows kept and flagged: worn but motionless` = "motionless_worn"))
  aff <- A[, .(n_activity_rows = .N, n_people_with_rows = uniqueN(person_id), n_people_affected = uniqueN(person_id[act_filtered]),
               n_flagged_motionless = sum(motionless_worn)), by = cohort]
  per <- A[, .(share = mean(act_filtered), n = .N), by = .(cohort, person_id)]
  q <- function(x) as.list(stats::setNames(100 * stats::quantile(x, c(0.5, 0.75, 0.9, 0.99), names = FALSE), c("p50_pct", "p75_pct", "p90_pct", "p99_pct")))
  spread <- per[, c(.(n_people = .N, mean_person_share_pct = 100 * mean(share), share_people_with_any_pct = 100 * mean(share > 0)), q(share)), by = cohort]
  list(share = st_suppress(tab, ST_VALUE_COLS, cfg), people = aff, spread = spread, person_share = per)
}

# ---- Part C: patterns and nesting ---------------------------------------------------------------------------------------------------------

# Share of days in each of the eight combinations before and after the filter, by cohort.
st_combos <- function(Dk, cfg) {
  one <- function(col, stage) rbindlist(lapply(names(ST_COMBO_LABELS), function(cb) {
    d <- Dk[, .(person_id, cohort, hit = get(col) == cb)]
    cbind(stage = stage, combo = cb, st_share_by(d, "hit"))
  }))
  t <- rbind(one("combo_before", "before the filter"), one("combo", "after the filter"))
  t[, `:=`(label = ST_COMBO_LABELS[combo], category = ST_COMBO_CATEGORY[combo])]
  st_suppress(t, ST_VALUE_COLS, cfg)
}

# Before the filter: how often heart rate is present when activity is missing, and sleep when heart rate is missing (both should be close to zero).
st_nesting_checks <- function(Dk, cfg) {
  a <- Dk[act_row == FALSE, .(person_id, cohort, y = hr_present)]; b <- Dk[hr_present == FALSE, .(person_id, cohort, y = sleep_main)]
  t <- rbind(cbind(check = "heart rate present when activity is missing (before the filter)", st_share_by(a, "y")),
             cbind(check = "sleep present when heart rate is missing (before the filter)", st_share_by(b, "y")))
  st_suppress(t, ST_VALUE_COLS, cfg)
}

# Share of analysed days (after the filter) in the unexpected combinations other than 'activity + sleep, heart rate missing', with the streams shifted
# by -1, 0, +1 days within the person: a convention mismatch (local date vs UTC) shows as a lower share at another lag.
st_date_check <- function(D, cfg) {
  X <- D[, .(person_id, reason, a = act_after, h = hr_present, s = sleep_main)]
  rows <- rbindlist(lapply(c("sleep", "heart rate"), function(st) rbindlist(lapply(-1:1, function(lag) {
    Y <- copy(X); col <- if (st == "sleep") "s" else "h"
    Y[, (col) := shift(get(col), -lag, type = "lag"), by = person_id]
    Y <- Y[reason == "analysed" & !is.na(get(col))]
    cb <- paste0(as.integer(Y$a), as.integer(Y$h), as.integer(Y$s))
    data.table(shifted_stream = st, lag_days = lag, n_days = nrow(Y), share_other_unexpected_pct = 100 * mean(cb %in% c("011", "010", "001")))
  }))))
  rows[, best_lag := lag_days[which.min(share_other_unexpected_pct)], by = shifted_stream][]
}

# Check of the minutes: share of heart-rate records whose minutes exceed the length of the day (1,440; 1,500 where clocks go back).
st_minutes_check <- function(D, cfg) {
  H <- D[hr_row == TRUE & !is.na(hr_min), .(person_id, cohort, exceed = hr_min > day_max)]
  list(by_cohort = st_suppress(st_share_by(H, "exceed"), ST_VALUE_COLS, cfg), overall = st_share_by(H[, .(person_id, cohort = "all", exceed)], "exceed"))
}

# The gate: parts D-G are only run when the date check and the minutes check pass.
st_gate <- function(D, cfg, minutes) {
  dc <- st_date_check(D, cfg)
  lag0 <- dc[lag_days == 0L]
  date_ok <- all(lag0$share_other_unexpected_pct <= 100 * cfg$unexpected_max_share) && all(dc$best_lag == 0L)
  mins <- minutes$overall$share_pct / 100
  data.table(check = c("other unexpected combinations <= 1% of days and best lag = 0", "heart-rate minutes above the length of the day <= 0.1%"),
             value_pct = c(max(lag0$share_other_unexpected_pct), 100 * mins), limit_pct = 100 * c(cfg$unexpected_max_share, cfg$minutes_exceed_max),
             passed = c(date_ok, mins <= cfg$minutes_exceed_max))
}

# ---- Part D: what state 3 is ---------------------------------------------------------------------------------------------------------------

# Steps on state-3 days (distribution) and, within the same type of day, the difference between a person's state-3 and state-1 days.
st_part_d <- function(Dk, cfg) {
  s3 <- Dk[state == 3L & !is.na(steps)]
  dist <- s3[, c(.(n_days = .N, n_people = uniqueN(person_id), share_zero_steps_pct = 100 * mean(steps == 0), mean_steps = mean(steps)), as.list(st_quant(steps))), by = cohort]
  list(distribution = st_suppress(dist, ST_VALUE_COLS, cfg), steps = s3[, .(cohort, steps)], compare = st_state3_compare(Dk, cfg))
}

# Person-level means by state and type of day, people with at least min_days_compare days of each state in that type; 'all' = the type-specific
# differences averaged with the weight of the person's state-3 days of each type.
st_person_state_means <- function(Dk) {
  d <- Dk[state %in% c(1L, 3L) & !is.na(steps)]
  m <- d[, .(n = .N, steps = mean(steps), mvpa = mean(mvpa_minutes, na.rm = TRUE)), by = .(person_id, cohort, dtype, state)]
  dcast(m, person_id + cohort + dtype ~ state, value.var = c("n", "steps", "mvpa"))
}

st_state3_compare <- function(Dk, cfg) {
  w <- st_person_state_means(Dk)
  w <- w[!is.na(n_1) & !is.na(n_3) & n_1 >= cfg$min_days_compare & n_3 >= cfg$min_days_compare]
  by_type <- rbindlist(lapply(c("steps", "mvpa"), function(mt) w[, {
    d <- get(paste0(mt, "_3")) - get(paste0(mt, "_1"))
    c(list(metric = mt), st_mean_ci(d), list(mean_state3 = mean(get(paste0(mt, "_3"))), mean_state1 = mean(get(paste0(mt, "_1")))))
  }, by = .(cohort, day_type = dtype)]))
  all_t <- rbindlist(lapply(c("steps", "mvpa"), function(mt) {
    d <- w[, .(diff = sum(n_3 * (get(paste0(mt, "_3")) - get(paste0(mt, "_1")))) / sum(n_3), m3 = weighted.mean(get(paste0(mt, "_3")), n_3), m1 = weighted.mean(get(paste0(mt, "_1")), n_3)), by = .(person_id, cohort)]
    d[, c(list(metric = mt, day_type = "all (weighted by state-3 days)"), st_mean_ci(diff), list(mean_state3 = mean(m3), mean_state1 = mean(m1))), by = cohort]
  }))
  out <- rbind(by_type, all_t, use.names = TRUE)
  out[, relative_diff_pct := 100 * mean / mean_state1]
  st_suppress(out, ST_VALUE_COLS, cfg)
}

# ---- Part E: composition ---------------------------------------------------------------------------------------------------------------------

# Share of days in each state, pooled over days (long follow-up weighs more) and as the average of each person's own share (everybody weighs the same).
st_comp_one <- function(Dk, level_col, cfg) {
  d <- Dk[, .(person_id, cohort, state, level = if (level_col == "all") "all" else as.character(get(level_col)))][!is.na(level)]
  pooled <- rbindlist(lapply(1:5, function(k) { d[, y := state == k]; cbind(state = k, st_share_by(d, "y", c("cohort", "level"))) }))
  per <- d[, .(n = .N, k1 = mean(state == 1L), k2 = mean(state == 2L), k3 = mean(state == 3L), k4 = mean(state == 4L), k5 = mean(state == 5L)), by = .(cohort, level, person_id)]
  avg <- rbindlist(lapply(1:5, function(k) cbind(state = k, per[, st_mean_ci(100 * get(paste0("k", k))), by = .(cohort, level)])))
  setnames(avg, c("mean", "low", "high"), c("share_pct", "low", "high"))
  pooled[, weighting := "all days (long follow-up weighs more)"]; avg[, weighting := "person average (every person weighs the same)"]
  out <- rbind(pooled, avg, fill = TRUE)
  out[, `:=`(stratum = level_col, state_label = ST_STATE_LABELS[as.character(state)])]
  st_suppress(out, ST_VALUE_COLS, cfg)
}

st_composition <- function(Dk, cfg) rbindlist(lapply(c("all", "age_group", "device_family"), function(v) st_comp_one(Dk, v, cfg)), fill = TRUE)

# Habitual night non-wear: the share of a person's worn days (states 1 and 4) that are state 4, people with enough worn days (description only: it uses the whole
# follow-up).
st_night_habit <- function(Dk, cfg) {
  p <- Dk[state %in% c(1L, 4L), .(worn_days = .N, share = mean(state == 4L)), by = .(cohort, person_id)][worn_days >= cfg$min_worn_days_habit]
  q <- function(x) as.list(stats::setNames(100 * stats::quantile(x, c(0.25, 0.5, 0.75), names = FALSE), c("p25", "p50", "p75")))
  tab <- p[, c(st_mean_ci(100 * share), q(share)), by = cohort]
  list(table = st_suppress(tab, ST_VALUE_COLS, cfg), persons = p)
}

# A seeded sample of people per cohort for the day-by-day picture (state 0 = day removed); never written to a table.
st_day_sample <- function(D, cfg) {
  ids <- D[, .(id = unique(person_id)), by = cohort][, .SD[with_seed(cfg$seed + as.integer(cohort[1L]), sample.int(.N, min(.N, cfg$sample_people)))], by = cohort]
  s <- D[person_id %in% ids$id, .(cohort, person_id, date, state = fifelse(reason == "analysed", state, 0L))]
  s[, day := as.integer(date - min(date)) + 1L, by = person_id]
  s[day <= cfg$day_plot_max][, row := .GRP, by = person_id][]
}

# ---- Part F: what 'night missing' contains ------------------------------------------------------------------------------------------------

# Nights (state 4 and, for comparison, state 1) whose previous day was also a worn day, with the heart-rate hours of both days.
st_nights <- function(D) {
  X <- D[, .(person_id, cohort, state = fifelse(reason == "analysed", state, 0L), hr_min)]
  X[, `:=`(prev_state = shift(state), prev_hr = shift(hr_min)), by = person_id]
  X[state %in% c(1L, 4L) & prev_state %in% c(1L, 4L) & !is.na(hr_min) & !is.na(prev_hr),
    .(person_id, cohort, state, hr_prev = prev_hr / 60, hr_now = hr_min / 60)]
}

# Share of nights where heart-rate wear was high on both days ('probably worn') at each cut-off, state 4 against state 1.
st_part_f <- function(D, cfg) {
  N <- st_nights(D)
  tab <- rbindlist(lapply(cfg$probably_worn_cuts, function(ct) {
    N[, hit := pmin(hr_prev, hr_now) >= ct]
    cbind(cutoff_hours = ct, st_share_by(N[, state_label := ST_STATE_LABELS[as.character(state)]], "hit", c("cohort", "state_label")))
  }))
  setnames(tab, "n_days", "n_nights")
  list(table = st_suppress(tab, ST_VALUE_COLS, cfg), nights = N)
}

# ---- Part G: partial wear on present days --------------------------------------------------------------------------------------------------

# Heart-rate wear on worn days (states 1 and 4): distribution and the share below each cut-off.
st_part_g <- function(Dk, cfg) {
  W <- Dk[state %in% c(1L, 4L) & !is.na(hr_min), .(person_id, cohort, hrh = hr_min / 60)]
  dist <- W[, c(.(n_days = .N, n_people = uniqueN(person_id)), as.list(st_quant(hrh, c(0.05, 0.1, 0.25, 0.5, 0.75, 0.95)))), by = cohort]
  below <- rbindlist(lapply(cfg$hr_cuts_present, function(ct) { W[, hit := hrh < ct]; cbind(cutoff_hours = ct, st_share_by(W, "hit")) }))
  list(distribution = st_suppress(dist, ST_VALUE_COLS, cfg), below = st_suppress(below, ST_VALUE_COLS, cfg), hours = W[, .(cohort, hrh)])
}

# ---- run ---------------------------------------------------------------------------------------------------------------------------------------

# Parts A-C and the checks (the gate): list of tables and plot data.
st_run_abc <- function(D, cfg) {
  Dk <- D[reason == "analysed"]
  a <- st_part_a(D, cfg); b <- st_part_b(Dk, cfg); mc <- st_minutes_check(D, cfg)
  list(flow = a$flow, people = a$people, activity_filter = b$share, activity_filter_people = b$people, activity_filter_spread = b$spread, person_share = b$person_share,
       combos = st_combos(Dk, cfg), nesting_checks = st_nesting_checks(Dk, cfg), date_check = st_date_check(D, cfg), minutes_check = mc$by_cohort,
       gate = st_gate(D, cfg, mc))
}

# Parts D-G.
st_run_dg <- function(D, cfg) {
  Dk <- D[reason == "analysed"]
  d <- st_part_d(Dk, cfg); nh <- st_night_habit(Dk, cfg); f <- st_part_f(D, cfg); g <- st_part_g(Dk, cfg)
  list(steps_state3 = d$distribution, steps_data = d$steps, state3_vs_state1 = d$compare, composition = st_composition(Dk, cfg), night_habit = nh$table,
       night_habit_persons = nh$persons, day_sample = st_day_sample(D, cfg), probably_worn = f$table, nights = f$nights, hr_wear = g$distribution, hr_wear_below = g$below, hr_hours = g$hours)
}

# ---- writing ------------------------------------------------------------------------------------------------------------------------------------

# Tables that go to <outdir>/tables as stream_<name>.tsv. Person-level and day-level data that only feed figures are never written.
ST_PLOT_ONLY <- c("person_share", "steps_data", "night_habit_persons", "day_sample", "nights", "hr_hours")

st_write_tables <- function(tabs, outdir) {
  dir.create(file.path(outdir, "tables"), recursive = TRUE, showWarnings = FALSE)
  keep <- setdiff(names(tabs), ST_PLOT_ONLY)
  for (nm in keep) if (is.data.frame(tabs[[nm]]) && nrow(tabs[[nm]])) fwrite(tabs[[nm]], file.path(outdir, "tables", paste0("stream_", nm, ".tsv")), sep = "\t")
  invisible(keep)
}

# Save figures as PNG (<outdir>/figures/<name>.png).
st_save_figures <- function(figs, outdir, dpi = 150) {
  dir.create(file.path(outdir, "figures"), recursive = TRUE, showWarnings = FALSE)
  size <- function(nm) switch(sub("^stream_([a-g]).*", "\\1", nm), a = c(9, 5), b = c(8, 6), c = c(11, 8), d = c(9, 6), e = c(11, 8), f = c(11, 8), g = c(8, 5), c(9, 6))
  for (nm in names(figs)) { sz <- size(nm); ggplot2::ggsave(file.path(outdir, "figures", paste0(nm, ".png")), figs[[nm]], width = sz[1], height = sz[2], dpi = dpi, bg = "white") }
  invisible(names(figs))
}

# ---- gate, one variant, methods text --------------------------------------------------------------------------------------------------------

# Stops when the date check (Part C) or the minutes check (Part G) failed: the remaining parts must not run on data that are not prepared correctly.
st_require_gate <- function(gate) {
  if (all(gate$passed)) return(invisible(TRUE))
  stop("the checks failed: ", paste(gate$check[!gate$passed], collapse = "; "), ". Fix the data preparation (date convention of the three tables, heart-rate zone minutes) ",
       "before running parts D-G; see stream_date_check and stream_minutes_check.", call. = FALSE)
}

# Parts A-C with the gate, then D-G, tables and figures written to outdir (used for the check without unknown-device days; the main run uses the cells of pipeline/run_stream.R).
st_run_variant <- function(D, cfg, pcfg, outdir) {
  abc <- st_run_abc(D, cfg); st_write_tables(abc, outdir); st_save_figures(st_figures(abc, NULL, pcfg, cfg), outdir)
  st_require_gate(abc$gate)
  dg <- st_run_dg(D, cfg); st_write_tables(dg, outdir); st_save_figures(st_figures(abc, dg, pcfg, cfg), outdir)
  invisible(list(abc = abc, dg = dg))
}

# The Methods text (also in STREAM.md), written next to the results.
st_methods_text <- function(cfg) c(
  "# Methods: stream states", "",
  "**People.** The people of the existing selection (no re-selection). **Days.** One row per person-day; a person's window runs from the first to the last record in any of the three streams.", "",
  "**Removed and counted (not treated as missing).** (1) Days when a stream was not active: the day falls before that stream's first record or after its last. (2) Days on a device without a heart-rate sensor, decided by the device used on that day (the tracker record of the person on or before the day, from the raw device export; the latest record wins). When the device of a day is unknown the day is kept if the person's usual device (the version with the most records) measures heart rate, otherwise it is removed; a person without any device record has no usual device, so such days are removed. A check repeats the analysis without ALL days of unknown device. Which versions have no heart-rate sensor is product knowledge (`no_hr_versions`, table stream_device_hr_lookup), not data. A day is removed for the first of these reasons, in the order stream not active, device without heart rate, device unknown.", "",
  "**Activity filter.** An activity row is set to missing when sedentary minutes cover the whole day (1,440 minutes; 1,380 where clocks go forward; the day length differs by time zone and is not known, so the shortest possible length of the date is used, and all other active minutes must be zero) and there is no heart-rate record that day. If sedentary minutes cover the whole day but heart rate was recorded, the row is kept and flagged (the watch was worn but recorded no movement).", "",
  "**Presence.** Activity: a row remains after the filter. Heart rate: a record with at least one minute. Sleep: a MAIN sleep record dated that day (naps do not count); Fitbit dates sleep by the wake-up day, so it describes the night before. **States:** 1 all present, 2 all missing, 3 non-wear (activity present, heart rate and sleep missing), 4 night missing (activity and heart rate present, sleep missing), 5 other. Worn days are states 1 and 4. Habitual night non-wear = the share of a person's worn days that are state 4, for people with at least 30 worn days, for description only (it is computed from the whole follow-up).", "",
  "**Heart rate switched off.** There is no separate rule for heart rate being switched off. The selection already required the longest heart-rate gap to be at most 14 days, so long switched-off periods were removed before this analysis.", "",
  "**Uncertainty.** Shares of days carry 95% intervals that account for repeated days of the same person (ratio estimator with people as clusters); averages of person-level values carry t intervals over people. No significance tests. Cells with 1-19 people are not shown. Dates are used as stored (local date versus UTC is checked in Part C). Cohort differences are unadjusted.", "",
  "**Gate.** Parts D-G run only if the unexpected combinations other than 'activity + sleep, heart rate missing' are at most 1% of days with the best date shift at 0, and at most 0.1% of heart-rate records have minutes above the length of the day. The steps of the previous day are never related to the states here.")
