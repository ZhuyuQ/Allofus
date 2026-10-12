# ==============================================================================
# gaps.R
# ------------------------------------------------------------------------------
# The core of the pipeline. Everything is derived from the SORTED OBSERVED DATES
# of each person - no calendar grid is ever built, so cost scales with the number
# of observed rows (about 10 million), not with follow-up days.
#
# Definitions (window (a): first to last observed record, per person and stream)
#   observed day : a row exists in the day table
#   gap          : a run of >= 1 consecutive calendar days with no row, lying
#                  strictly between two observed days (so leading/trailing
#                  missingness is not visible by design)
#   gap length   : number of missing days in the run
#
# Worked example (O = observed, M = missing):   day  1 2 3 4 5 6 7 8 9 10
#                                               --   O O O M M O O M O O
#   gaps: days 4-5 (length 2) and day 8 (length 1); window = 10 days, 7 observed.
# ==============================================================================

# ---- gap runs ------------------------------------------------------------------

# One row per gap. `days` must have person_id and date (one row per person-day).
# Columns: person_id, gap_index, prev_obs_date (last observed day BEFORE the gap),
#          gap_start, gap_end, next_obs_date (first observed day AFTER the gap),
#          gap_length.
get_gap_runs <- function(days) {
  d <- days[, .(person_id, date)]
  setorder(d, person_id, date)
  d[, next_date := shift(date, type = "lead"), by = person_id]
  g <- d[!is.na(next_date) & as.integer(next_date - date) > 1L]
  g[, `:=`(prev_obs_date = date, next_obs_date = next_date)]
  g[, `:=`(gap_start  = prev_obs_date + 1L,
           gap_end    = next_obs_date - 1L,
           gap_length = as.integer(next_obs_date - prev_obs_date) - 1L)]
  g[, gap_index := seq_len(.N), by = person_id]
  g[, .(person_id, gap_index, prev_obs_date, gap_start, gap_end, next_obs_date, gap_length)]
}

# One row per MISSING DAY (only the missing days are enumerated, never the grid).
expand_gap_days <- function(gaps) {
  idx <- rep.int(seq_len(nrow(gaps)), gaps$gap_length)
  out <- gaps[idx, .(person_id, gap_index, gap_length)]
  out[, date := gaps$gap_start[idx] + (sequence(gaps$gap_length) - 1L)]
  out[]
}

# ---- person-level summaries ------------------------------------------------------

# first/last observed date, observed days and window length per person.
person_windows <- function(days) {
  w <- days[, .(first_date = min(date), last_date = max(date), n_observed = .N),
            by = person_id]
  w[, window_days := as.integer(last_date - first_date) + 1L]
  w[]
}

# Person-level missingness summary (one row per person), the unit of all cohort
# comparisons. Persons with no gap get n_gaps = 0 and gap_max = 0.
summarise_person_gaps <- function(days, gaps) {
  w <- person_windows(days)
  gs <- gaps[, .(n_gaps = .N, n_missing_gaps = sum(gap_length),
                 gap_mean = mean(gap_length), gap_median = as.numeric(median(gap_length)),
                 gap_max = max(gap_length),
                 n_gaps_ge2 = sum(gap_length >= 2L), n_gaps_ge7 = sum(gap_length >= 7L)),
             by = person_id]
  out <- merge(w, gs, by = "person_id", all.x = TRUE)
  zero_cols <- c("n_gaps", "n_missing_gaps", "gap_max", "n_gaps_ge2", "n_gaps_ge7")
  for (cl in zero_cols) set(out, which(is.na(out[[cl]])), cl, 0)
  out[, `:=`(n_missing = window_days - n_observed,
             pct_missing = 100 * (window_days - n_observed) / window_days)]

  # Internal consistency: the missing days implied by the window must equal the
  # sum of the gap lengths. A mismatch means duplicated or unsorted dates.
  if (any(out$n_missing != out$n_missing_gaps)) {
    stop("inconsistent gap accounting: window - observed != sum(gap lengths)", call. = FALSE)
  }

  # longest run of consecutive observed days
  d <- days[, .(person_id, date)]
  setorder(d, person_id, date)
  d[, seg := cumsum(c(FALSE, diff(as.integer(date)) > 1L)), by = person_id]
  streak <- d[, .N, by = .(person_id, seg)][, .(longest_streak = max(N)), by = person_id]
  merge(out, streak, by = "person_id", all.x = TRUE)[]
}

# ---- the risk set for gap ONSET ---------------------------------------------------

# One row per OBSERVED day that is not the person's last day ("prior day" d).
# The row asks: "does a gap start on day d+1?"   (see the worked example above:
# days 1,2,3,6,7,9 are rows; days 4,5,8 have no data and day 10 is the last day).
#
# Columns
#   date            prior day d (the day whose activity is the predictor)
#   gap_len_next    0 if day d+1 is observed, otherwise the length of the gap that starts there
#   y_ge<k>         1 if a gap of at least k days starts on d+1   (k in cfg$onset)
#   prev_gap_len    length of the most recent gap that ended on or before d (0 if none)
#   any_prior_gap   TRUE if the person has had a gap before d
#   days_since_return  days since the person came back from that gap (capped; 0 if none)
#   study_day       days since the person's first record, at d (1 = first record)
#   dow_t, year_t   weekday (Sun..Sat) and calendar year of the predicted day
#
# History columns use only information known at day d (model A).
build_risk_set <- function(days, cfg) {
  ks  <- sort(unique(c(cfg$onset$primary_k, cfg$onset$secondary_k)))
  cap <- cfg$onset$days_since_cap

  rs <- days[, .(person_id, date)]
  setorder(rs, person_id, date)
  # one assignment per statement: newer data.table versions reject mixing per-day and
  # per-person lengths inside a single grouped `:=`
  rs[, next_date  := shift(date, type = "lead"), by = person_id]
  rs[, prev_date  := shift(date, type = "lag"),  by = person_id]
  rs[, first_date := date[1L],                   by = person_id]

  # gap that ENDED right before this observed day (0 if none)
  rs[, gap_before := fifelse(is.na(prev_date), 0L, as.integer(date - prev_date) - 1L)]
  # carry the last completed gap forward within the person
  rs[, `:=`(prev_gap_len = nafill(fifelse(gap_before > 0L, gap_before, NA_integer_), type = "locf"),
            return_day   = nafill(fifelse(gap_before > 0L, as.integer(date), NA_integer_), type = "locf")),
     by = person_id]
  rs[, any_prior_gap := !is.na(prev_gap_len)]
  rs[, days_since_return := fifelse(any_prior_gap,
                                    pmin(as.integer(date) - return_day, cap), 0L)]
  rs[is.na(prev_gap_len), prev_gap_len := 0L]

  rs <- rs[!is.na(next_date)]            # the last observed day has no "next day"
  rs[, gap_len_next := as.integer(next_date - date) - 1L]
  for (k in ks) set(rs, j = paste0("y_ge", k), value = as.integer(rs$gap_len_next >= k))

  rs[, `:=`(t_date    = date + 1L,
            study_day = as.integer(date - first_date) + 1L)]
  rs[, `:=`(dow_t  = factor(wday(t_date), levels = 1:7,
                            labels = c("Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat")),
            year_t = year(t_date))]
  rs[, c("next_date", "prev_date", "first_date", "gap_before", "return_day", "t_date") := NULL]
  setkey(rs, person_id, date)
  rs[]
}

# Which rows lie up to and including each person's FIRST gap onset? (logical vector, no copy)
# Used by the first-gap-only sensitivity: removes any dependence on earlier gaps; at most one
# event per person.
first_gap_rows <- function(rs) {
  fo <- rs[gap_len_next > 0L, .(first_onset = min(date)), by = person_id]
  first_onset <- fo$first_onset[match(rs$person_id, fo$person_id)]
  is.na(first_onset) | rs$date <= first_onset
}

restrict_first_gap <- function(rs) rs[first_gap_rows(rs)]
