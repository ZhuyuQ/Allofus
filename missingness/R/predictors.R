# ==============================================================================
# predictors.R
# ------------------------------------------------------------------------------
# Build the prior-day predictors for the gap-onset models (and the gap-level models B and C).
#
# For every predictor p (e.g. log steps) we create, per observed day d of a person:
#   x_p   the transformed value
#   m_p   the person's USUAL level on day d (the baseline, see below)
#   w_p   WITHIN-person deviation:  (x - m) / sd_w      "how much higher than THIS person's usual was day d"
#         (+1 = one within-person SD above usual; the primary predictor)
#   b_p   BETWEEN-person level:     (m - mu_b) / sd_b   "is this person generally higher or lower than other people"
#   tw_p  like w_p but with the TRAILING baseline (mean of the previous 14 observed days; sensitivity)
# so that x = m + sd_w * w and m = mu_b + sd_b * b: the two parts add up to the day's value and are separate effects.
#
# standardization = "pre_gap" (primary, cfg$standardization$method)
#   m = mean of the person's observed days STRICTLY BEFORE d (expanding history; needs >= baseline_min_n of them, else NA).
#   Only data observed before day d - hence before any gap that starts on d + 1 - enter the person's usual, so the
#   within- and between-person z-scores carry no information from the future of that person.
#   sd_w = root mean square of (x - m) over all days that have a baseline ("pooled within-person SD about the prior usual");
#   mu_b, sd_b = mean and SD of m over those same days (the reference distribution of the between-person level).
#   sd_w, mu_b, sd_b are one constant per predictor and stream: they only rescale, they never move a person's own
#   values, and they are stored in attr(result, "scales").
# standardization = "full_period" (sensitivity only)
#   m = mean over ALL of the person's observed days (includes days after d: look-ahead); sd_w = pooled within-person SD,
#   mu_b, sd_b = mean and SD of the person means over persons. This was the specification before the pre-gap change.
# ==============================================================================

# Mean of the observed values BEFORE each position (the position itself excluded; NA values ignored), NA until at least
# `min_n` observed values are available. Cumulative sums: linear time, no rolling-window function.
prior_mean <- function(x, min_n) {
  ok <- !is.na(x)
  cs <- cumsum(ifelse(ok, x, 0)); cn <- cumsum(ok)
  s_prev <- c(0, cs[-length(cs)]); n_prev <- c(0L, cn[-length(cn)])
  out <- s_prev / n_prev
  out[n_prev < min_n] <- NA_real_
  out
}

# Same, but the position itself is included (the person's usual as known at the END of that day).
through_mean <- function(x, min_n) {
  ok <- !is.na(x)
  out <- cumsum(ifelse(ok, x, 0)) / cumsum(ok)
  out[cumsum(ok) < min_n] <- NA_real_
  out
}

# Mean of the previous `n` observations (today excluded; NA values ignored), NA until at
# least `min_n` non-missing values are available. Cumulative sums, so no rolling-window
# function (whose NA/adaptive options differ between data.table versions) is needed.
trailing_mean_excl_current <- function(x, n, min_n) {
  k  <- seq_along(x)
  S  <- c(0, cumsum(ifelse(is.na(x), 0, x)))        # S[j + 1] = sum of the first j values
  C  <- c(0, cumsum(!is.na(x)))                      # C[j + 1] = count of non-NA among the first j
  hi <- k - 1L                                       # last position in the window (today excluded)
  lo <- pmax(0L, hi - n)                             # positions before the window
  cnt <- C[hi + 1L] - C[lo + 1L]
  out <- (S[hi + 1L] - S[lo + 1L]) / cnt
  out[hi < 1L | cnt < min_n] <- NA_real_
  out
}

# `days`       : day table of ONE stream (keyed person_id, date)
# `pred_names` : predictor names (from cfg$predictors) that belong to this stream
# trailing = FALSE skips the trailing-baseline columns (tw_); keep_x = FALSE drops the transformed
# raw value (x_) once w_ and b_ are made. Both save memory; the models only need w_ and b_.
# standardization : "pre_gap" (primary) or "full_period" (sensitivity); see the header of this file.
# keep_baseline = TRUE (pre_gap only) also keeps m_<p> = the person's usual THROUGH day d (day d included): the
#   baseline that is known when a gap starts on d + 1; the gap-level models standardize the days AFTER a gap with it.
predictor_table <- function(days, pred_names, cfg, trailing = TRUE, keep_x = TRUE,
                            standardization = cfg$standardization$method %||% "pre_gap", keep_baseline = FALSE) {
  stopifnot(standardization %in% c("pre_gap", "full_period"))
  min_n <- cfg$standardization$baseline_min_n %||% 7L
  dt <- days[, c("person_id", "date"), with = FALSE]
  scales <- vector("list", length(pred_names))

  for (i in seq_along(pred_names)) {
    p    <- pred_names[i]
    spec <- cfg$predictors[[p]]
    assert_columns(days, spec$col, paste("day table for predictor", p))
    x_nm <- paste0("x_", p)
    dt[, (x_nm) := apply_transform(days[[spec$col]], spec$transform)]

    if (standardization == "pre_gap") {
      # usual = mean of the earlier observed days of the person; constants from the days that have such a baseline
      dt[, m_i := prior_mean(get(x_nm), min_n), by = person_id]
      dev  <- dt[[x_nm]] - dt$m_i
      sd_w <- sqrt(mean(dev^2, na.rm = TRUE))
      mu_b <- mean(dt$m_i, na.rm = TRUE)
      sd_b <- stats::sd(dt$m_i, na.rm = TRUE)
      if (keep_baseline) dt[, (paste0("m_", p)) := through_mean(get(x_nm), min_n), by = person_id]
    } else {
      # usual = mean over the person's whole observation period (look-ahead); constants over days / persons
      dt[, m_i := mean(get(x_nm), na.rm = TRUE), by = person_id]
      dt[is.nan(m_i), m_i := NA_real_]
      dev  <- dt[[x_nm]] - dt$m_i
      n_ok <- sum(!is.na(dev))
      n_p  <- dt[!is.na(m_i), uniqueN(person_id)]
      sd_w <- sqrt(sum(dev^2, na.rm = TRUE) / max(n_ok - n_p, 1L))
      pm   <- unique(dt[!is.na(m_i), .(person_id, m_i)])
      sd_b <- stats::sd(pm$m_i)
      mu_b <- mean(pm$m_i)
    }

    dt[, (paste0("w_", p)) := (get(x_nm) - m_i) / sd_w]
    dt[, (paste0("b_", p)) := (m_i - mu_b) / sd_b]

    # trailing baseline: mean of the previous `trailing_n` observed days (excluding
    # today), needing at least `trailing_min_n` of them
    if (trailing) {
      dt[, tb := trailing_mean_excl_current(get(x_nm), cfg$onset$trailing_n, cfg$onset$trailing_min_n),
         by = person_id]
      dt[, (paste0("tw_", p)) := (get(x_nm) - tb) / sd_w]
      dt[, tb := NULL]
    }
    dt[, m_i := NULL]
    if (!keep_x) dt[, (x_nm) := NULL]
    scales[[i]] <- data.table(predictor = p, sd_within = sd_w, sd_between = sd_b,
                              mean_of_person_means = mu_b, standardization = standardization)
  }
  setattr(dt, "scales", rbindlist(scales))
  dt[]
}

# Evaluate a partial-wear rule (an R expression as text, see cfg$onset$partial_wear)
# on a day table; returns a logical vector (NA treated as FALSE).
partial_wear_flag <- function(days, rule) {
  out <- tryCatch(eval(str2lang(rule), envir = days), error = function(e)
    stop(sprintf("partial-wear rule '%s' cannot be evaluated on this day table: %s", rule, conditionMessage(e)), call. = FALSE))
  out[is.na(out)] <- FALSE
  out
}

# ---- one self-contained bundle per stream -----------------------------------------

# Everything the onset and gap-level models need for one stream, computed ONCE:
#   gaps   : one row per gap
#   rs     : the risk set (one row per observed day except the last) with outcomes,
#            history, calendar variables, own predictors, activity predictors,
#            partial-wear flags, person covariates and the person's missing days in the run-in
#   scales : within/between SDs used for scaling
#
# `observed_filter` (optional, e.g. "has_main_sleep") restricts which sleep/hr/activity
# rows count as "observed" (used for the main-sleep-only sensitivity).
# trailing = TRUE also builds the trailing-baseline columns (only the sensitivity analyses use them).
# standardization: "pre_gap" (primary) or "full_period" (the full-period sensitivity), see the header of this file.
prepare_stream_bundle <- function(stream, day_tables, persons, cfg, observed_filter = NULL, trailing = FALSE,
                                  standardization = cfg$standardization$method %||% "pre_gap") {
  days_full <- day_tables[[stream]]
  days <- if (is.null(observed_filter)) days_full
          else days_full[eval(str2lang(observed_filter), envir = days_full) %in% TRUE]

  gaps <- get_gap_runs(days)
  rs   <- build_risk_set(days, cfg)

  # Columns are added to `rs` BY REFERENCE (update join), so the big table is never copied.
  add_columns <- function(tab, key) {
    cols <- setdiff(names(tab), key)
    rs[tab, on = key, (cols) := mget(paste0("i.", cols))]
    invisible(NULL)
  }

  # own-stream predictors (computed on the stream's observed days)
  own_names <- cfg$own_predictors[[stream]]
  own_tab   <- predictor_table(days, own_names, cfg, trailing = trailing, keep_x = FALSE, standardization = standardization)
  scales    <- attr(own_tab, "scales")
  own_tab[, pw_own := partial_wear_flag(days, cfg$onset$partial_wear[[stream]])]
  add_columns(own_tab, c("person_id", "date"))

  # activity predictors for the cross-stream models
  if (stream != "activity") {
    act_days <- day_tables[["activity"]]
    act_tab  <- predictor_table(act_days, cfg$activity_predictors, cfg, trailing = trailing, keep_x = FALSE, standardization = standardization)
    act_tab[, pw_act := partial_wear_flag(act_days, cfg$onset$partial_wear$activity)]
    add_columns(act_tab, c("person_id", "date"))
    scales <- rbind(scales, attr(act_tab, "scales"))
    rm(act_tab, act_days)
  } else {
    rs[, pw_act := pw_own]             # activity: own stream IS the activity stream
  }
  rm(own_tab)

  # person covariates
  cov <- person_covariates(persons, cfg)                 # cohort, age10, sex2, CCI and device family (covariates.R)
  add_columns(cov, "person_id")
  add_columns(runin_person_table(days, cfg$runin$days)[, .(person_id, runin_missing)], "person_id")      # missing days in the first days of follow-up (runin.R); only the run-in models use it
  setkey(rs, person_id, date)

  list(stream = stream, gaps = gaps, rs = rs, scales = scales,
       n_persons = uniqueN(rs$person_id), trailing = trailing, standardization = standardization)
}
