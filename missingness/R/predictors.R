# ==============================================================================
# predictors.R
# ------------------------------------------------------------------------------
# Build the prior-day predictors for the gap-onset models.
#
# For every predictor p (e.g. log steps) we create, per observed day:
#   x_p   the transformed value
#   w_p   WITHIN-person deviation:   (x - person's own mean) / pooled within-person SD
#         "how unusual was yesterday for THIS person" (primary predictor)
#   b_p   BETWEEN-person level:      (person's mean - average person) / SD of person means
#         "is this person generally high or low" (keeps the two effects apart)
#   tw_p  same as w_p but using the person's TRAILING baseline (mean of the previous
#         14 observed days) instead of the whole-follow-up mean (sensitivity: no look-ahead)
#
# Scaling constants are computed once on the pooled data of the stream and stored in
# attr(result, "scales"), so odds ratios are "per 1 within-person SD" everywhere.
# ==============================================================================

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
predictor_table <- function(days, pred_names, cfg, trailing = TRUE, keep_x = TRUE) {
  dt <- days[, c("person_id", "date"), with = FALSE]
  scales <- vector("list", length(pred_names))

  for (i in seq_along(pred_names)) {
    p    <- pred_names[i]
    spec <- cfg$predictors[[p]]
    assert_columns(days, spec$col, paste("day table for predictor", p))
    x_nm <- paste0("x_", p)
    dt[, (x_nm) := apply_transform(days[[spec$col]], spec$transform)]

    # person mean over the person's observed days
    dt[, m_i := mean(get(x_nm), na.rm = TRUE), by = person_id]
    dt[is.nan(m_i), m_i := NA_real_]

    # pooled within-person SD and between-person SD (one constant per predictor)
    dev  <- dt[[x_nm]] - dt$m_i
    n_ok <- sum(!is.na(dev))
    n_p  <- dt[!is.na(m_i), uniqueN(person_id)]
    sd_w <- sqrt(sum(dev^2, na.rm = TRUE) / max(n_ok - n_p, 1L))
    pm   <- unique(dt[!is.na(m_i), .(person_id, m_i)])
    sd_b <- stats::sd(pm$m_i)
    mu_b <- mean(pm$m_i)

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
                              mean_of_person_means = mu_b)
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
#            partial-wear flags and person covariates
#   scales : within/between SDs used for scaling
#
# `observed_filter` (optional, e.g. "has_main_sleep") restricts which sleep/hr/activity
# rows count as "observed" (used for the main-sleep-only sensitivity).
# trailing = TRUE also builds the trailing-baseline columns (only the sensitivity analyses use them).
prepare_stream_bundle <- function(stream, day_tables, persons, cfg, observed_filter = NULL, trailing = FALSE) {
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
  own_tab   <- predictor_table(days, own_names, cfg, trailing = trailing, keep_x = FALSE)
  scales    <- attr(own_tab, "scales")
  own_tab[, pw_own := partial_wear_flag(days, cfg$onset$partial_wear[[stream]])]
  add_columns(own_tab, c("person_id", "date"))

  # activity predictors for the cross-stream models
  if (stream != "activity") {
    act_days <- day_tables[["activity"]]
    act_tab  <- predictor_table(act_days, cfg$activity_predictors, cfg, trailing = trailing, keep_x = FALSE)
    act_tab[, pw_act := partial_wear_flag(act_days, cfg$onset$partial_wear$activity)]
    add_columns(act_tab, c("person_id", "date"))
    scales <- rbind(scales, attr(act_tab, "scales"))
    rm(act_tab, act_days)
  } else {
    rs[, pw_act := pw_own]             # activity: own stream IS the activity stream
  }
  rm(own_tab)

  # person covariates
  cov <- persons[, .(person_id, cohort, sex2 = factor(as.character(sex), levels = c("Female", "Male")),
                     age10 = (age_first - 60) / 10)]
  add_columns(cov, "person_id")
  setkey(rs, person_id, date)

  list(stream = stream, gaps = gaps, rs = rs, scales = scales,
       n_persons = uniqueN(rs$person_id), trailing = trailing)
}
