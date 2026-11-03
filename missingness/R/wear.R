# ==============================================================================
# wear.R
# ------------------------------------------------------------------------------
# Wear time of the LAST recorded day before a gap (d), and the table that follows the rows of the wear sensitivity analyses.
#
# Measure. One wear measure for every stream: the heart-rate minutes recorded on the calendar day (`total_hr_min` = sum of the four zone minutes, at most
# 1440), cfg$wear. The activity table has no wear time of its own (Fitbit counts non-wear as sedentary), so heart-rate minutes are used for the activity,
# sleep and heart-rate streams alike. A day with NO heart-rate record has no wear: it is unknown, excluded and counted, never treated as a short day.
# The minutes are a lower bound of the time worn (a worn watch can lose the signal).
#
# Days. In a stream's risk set (every observed day except the person's last) a day d is
#   "last"     a gap starts on d + 1 (the last day with a record before the gap), or
#   "ordinary" d + 1 is observed too. Ordinary days include the first day back after a gap (the return day, often a partial day).
# Description (stage "wear", no model): the distribution of the wear of both kinds of day by stream and cohort, the share under each threshold of
# cfg$wear$hours with a person-clustered interval, and the difference last - ordinary. The sensitivity analyses wear_ge<h>h are in sensitivity.R.
# ==============================================================================

WEAR_GROUP_LABELS <- c(last = "Last day before a gap (d)", ordinary = "Ordinary day (a record on d + 1)")
WEAR_STREAM_LABELS <- c(activity = "Activity", sleep = "Sleep", hr = "Heart rate")

# One row per risk-set day of a stream with the heart-rate minutes recorded on it.
#   days   : the stream's day table (person_id, date)     wear : person_id, date, wear_min (heart-rate minutes; absent day = no record)
#   persons: person_id, cohort
wear_days <- function(days, wear, persons, stream) {
  d <- days[, .(person_id, date)]
  setorder(d, person_id, date)
  d[, next_date := shift(date, type = "lead"), by = person_id]
  d <- d[!is.na(next_date)]                                  # the last day with a record has no "next day"
  d[, `:=`(last_day = as.integer(next_date - date) > 1L, next_date = NULL)]
  d[wear, wear_min := i.wear_min, on = c("person_id", "date")]
  d[as.data.table(persons)[, .(person_id, cohort)], cohort := i.cohort, on = "person_id"]
  d[, stream := stream]
  d[]
}

# Distribution, summary and contrast of one stream (the output of wear_days).
#   distribution : stream, cohort, group, wear_min (whole minutes, 0..max_min), n, cum_share  (the ECDF, one row per minute value present)
#   summary      : stream, cohort, group, n_days, n_people, n_no_hr_record, quartiles of the hours, and per threshold the share below it
#                  (% of days with a record, 95% interval with people as clusters)
#   contrast     : stream, cohort, quantity, last, ordinary, difference (last - ordinary), 95% interval, p (cluster-robust, CR1), p_holm (over the cohorts)
wear_describe_stream <- function(x, cfg) {
  hrs <- cfg$wear$hours; cap <- cfg$wear$max_min
  x <- copy(x)[, group := fifelse(last_day, "last", "ordinary")]
  known <- x[!is.na(wear_min)]
  dist <- known[, .(n = .N), by = .(stream, cohort, group, wear_min = pmin(pmax(round(wear_min), 0), cap))][order(stream, cohort, group, wear_min)]
  dist[, cum_share := cumsum(n) / sum(n), by = .(stream, cohort, group)]
  loc <- x[, {
    k <- wear_min[!is.na(wear_min)]
    q <- if (length(k)) stats::quantile(k, c(0.25, 0.5, 0.75), names = FALSE) / 60 else rep(NA_real_, 3L)
    .(n_days = .N, n_people = uniqueN(person_id), n_no_hr_record = sum(is.na(wear_min)), n_over_max = sum(k > cap),
      q1_h = q[1L], median_h = q[2L], q3_h = q[3L])
  }, by = .(stream, cohort, group)]
  below <- rbindlist(lapply(hrs, function(h) known[, {
    ci <- cluster_mean_ci(as.numeric(wear_min < h * 60), person_id)
    .(threshold_h = h, n_known = .N, n_below = sum(wear_min < h * 60), pct_below = 100 * ci$mean, pct_low = 100 * ci$lo, pct_high = 100 * ci$hi)
  }, by = .(stream, cohort, group)]))
  contr <- rbindlist(lapply(c(hrs, NA_real_), function(h) known[, {
    v <- if (is.na(h)) wear_min / 60 else 100 * as.numeric(wear_min < h * 60)
    r <- day_contrast(v, as.integer(last_day), person_id)
    if (is.null(r)) NULL
    else .(quantity = if (is.na(h)) "mean wear (hours)" else sprintf("%% of days under %g h", h),
           last = r$mean_gap, ordinary = r$mean_no_gap, difference = r$difference, low = r$low, high = r$high, p = r$p,
           n_last = sum(last_day), n_ordinary = sum(!last_day))
  }, by = .(stream, cohort)]))
  if (nrow(contr)) contr[, p_holm := stats::p.adjust(p, "holm"), by = .(stream, quantity)]
  list(distribution = dist, summary = merge(loc, below, by = c("stream", "cohort", "group")), contrast = contr)
}

# ECDF of the wear of the last day (solid) and of ordinary days (dashed): streams side by side, cohorts in colour.
plot_wear_last_day <- function(dist, cfg) {
  d <- copy(dist)
  d[, `:=`(cohort = factor(cohort, levels = cfg$cohorts), hours = wear_min / 60,
           stream = factor(stream, levels = names(WEAR_STREAM_LABELS), labels = WEAR_STREAM_LABELS),
           group = factor(group, levels = names(WEAR_GROUP_LABELS), labels = WEAR_GROUP_LABELS))]
  ggplot(d, aes(x = hours, y = cum_share, colour = cohort, linetype = group)) +
    geom_vline(xintercept = cfg$wear$hours, colour = "grey55", linewidth = 0.4, linetype = "dotted") +
    geom_step(linewidth = 0.8) +
    facet_wrap(~stream, nrow = 1) +
    scale_colour_manual(values = cfg$palette$cohort) +
    scale_linetype_manual(values = c("solid", "dashed")) +
    scale_x_continuous(breaks = seq(0, 24, 4), limits = c(0, 24), expand = expansion(mult = c(0.01, 0.01))) +
    scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"), limits = c(0, 1)) +
    labs(title = "Wear Time on the Last Day Before a Gap",
         subtitle = ov_wrap(sprintf("Heart-rate hours recorded on day d (the same measure for all streams). Solid = d is the last day with a record before a gap; dashed = ordinary days (a record on d + 1, return days included). Y = share of days with at most x hours. Dotted lines: %s h. Days with no heart-rate record are not drawn (counted in wear_last_day.tsv).",
                            paste(cfg$wear$hours, collapse = " and ")), 120),
         x = "Heart-rate hours recorded on day d", y = "Cumulative share of days", colour = "Cohort", linetype = NULL) +
    theme_bw(base_size = 11) +
    theme(strip.text = element_text(face = "bold"), legend.position = "bottom", legend.box = "vertical")
}

# Rows and events of the day-d model and of the wear analyses, so the effect of the filter on the data is visible next to the estimates:
# the primary rows (predictors of d - 1), every row with complete predictors on d, the same rows by wear class of d, and the rows kept by each analysis.
# mean_w = mean of the primary term on those rows (the within-person z-score of steps; a filter on full days shifts it above the primary model's value).
wear_filter_flow <- function(bundle, spec, cfg, term = cfg$primary_term) {
  spec <- as.list(spec); spec$predictors <- list(spec$predictors[[1L]])
  hrs <- cfg$wear$hours; ycol <- paste0("y_ge", spec$k)
  p <- build_onset_data(bundle, spec, cfg)
  d <- build_onset_data(bundle, spec, cfg, list(day = "day_before_gap", keep_wear = TRUE))
  row <- function(label, x) data.table(stream = bundle$stream, rows = label, n_rows = nrow(x), n_people = uniqueN(x$person_id),
                                       n_events = sum(x[[ycol]]), event_pct = if (nrow(x)) 100 * mean(x[[ycol]]) else NA_real_,
                                       mean_w = if (nrow(x)) mean(x[[term]]) else NA_real_)
  lab <- c(paste0("under ", hrs[1L], " h"), if (length(hrs) > 1L) paste0(hrs[-length(hrs)], " to under ", hrs[-1L], " h"), paste0(hrs[length(hrs)], " h or more"))
  bin <- cut(d$wear_min, c(-Inf, hrs * 60, Inf), labels = lab, right = FALSE)
  out <- c(list(row("Primary model rows (predictors of the last full day d - 1)", p),
                row("Day d: every row with complete predictors", d),
                row("  no heart-rate record on d (wear unknown, excluded)", d[is.na(wear_min)])),
           lapply(lab, function(l) row(paste0("  wear on d ", l), d[which(bin == l)])),
           lapply(hrs, function(h) row(sprintf("Kept by wear_ge%gh (wear on d >= %g h)", h, h), d[which(wear_min >= h * 60)])))
  rbindlist(out)
}
