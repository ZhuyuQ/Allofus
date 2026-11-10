# ==============================================================================
# plots_stream_states.R - figures of the stream-state analysis (stream_states.R). Cohort colours as in the other figures; estimates are drawn with their
# person-clustered 95% intervals; cells with 1-19 people are not drawn (their values are NA). One figure per question, named stream_<part>_<what>.
# ==============================================================================

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

ST_STATE_COLORS <- c(`All present` = "#2E8B57", `All missing` = "#7F7F7F", `Non-wear` = "#E69F00", `Night missing` = "#0072B2", Other = "#CC3311")
ST_REASON_COLORS <- c(`analysed days` = "#2E8B57", `stream not active` = "#999999", `device without heart rate` = "#E69F00", `device unknown, removed` = "#CC3311")

st_theme <- function() theme_bw(base_size = 10) + theme(strip.text = element_text(face = "bold"), legend.position = "bottom")

# A dodged point-and-interval layer for the share tables (x = estimate).
st_pointrange <- function(width = 0.6) geom_pointrange(aes(xmin = low, xmax = high), position = position_dodge(width), size = 0.4)

st_cohorts <- function(d, pcfg) { d <- copy(d); d[, cohort := factor(cohort, levels = pcfg$cohorts)]; d[] }

# ---- Part A -------------------------------------------------------------------------------------------------------------------------------

# Flow from all days in the windows to the days analysed, as shares of the window days (stacked).
plot_st_flow <- function(flow, pcfg) {
  keep <- c(`analysed days` = "analysed days", `removed: a stream not active` = "stream not active", `removed: device without heart rate` = "device without heart rate",
            `removed: device unknown (usual device without heart rate, or none)` = "device unknown, removed")
  d <- st_cohorts(flow[step %in% names(keep) & !is.na(share_pct)], pcfg)[, reason := factor(keep[step], levels = names(ST_REASON_COLORS))]
  ggplot(d, aes(cohort, share_pct, fill = reason)) + geom_col(width = 0.7) + scale_fill_manual(values = ST_REASON_COLORS) +
    labs(title = "Part A. From All Days in the Windows to the Days Analysed", x = NULL, y = "% of the days in the windows", fill = NULL,
         caption = "A day is removed for the first of these reasons, in this order. Window = first to last record in any stream.") + st_theme()
}

# Share of days removed by each reason, by cohort, with person-clustered intervals.
plot_st_reasons <- function(flow, pcfg) {
  keep <- c(`removed: a stream not active` = "stream not active", `removed: device without heart rate` = "device without heart rate",
            `removed: device unknown (usual device without heart rate, or none)` = "device unknown, removed",
            `  of which device unknown, kept (usual device measures heart rate)` = "device unknown, kept (check)")
  d <- st_cohorts(flow[step %in% names(keep) & !is.na(share_pct)], pcfg)[, reason := factor(keep[step], levels = rev(keep))]
  ggplot(d, aes(share_pct, reason, color = cohort)) + st_pointrange() + scale_color_manual(values = pcfg$palette$cohort) +
    labs(title = "Part A. Share of Days Removed by Each Reason", x = "% of the days in the windows (95% interval, people as clusters)", y = NULL, color = "Cohort") + st_theme()
}

# ---- Part B -------------------------------------------------------------------------------------------------------------------------------

plot_st_filter <- function(tab, pcfg) {
  d <- st_cohorts(tab[!is.na(share_pct)], pcfg)
  ggplot(d, aes(share_pct, cohort, color = cohort)) + st_pointrange() + facet_wrap(~step, ncol = 1, scales = "free_x") + scale_color_manual(values = pcfg$palette$cohort, guide = "none") +
    labs(title = "Part B. Activity Rows Affected by the Filter", x = "% of activity rows (95% interval, people as clusters)", y = NULL) + st_theme()
}

# Spread across people of the share of a person's activity rows set to missing (ECDF; people with at least one activity row).
plot_st_filter_person <- function(per, pcfg) {
  d <- st_cohorts(per, pcfg)
  ggplot(d, aes(100 * share, color = cohort)) + stat_ecdf(linewidth = 0.8) + scale_x_sqrt() + scale_color_manual(values = pcfg$palette$cohort) +
    labs(title = "Part B. Spread of the Filtered Share Across People", x = "% of a person's activity rows set to missing (square-root scale)", y = "Cumulative share of people", color = "Cohort") + st_theme()
}

# ---- Part C -------------------------------------------------------------------------------------------------------------------------------

# Share of days in each of the eight combinations, before and after the filter.
plot_st_combos <- function(combos, pcfg) {
  d <- st_cohorts(combos[!is.na(share_pct) & share_pct > 0], pcfg)
  d[, label := factor(label, levels = unique(ST_COMBO_LABELS))][, stage := factor(stage, levels = c("before the filter", "after the filter"))]
  ggplot(d, aes(share_pct, label, color = stage, shape = category)) + geom_pointrange(aes(xmin = pmax(low, 1e-4), xmax = high), position = position_dodge(0.6), size = 0.35) +
    facet_wrap(~cohort, ncol = 2) + scale_x_log10() + scale_color_manual(values = c("#7F7F7F", "#0072B2")) + scale_y_discrete(limits = rev(levels(d$label))) +
    labs(title = "Part C. Combinations of the Three Streams", x = "% of days (log scale; 95% interval, people as clusters)", y = NULL, color = NULL, shape = NULL,
         caption = "Expected = activity >= heart rate >= sleep (nested). 'Activity + sleep, heart rate missing' is shown separately from the other unexpected combinations.") + st_theme()
}

plot_st_checks <- function(checks, pcfg) {
  d <- st_cohorts(checks[!is.na(share_pct)], pcfg)
  ggplot(d, aes(share_pct, cohort, color = cohort)) + st_pointrange() + facet_wrap(~check, ncol = 1, scales = "free_x") + scale_color_manual(values = pcfg$palette$cohort, guide = "none") +
    labs(title = "Part C. Nesting Checks Before the Filter (both should be close to zero)", x = "% of days (95% interval, people as clusters)", y = NULL) + st_theme()
}

# ---- Part D -------------------------------------------------------------------------------------------------------------------------------

plot_st_steps_state3 <- function(steps, pcfg) {
  d <- st_cohorts(steps, pcfg)
  ggplot(d, aes(steps, color = cohort)) + stat_ecdf(linewidth = 0.8, n = 500) +
    scale_x_continuous(trans = scales::pseudo_log_trans(base = 10), breaks = c(0, 100, 1000, 10000)) + scale_color_manual(values = pcfg$palette$cohort) +
    labs(title = "Part D. Steps on State-3 Days (activity present, heart rate and sleep missing)", x = "Steps (pseudo-log scale)", y = "Cumulative share of state-3 days", color = "Cohort") + st_theme()
}

# Mean within-person difference (state 3 minus state 1 days, same type of day) with its interval over people.
plot_st_state3_diff <- function(cmp, pcfg) {
  d <- st_cohorts(cmp[!is.na(mean)], pcfg)[, metric := factor(metric, c("steps", "mvpa"), c("Steps per day", "MVPA minutes per day"))]
  ggplot(d, aes(mean, cohort, color = day_type)) + geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") + st_pointrange() + facet_wrap(~metric, scales = "free_x") +
    scale_color_brewer(palette = "Dark2") + labs(title = "Part D. State-3 Days Compared With the Person's Own State-1 Days", x = "Mean difference, state 3 minus state 1 (95% interval over people)", y = NULL, color = "Type of day",
         caption = "Within the same type of day (weekday / weekend); people with at least three days of each state. Unadjusted.") + st_theme()
}

# ---- Part E -------------------------------------------------------------------------------------------------------------------------------

# State shares by cohort: all days (long follow-up weighs more) and person average (every person weighs the same); stratum = all, age_group or device_family.
plot_st_composition <- function(comp, by, pcfg) {
  d <- st_cohorts(comp[stratum == by & !is.na(share_pct)], pcfg)[, state_label := factor(state_label, levels = ST_STATE_LABELS)]
  p <- ggplot(d, aes(share_pct, state_label, color = cohort)) + st_pointrange() + scale_y_discrete(limits = rev(ST_STATE_LABELS)) + scale_color_manual(values = pcfg$palette$cohort) +
    labs(title = sprintf("Part E. Share of Days in Each State%s", if (by == "all") "" else paste(" by", sub("_", " ", by))), x = "% of days (95% interval)", y = NULL, color = "Cohort",
         caption = "Unadjusted: cohorts differ in age, device and length of follow-up; read as patterns, not effects.") + st_theme()
  if (by == "all") p + facet_wrap(~weighting) else p + facet_grid(level ~ weighting, scales = "free_x")
}

# Day-by-day states of a random sample of people per cohort (one row per person; grey = day removed).
plot_st_day_by_day <- function(s, pcfg) {
  d <- st_cohorts(s, pcfg)[, state_label := factor(fifelse(state == 0L, "removed", ST_STATE_LABELS[as.character(state)]), levels = c("removed", ST_STATE_LABELS))]
  ggplot(d, aes(day, row, fill = state_label)) + geom_raster() + facet_grid(cohort ~ ., scales = "free_y", space = "free_y") +
    scale_fill_manual(values = c(removed = "white", ST_STATE_COLORS), drop = FALSE) + scale_y_reverse() +
    labs(title = "Part E. Day-by-Day States of a Random Sample of People", x = "Day since the start of the window", y = "People (one row each)", fill = NULL,
         caption = "Individual-level picture: check the disclosure rules before it leaves the workspace.") + st_theme() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
}

plot_st_night_habit <- function(persons, pcfg) {
  d <- st_cohorts(persons, pcfg)
  ggplot(d, aes(cohort, 100 * share, fill = cohort)) + geom_violin(alpha = 0.5, color = NA, scale = "width") + geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
    scale_fill_manual(values = pcfg$palette$cohort, guide = "none") +
    labs(title = "Part E. Habitual Night Non-Wear", x = NULL, y = "% of a person's worn days that are state 4 (night missing)",
         caption = "People with enough worn days; computed from the whole follow-up, so for description only.") + st_theme()
}

# ---- Part F -------------------------------------------------------------------------------------------------------------------------------

# Heart-rate hours of the day before and of the day itself, for state-4 nights next to state-1 nights.
plot_st_nights_hr <- function(nights, pcfg, min_people = 20L) {
  d <- st_cohorts(nights[, .(cohort, person_id, state, x = pmin(floor(hr_prev), 23), y = pmin(floor(hr_now), 23))], pcfg)
  d[, state_label := factor(ST_STATE_LABELS[as.character(state)], levels = c("All present", "Night missing"))]
  b <- d[, .(n = .N, people = uniqueN(person_id)), by = .(cohort, state_label, x, y)][, share := n / sum(n), by = .(cohort, state_label)][people < min_people, share := NA_real_]
  ggplot(b, aes(x, y, fill = share)) + geom_tile() + facet_grid(state_label ~ cohort) + scale_fill_viridis_c(trans = "sqrt", na.value = "grey85", labels = scales::percent) +
    labs(title = "Part F. Heart-Rate Wear on the Day Before and on the Day Itself", x = "Heart-rate hours, day before", y = "Heart-rate hours, the day", fill = "share of nights",
         caption = "State-4 nights next to state-1 nights; only nights where both days were worn days. Grey: fewer than 20 people.") + st_theme()
}

plot_st_probably_worn <- function(tab, pcfg) {
  d <- st_cohorts(tab[!is.na(share_pct)], pcfg)[, cutoff := factor(paste(cutoff_hours, "hours on both days"))]
  ggplot(d, aes(share_pct, cohort, color = state_label)) + st_pointrange() + facet_wrap(~cutoff, ncol = 1) + scale_color_manual(values = ST_STATE_COLORS) +
    labs(title = "Part F. Nights That Were Probably Worn", x = "% of nights with heart-rate wear at least this high on both days (95% interval)", y = NULL, color = NULL,
         caption = "High wear on both days makes night wear likely but does not prove it.") + st_theme()
}

# ---- Part G -------------------------------------------------------------------------------------------------------------------------------

plot_st_hr_wear <- function(hours, below, pcfg, cuts = c(10, 20)) {
  d <- st_cohorts(hours, pcfg)
  p1 <- ggplot(d, aes(hrh, color = cohort)) + stat_ecdf(linewidth = 0.8, n = 300) + geom_vline(xintercept = cuts, linetype = "dotted", color = "grey40") +
    scale_color_manual(values = pcfg$palette$cohort) + scale_x_continuous(breaks = seq(0, 24, 4)) +
    labs(title = "Part G. Heart-Rate Wear on Worn Days (states 1 and 4)", x = "Heart-rate hours recorded", y = "Cumulative share of days", color = "Cohort") + st_theme()
  p1
}

plot_st_hr_below <- function(below, pcfg) {
  d <- st_cohorts(below[!is.na(share_pct)], pcfg)[, cutoff := factor(paste("below", cutoff_hours, "hours"))]
  ggplot(d, aes(share_pct, cohort, color = cohort)) + st_pointrange() + facet_wrap(~cutoff, ncol = 1, scales = "free_x") + scale_color_manual(values = pcfg$palette$cohort, guide = "none") +
    labs(title = "Part G. Worn Days With Little Heart-Rate Wear", x = "% of worn days (95% interval, people as clusters)", y = NULL) + st_theme()
}

plot_st_minutes <- function(tab, pcfg, limit_pct = 0.1) {
  d <- st_cohorts(tab[!is.na(share_pct)], pcfg)
  ggplot(d, aes(share_pct, cohort, color = cohort)) + geom_vline(xintercept = limit_pct, linetype = "dashed", color = "grey40") + st_pointrange() +
    scale_color_manual(values = pcfg$palette$cohort, guide = "none") +
    labs(title = "Part G. Heart-Rate Records Above the Length of the Day (check)", x = "% of heart-rate records (dashed line = the limit; should be almost zero)", y = NULL) + st_theme()
}

# ---- all figures of a run -----------------------------------------------------------------------------------------------------------------

# The figures of parts A-C (abc: the result of st_run_abc, with `person_share` for part B) and, when given, D-G (dg: st_run_dg); a named list of ggplots.
st_figures <- function(abc, dg = NULL, pcfg, cfg) {
  out <- list(stream_a_flow = plot_st_flow(abc$flow, pcfg), stream_a_reasons = plot_st_reasons(abc$flow, pcfg),
              stream_b_filter_share = plot_st_filter(abc$activity_filter, pcfg), stream_b_filter_people = plot_st_filter_person(abc$person_share, pcfg),
              stream_c_combos = plot_st_combos(abc$combos, pcfg), stream_c_checks = plot_st_checks(abc$nesting_checks, pcfg),
              stream_g_minutes_check = plot_st_minutes(abc$minutes_check, pcfg, 100 * cfg$minutes_exceed_max))
  if (is.null(dg)) return(out)
  c(out, list(stream_d_steps_state3 = plot_st_steps_state3(dg$steps_data, pcfg), stream_d_state3_vs_state1 = plot_st_state3_diff(dg$state3_vs_state1, pcfg),
              stream_e_composition = plot_st_composition(dg$composition, "all", pcfg), stream_e_composition_age = plot_st_composition(dg$composition, "age_group", pcfg),
              stream_e_composition_device = plot_st_composition(dg$composition, "device_family", pcfg), stream_e_day_by_day = plot_st_day_by_day(dg$day_sample, pcfg),
              stream_e_night_habit = plot_st_night_habit(dg$night_habit_persons, pcfg), stream_f_nights_hr = plot_st_nights_hr(dg$nights, pcfg, cfg$min_people),
              stream_f_probably_worn = plot_st_probably_worn(dg$probably_worn, pcfg), stream_g_hr_wear = plot_st_hr_wear(dg$hr_hours, dg$hr_wear_below, pcfg, cfg$hr_cuts_present),
              stream_g_hr_below = plot_st_hr_below(dg$hr_wear_below, pcfg)))
}
