# ==============================================================================
# plots_weekday_rhythm.R
# ------------------------------------------------------------------------------
# Figures of the weekly rhythm of gap starts and gap ends (weekday_rhythm.R). Shared rules: weekday axis Monday to Sunday with Saturday and Sunday
# shaded; columns are streams (activity, sleep, heart rate); gap start and gap end keep one colour; ratios on a log axis with a reference line at 1;
# intervals follow cfg$plot$interval_style. Cells with 1-19 people are not drawn (grey); cells with few days at risk are grey too.
# ==============================================================================

WR_COLORS <- c(`gap start` = "#2C3E50", `gap end` = "#16A085")
WR_STREAMS <- c(activity = "Activity", sleep = "Sleep", hr = "Heart rate")
WR_SHADE <- annotate("rect", xmin = 5.5, xmax = 7.5, ymin = -Inf, ymax = Inf, fill = "grey85", alpha = 0.5)

# Weekday index (0-6) -> factor Mon..Sun; stream and outcome factors in the fixed order.
wr_label <- function(d) {
  d <- copy(d)
  if ("weekday" %in% names(d) && is.numeric(d$weekday)) d[, weekday := factor(WR_WEEKDAYS[weekday + 1L], levels = WR_WEEKDAYS)]
  if ("stream" %in% names(d)) d[, stream := factor(WR_STREAMS[as.character(stream)], levels = WR_STREAMS)]
  d[]
}

wr_theme <- function() theme_bw(base_size = 10) + theme(strip.text = element_text(face = "bold"), legend.position = "bottom", axis.text.x = element_text(angle = 0))

# P-value text for a test row.
wr_ptext <- function(p) ifelse(is.na(p), "NA", ifelse(p < 1e-4, "< 0.0001", sprintf("%.4f", p)))

# Figure 1 (H1): predicted % of gap starts / gap ends by weekday (S1, E1), crude rates behind, dashed weekly average.
plot_wr_profile <- function(profile, summ, tests, crude, cfg) {
  style <- cfg$plot$interval_style %||% "se"
  key <- function(x) x[fit_id %in% c("S1", "E1")][, outcome := fifelse(fit_id == "S1", "gap start", "gap end")]
  d <- wr_label(key(profile)); add_interval_columns(d, "pct", "pct_low", "pct_high", style = style)
  s <- wr_label(key(summ)); tt <- tests[id %in% c("H1s", "H1e"), .(stream, outcome, p_holm)]
  s <- merge(s, tt[, .(stream = factor(WR_STREAMS[stream], levels = WR_STREAMS), outcome, p_holm)], by = c("stream", "outcome"), all.x = TRUE)
  s[, label := sprintf("peak-to-trough %.2f (%.2f-%.2f)\nH1 Holm p %s", peak_to_trough, ptr_low, ptr_high, wr_ptext(p_holm))]
  cr <- wr_label(crude[cohort == "All" & !is.na(rate_pct)])
  ggplot(d, aes(x = weekday, y = est_)) + WR_SHADE +
    geom_hline(data = s, aes(yintercept = average_pct), linetype = "dashed", color = "grey40") +
    geom_point(data = cr, aes(y = rate_pct), shape = 1, size = 2, color = "grey35") +
    interval_layer(style, horizontal = FALSE, color_var = "outcome", width = 0.4) +
    geom_text(data = s, aes(x = 1, y = Inf, label = label), hjust = 0, vjust = 1.1, size = 2.8, inherit.aes = FALSE) +
    facet_grid(outcome ~ stream, scales = "free_y") +
    scale_color_manual(values = WR_COLORS, guide = "none") + scale_fill_manual(values = WR_COLORS, guide = "none") +
    labs(title = "Weekly Rhythm of Gap Starts and Gap Ends", x = NULL, y = "Predicted % (population average)", caption = interval_caption(style)) + wr_theme()
}

# Figure 2, panel A: heat map of weekday x days since return (log2 ratio to the overall rate; grey = suppressed or few days) with the crude % by weekday on top
# and by days since return on the right.
wr_heat_panels <- function(ce) {
  lim <- max(abs(ce$fill), 0.5, na.rm = TRUE)
  heat <- ggplot(ce, aes(weekday, days_since_return, fill = fill)) + geom_tile(color = "white", linewidth = 0.2) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-lim, lim), na.value = "grey80",
                         name = "ratio to the overall rate", breaks = c(-lim, 0, lim), labels = function(x) sprintf("%.2f", 2^x)) +
    labs(x = NULL, y = "Days since return") + theme_minimal(base_size = 9) + theme(legend.position = "bottom", panel.grid = element_blank())
  bar <- function(d, x) ggplot(d, aes(.data[[x]], r * 100)) + geom_col(fill = WR_COLORS[["gap start"]], width = 0.7) + labs(x = NULL, y = "% gap starts") + theme_minimal(base_size = 8)
  top <- bar(ce[, .(r = sum(n_events, na.rm = TRUE) / sum(n_days, na.rm = TRUE)), by = weekday], "weekday")
  right <- bar(ce[, .(r = sum(n_events, na.rm = TRUE) / sum(n_days, na.rm = TRUE)), by = days_since_return], "days_since_return") + coord_flip() + theme(axis.text.y = element_blank())
  gridExtra::arrangeGrob(top, grid::nullGrob(), heat, right, layout_matrix = rbind(c(1, 2), c(3, 4)), widths = c(4, 1), heights = c(1, 4))
}

# Figure 2, panel B: ratio to the average day without (S1) and with (S3) days since return, joined by a line.
wr_dumbbell_panel <- function(pr, sh, h2) {
  lv <- c("weekday only (S1)", "weekday + days since return (S3)")
  pr <- copy(pr)[, model := factor(fifelse(fit_id == "S1", lv[1L], lv[2L]), levels = lv)]
  ggplot(pr, aes(ratio_to_average_day, weekday)) + geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
    geom_line(aes(group = weekday), color = "grey60") + geom_point(aes(color = model), size = 2.5) + scale_x_log10() + scale_y_discrete(limits = rev(WR_WEEKDAYS)) +
    scale_color_manual(values = c("grey50", WR_COLORS[["gap start"]]), name = NULL) +
    labs(x = "Ratio to the average day", y = NULL, subtitle = sprintf("share of the pattern explained by days since return: %s; H2 Holm p %s",
         ifelse(nrow(sh) > 0L, sprintf("%.0f%%", 100 * sh$shrinkage), "NA"), wr_ptext(h2))) + wr_theme() + theme(plot.subtitle = element_text(size = 7))
}

# Figure 2, panel C: the days-since-return curve without (S2) and with (S3) the weekday.
wr_curve_panel <- function(cv, lv) {
  cv <- copy(cv)[, level := factor(level, levels = lv)][, model := fifelse(fit_id == "S2", "days since return only (S2)", "weekday + days since return (S3)")]
  cols <- c("grey50", WR_COLORS[["gap start"]])
  ggplot(cv, aes(level, pct, group = model, color = model, fill = model)) + geom_ribbon(aes(ymin = pct_low, ymax = pct_high), alpha = 0.2, color = NA) + geom_line() +
    scale_color_manual(values = cols, name = NULL) + scale_fill_manual(values = cols, name = NULL) +
    labs(x = "Days since return", y = "Predicted % of gap starts") + wr_theme() + theme(axis.text.x = element_text(size = 6))
}

# Figure 2 (H2), one stream (s): A heat map with margins, B weekday ratios before / after days since return, C the days-since-return curve.
plot_wr_calendar_vs_cycle <- function(s, cells, profile, curve, shrink, tests, cfg) {
  ce <- wr_label(cells[stream == s]); lv <- levels(ce$days_since_return)
  ce[, `:=`(fill = fifelse(suppressed | low_exposure, NA_real_, log2(ratio_to_overall)), days_since_return = factor(days_since_return, levels = rev(lv)))]
  B <- wr_dumbbell_panel(wr_label(profile[stream == s & fit_id %in% c("S1", "S3")]), shrink[stream == s], tests[id == "H2" & stream == s, p_holm])
  C <- wr_curve_panel(curve[stream == s & fit_id %in% c("S2", "S3")], lv)
  out <- gridExtra::arrangeGrob(wr_heat_panels(ce), gridExtra::arrangeGrob(B, C, ncol = 2), nrow = 2, heights = c(1.5, 1),
                                top = grid::textGrob(sprintf("Calendar Routine or Charging Cycle? (%s stream)", WR_STREAMS[[s]]), x = 0.01, hjust = 0, gp = grid::gpar(fontface = "bold")))
  attr(out, "size") <- c(11, 11)
  out
}

# Figure 3 (H3), one stream: personal first-harmonic rhythm (a + u_c, b + u_s); angle = peak weekday, radius = amplitude. Binned BLUPs (bins under
# min_people suppressed), the shared rhythm (arrow), the model-based 95% region (ellipse).
plot_wr_personal <- function(s, personal, bins, ellipse, tests) {
  p <- personal[personal$stream == s]; b <- bins[bins$stream == s & !is.na(bins$n)]; el <- ellipse[ellipse$stream == s]
  r <- max(abs(c(b$x, b$y, el$x, el$y)), 0.1) * 1.15; ang <- 2 * pi * 0:6 / 7
  grid_c <- rbindlist(lapply(seq(r / 3, r, length.out = 3), function(rr) data.table(x = rr * cos(seq(0, 2 * pi, length.out = 101)), y = rr * sin(seq(0, 2 * pi, length.out = 101)), rr = rr)))
  h3 <- tests[id == "H3" & stream == s, p_holm]
  ggplot() + geom_path(data = grid_c, aes(x, y, group = rr), color = "grey85") +
    geom_segment(data = data.table(ang), aes(x = 0, y = 0, xend = r * cos(ang), yend = r * sin(ang)), color = "grey85") +
    geom_tile(data = b, aes(x, y, fill = n), width = unique(b$width)[1L], height = unique(b$width)[1L]) +
    geom_path(data = el, aes(x, y), color = "black", linetype = "dashed") +
    geom_segment(aes(x = 0, y = 0, xend = p$shared_a, yend = p$shared_b), arrow = grid::arrow(length = grid::unit(0.15, "cm")), color = "#C0392B", linewidth = 0.9) +
    geom_text(data = data.table(ang, lab = WR_WEEKDAYS), aes(x = 1.08 * r * cos(ang), y = 1.08 * r * sin(ang), label = lab), size = 3) +
    scale_fill_gradient(low = "#D6EAF8", high = "#1B4F72", name = "people per bin") + coord_fixed(xlim = c(-1.2 * r, 1.2 * r), ylim = c(-1.2 * r, 1.2 * r)) +
    labs(title = sprintf("Personal Weekly Rhythm of Gap Starts (%s stream)", WR_STREAMS[[s]]),
         subtitle = sprintf("person peak-to-trough, 10th-90th percentile: %.2f-%.2f; H3 Holm p %s", p$person_ptr_p10, p$person_ptr_p90, wr_ptext(h3)), x = NULL, y = NULL) +
    theme_void(base_size = 10) + theme(plot.subtitle = element_text(size = 8), legend.position = "bottom")
}

# Figure 4 (H4): the weekday profile per cohort (S4 starts, E4 ends); ratio to the cohort's own average day, or % when level = TRUE.
plot_wr_cohort <- function(profile, tests, cfg, level = FALSE) {
  style <- cfg$plot$interval_style %||% "se"
  d <- wr_label(profile[fit_id %in% c("S4", "E4") & !is.na(cohort)])[, outcome := fifelse(fit_id == "S4", "gap start", "gap end")]
  est <- if (level) "pct" else "ratio_to_average_day"; lo <- if (level) "pct_low" else "ratio_low"; hi <- if (level) "pct_high" else "ratio_high"
  d[, cohort := factor(cohort, levels = cfg$cohorts)]
  add_interval_columns(d, est, lo, hi, log = !level, style = style)
  lab <- merge(unique(d[, .(stream, outcome)]), tests[id %in% c("H4s", "H4e"), .(stream = factor(WR_STREAMS[stream], levels = WR_STREAMS), outcome, p_holm)], by = c("stream", "outcome"), all.x = TRUE)
  p <- ggplot(d, aes(weekday, est_, color = cohort, group = cohort)) + WR_SHADE + geom_line(position = position_dodge(0.5), alpha = 0.6) +
    interval_layer(style, horizontal = FALSE, color_var = "cohort", position = position_dodge(0.5)) +
    geom_text(data = lab, aes(x = 1, y = Inf, label = paste("H4 Holm p", wr_ptext(p_holm))), hjust = 0, vjust = 1.2, size = 2.8, inherit.aes = FALSE) +
    facet_grid(outcome ~ stream, scales = "free_y") + scale_color_manual(values = cfg$palette$cohort) + scale_fill_manual(values = cfg$palette$cohort) +
    labs(title = if (level) "Weekly Rhythm by Cohort (Level)" else "Weekly Rhythm by Cohort (Shape)", x = NULL,
         y = if (level) "Predicted % (population average)" else "Ratio to the cohort's own average day", color = "Cohort", fill = "Cohort", caption = interval_caption(style)) + wr_theme()
  if (!level) p <- p + scale_y_log10() + geom_hline(yintercept = 1, linetype = "dashed", color = "grey40")
  p
}

# Figure 5: sensitivity of the weekday ratios; one row per analysis, one narrow column per weekday, colour = stream.
plot_wr_sensitivity <- function(profile, outcome = c("start", "end"), cfg) {
  outcome <- match.arg(outcome); style <- cfg$plot$interval_style %||% "se"
  ids <- if (outcome == "start") c(main = "S1", holidays = "S1_holidays", before_first_gap = "S1_before_first_gap", activity_same_rows = "S1_activity_same_rows", activity = "S1_activity")
         else c(main = "E1", holidays = "E1_holidays")
  d <- wr_label(profile[fit_id %in% ids])[, analysis := factor(names(ids)[match(fit_id, ids)], levels = rev(names(ids)))]
  add_interval_columns(d, "ratio_to_average_day", "ratio_low", "ratio_high", log = TRUE, style = style)
  ggplot(d, aes(est_, analysis, color = stream)) + geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
    interval_layer(style, horizontal = TRUE, color_var = "stream", position = position_dodge(0.6), width = 0.4) +
    facet_grid(. ~ weekday) + scale_x_log10(breaks = c(0.5, 1, 2)) + scale_color_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2") +
    labs(title = sprintf("Weekday Ratios Under Sensitivity Analyses (Gap %s)", if (outcome == "start") "Starts" else "Ends"), x = "Ratio to the average day (log scale)", y = NULL,
         color = "Stream", fill = "Stream", caption = interval_caption(style)) + wr_theme()
}

# Supplement: share of gaps of each length by start weekday (rows add to 100%), one panel per stream.
plot_wr_gap_length <- function(gap_length) {
  d <- wr_label(gap_length)
  ggplot(d, aes(len, weekday, fill = share_pct)) + geom_tile(color = "white") + facet_wrap(~stream, nrow = 1) + scale_y_discrete(limits = rev(WR_WEEKDAYS)) +
    scale_fill_gradient(low = "white", high = "#1B4F72", na.value = "grey80", name = "% of gaps starting that weekday") +
    labs(title = "Gap Length by Start Weekday", x = "Gap length (days)", y = "Start weekday") + wr_theme()
}

# Supplement: days at risk per weekday x days since return cell, one stream.
plot_wr_dsr_counts <- function(cells, s) {
  d <- wr_label(cells[cells$stream == s])
  ggplot(d, aes(weekday, days_since_return, fill = n_days)) + geom_tile(color = "white") + scale_y_discrete(limits = rev(levels(d$days_since_return))) +
    scale_fill_viridis_c(trans = "log10", na.value = "grey80", name = "days at risk") +
    labs(title = sprintf("Days at Risk by Weekday and Days Since Return (%s stream)", WR_STREAMS[[s]]), x = NULL, y = "Days since return") + wr_theme()
}

# Supplement: crude gap-start rate by days since return for each device family (families with enough people), one panel per stream.
plot_wr_dsr_by_device <- function(by_device) {
  d <- wr_label(by_device[!is.na(rate_pct)])
  ggplot(d, aes(days_since_return, rate_pct, group = device, color = device)) + geom_line() + geom_point(size = 1) + facet_wrap(~stream, nrow = 1, scales = "free_y") +
    labs(title = "Gap Starts by Days Since Return and Device Family (crude)", x = "Days since return", y = "% of days followed by a gap", color = "Device family") +
    wr_theme() + theme(axis.text.x = element_text(size = 6, angle = 90))
}


# The figures of the stage as a plan: the tables each needs and how to draw it. A figure whose tables are missing or empty is not drawn.
wr_figure_plan <- function(wr, s, cfg) {
  list(
    weekday_rhythm_profile = list(need = c("profile", "tests"), make = function() plot_wr_profile(wr$profile, wr$summary, wr$tests, wr$crude, cfg)),
    weekday_rhythm_cohort = list(need = "profile_cohort", make = function() plot_wr_cohort(wr$profile, wr$tests, cfg)),
    weekday_rhythm_cohort_level = list(need = "profile_cohort", make = function() plot_wr_cohort(wr$profile, wr$tests, cfg, TRUE)),
    weekday_rhythm_sensitivity_start = list(need = "profile", make = function() plot_wr_sensitivity(wr$profile, "start", cfg)),
    weekday_rhythm_sensitivity_end = list(need = "profile_end", make = function() plot_wr_sensitivity(wr$profile, "end", cfg)),
    weekday_rhythm_gap_length = list(need = "gap_length", make = function() plot_wr_gap_length(wr$gap_length)),
    weekday_rhythm_dsr_by_device = list(need = "dsr_by_device", make = function() plot_wr_dsr_by_device(wr$dsr_by_device)))
}

# The figures of one stream (the same plan form).
wr_stream_figure_plan <- function(wr, s, cfg) {
  stats::setNames(list(
    list(need = c("dsr_cells", "dsr_curve", "shrinkage"), make = function() plot_wr_calendar_vs_cycle(s, wr$dsr_cells, wr$profile, wr$dsr_curve, wr$shrinkage, wr$tests, cfg)),
    list(need = "personal", make = function() plot_wr_personal(s, wr$personal, wr$blup_bins, wr$ellipse, wr$tests)),
    list(need = "dsr_cells", make = function() plot_wr_dsr_counts(wr$dsr_cells, s))),
    paste0(c("weekday_rhythm_calendar_vs_cycle_", "weekday_rhythm_personal_", "weekday_rhythm_dsr_counts_"), s))
}

# Draw every planned figure whose tables exist (and, for per-stream tables, hold rows of the stream).
wr_draw_plan <- function(plan, wr) {
  has <- function(n) { x <- wr[[n]]; !is.null(x) && nrow(x) > 0L }
  ok <- vapply(plan, function(p) all(vapply(p$need, has, logical(1))), logical(1))
  lapply(plan[ok], function(p) p$make())
}

# All figures of the stage that the collected results allow, by output name.
wr_figures <- function(wr, streams, cfg) {
  wr$profile_cohort <- wr$profile[fit_id %in% c("S4", "E4")]; wr$profile_end <- wr$profile[fit_id %in% c("E1", "E1_holidays")]
  per_stream <- lapply(intersect(streams, unique(wr$profile$stream)), function(s) wr_draw_plan(wr_stream_figure_plan(wr, s, cfg), wr))
  c(wr_draw_plan(wr_figure_plan(wr, NULL, cfg), wr), unlist(per_stream, recursive = FALSE))
}
