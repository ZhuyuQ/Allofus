# ==============================================================================
# plots.R
# ------------------------------------------------------------------------------
# All figures are ggplot2 and keep the look of the earlier plotting code:
#   * violin + white boxplot, facet_wrap(scales = "free_y"), theme_bw(base_size = 11),
#     bold strips, 30-degree x labels, legend at the bottom      (value / missingness plots)
#   * theme_minimal(base_size = 13), legend on top, bold centred title (demographic style)
#   * the cohort / sex / age / CCI palettes in cfg$palette, italic "KW p" labels,
#     Dunn tests with Bonferroni adjustment, ggsave(dpi = 150)
# Every function RETURNS a ggplot object; saving is done by save_plot().
# ==============================================================================

suppressPackageStartupMessages(library(ggplot2))

# ---- themes and small helpers -----------------------------------------------------------

theme_violin_panel <- function() {
  theme_bw(base_size = 11) +
    theme(strip.text = element_text(face = "bold"),
          axis.text.x = element_text(angle = 30, hjust = 1),
          legend.position = "bottom")
}

theme_minimal_top <- function() {
  theme_minimal(base_size = 13) +
    theme(legend.position = "top", plot.title = element_text(face = "bold", hjust = 0.5),
          strip.text = element_text(face = "bold"))
}

kw_label <- function(p) paste0("KW p", ifelse(p < 0.001, "<0.001", paste0("=", round(p, 3))))

nice_label <- function(x) {
  out <- tools::toTitleCase(gsub("_+", " ", x))
  gsub("\\bHr\\b", "HR", gsub("\\bCci\\b", "CCI", gsub("\\bMvpa\\b", "MVPA", out)))
}

# Kruskal-Wallis p-value of `value` across `cohort` (NA if fewer than two groups).
kw_p <- function(value, cohort) {
  ok <- !is.na(value)
  if (uniqueN(cohort[ok]) < 2L) return(NA_real_)
  stats::kruskal.test(value[ok] ~ factor(cohort[ok]))$p.value
}

save_plot <- function(p, path, width = 7, height = 5, cfg = make_config()) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggsave(p, filename = path, width = width, height = height, dpi = cfg$plot$dpi)
  invisible(path)
}

# ---- intervals of ESTIMATES: box-and-whisker style ------------------------------------------------------------------
#
# Every figure that shows an estimate with its uncertainty (forest plots, cohort estimates, predicted probabilities,
# transition probabilities) draws ONE BOX PER ESTIMATE:
#     line in the box = the estimate     box = the central 50% interval     whiskers = the 95% interval
# The box is NOT the distribution of the data: it is the uncertainty of the same estimate. It is obtained from the reported
# 95% interval by scaling each side of it with qnorm(0.75) / qnorm(0.975) = 0.344 (exact for a normal sampling distribution,
# on the log scale for odds ratios; an approximation for the percentile intervals of probabilities and bootstrap estimates).
# cfg$plot$interval_style: "box" (default), "se" (one dot and +/- 1 SD = standard error), "ci95" (one dot and the 95% interval).
INTERVAL_K <- stats::qnorm(0.75) / stats::qnorm(0.975)

# Adds the columns est_, lo_, hi_, q25_, q75_ to `d`, in place: the estimate, the ends of the whiskers / bar and the box.
add_interval_columns <- function(d, est, lo, hi, log = FALSE, scale = 1, style = "box") {
  e <- d[[est]]; l <- d[[lo]]; h <- d[[hi]]
  f <- if (log) base::log else identity; g <- if (log) base::exp else identity
  q25 <- g(f(e) - (f(e) - f(l)) * INTERVAL_K); q75 <- g(f(e) + (f(h) - f(e)) * INTERVAL_K)
  if (identical(style, "se")) {                                    # one dot and +/- 1 SD (SE), each side separately (percentile intervals can be asymmetric)
    se_lo <- (f(e) - f(l)) / stats::qnorm(0.975); se_hi <- (f(h) - f(e)) / stats::qnorm(0.975)
    l <- g(f(e) - se_lo); h <- g(f(e) + se_hi)
  }
  set(d, j = "est_", value = scale * e); set(d, j = "lo_", value = scale * l); set(d, j = "hi_", value = scale * h)
  set(d, j = "q25_", value = scale * q25); set(d, j = "q75_", value = scale * q75)
  invisible(d)
}

# The layer that draws what add_interval_columns() prepared.
#   horizontal = TRUE for estimates along x (forest plots), FALSE along y; color_var: a column to colour by (NULL: one colour).
interval_layer <- function(style = "box", horizontal = TRUE, color = "#2C3E50", color_var = NULL, width = 0.5,
                           position = "identity", show.legend = NA) {
  by <- if (!is.null(color_var)) aes(color = .data[[color_var]], fill = .data[[color_var]]) else NULL
  fixed <- if (is.null(color_var)) list(color = color, fill = color) else list()      # an explicit NULL colour would break the mapping
  if (!identical(style, "box")) {
    m <- if (horizontal) aes(xmin = lo_, xmax = hi_) else aes(ymin = lo_, ymax = hi_)
    if (!is.null(by)) m <- utils::modifyList(m, aes(color = .data[[color_var]]))
    pt <- if (is.null(color_var)) list(color = color) else list()
    return(do.call(geom_pointrange, c(list(mapping = m, size = 0.5, position = position, show.legend = show.legend), pt)))
  }
  m <- if (horizontal) aes(xmin = lo_, xlower = q25_, xmiddle = est_, xupper = q75_, xmax = hi_) else
                       aes(ymin = lo_, lower = q25_, middle = est_, upper = q75_, ymax = hi_)
  if (!is.null(by)) m <- utils::modifyList(m, by)
  do.call(geom_boxplot, c(list(mapping = m, stat = "identity", width = width, alpha = 0.3, fatten = 3, linewidth = 0.5,
                               position = position, show.legend = show.legend, orientation = if (horizontal) "y" else "x"), fixed))
}

interval_caption <- function(style = "box") {
  if (identical(style, "ci95")) "Dot: estimate; bar: 95% interval"
  else if (identical(style, "se")) "Dot: estimate; bar: +/- 1 SD of the estimate (standard error)"
  else "Line: estimate; box: 50% interval; whiskers: 95% interval (uncertainty of the estimate, not the spread of the data)"
}

# ---- layout of violin figures ----------------------------------------------------------------------------------------
#
# A figure with many violins needs room to be read, both along x (the groups of one panel) and along y (the spread of
# one violin). With `n_per_facet` violins in the fullest panel and `n_facets` panels, the TOTAL is n_per_facet x n_facets:
#   total <= 12     the panels keep the layout they had (side by side), no size is changed;
#   13 .. 24        the panels are placed in 2 ROWS;     more than 24: in 3 ROWS (never more rows than panels);
# and the image size follows from the layout: every panel gets 0.45 inch per violin (at least 3.5 inches, at least 9 inches
# when a row holds a single panel) and 3.8 inches of height per row. The size is stored in attr(plot, "size") = c(width,
# height) and used when the figure is saved.
violin_layout <- function(n_per_facet, n_facets, default_ncol) {
  total <- n_per_facet * n_facets
  if (total <= 12L || n_facets < 2L) return(list(nrow = NULL, ncol = default_ncol, size = NULL, total = total))
  rows <- min(n_facets, if (total <= 24L) 2L else 3L)
  ncol <- ceiling(n_facets / rows)
  panel_w <- max(3.5, 0.45 * n_per_facet); if (ncol == 1L) panel_w <- max(panel_w, 9)
  list(nrow = rows, ncol = ncol, size = c(min(ncol * panel_w + 1.2, 24), rows * 3.8 + 1.8), total = total)
}

# ---- values by cohort -----------------------------------------------------------------------

# Violin + box of person-level means, one panel per metric, with the KW p-value.
# `long` needs columns: cohort, metric, value.
plot_metric_violin <- function(long, cfg, title, ylab = "Mean per day", ncol = 2L) {
  d <- long[!is.na(value)]
  d[, `:=`(cohort = factor(cohort, levels = cfg$cohorts), metric_label = nice_label(metric))]
  kw <- d[, .(p = kw_p(value, cohort)), by = metric_label][!is.na(p)][, label := kw_label(p)]
  lay <- violin_layout(d[, uniqueN(cohort), by = metric_label][, max(V1)], uniqueN(d$metric_label), ncol)
  p <- ggplot(d, aes(x = cohort, y = value, fill = cohort)) +
    geom_violin(color = NA, alpha = 0.8) +
    geom_boxplot(width = 0.15, outlier.size = 0.3, outlier.alpha = 0.6, fill = "white", color = "gray30") +
    geom_text(data = kw, aes(x = 2.5, y = Inf, label = label), inherit.aes = FALSE,
              vjust = 1.5, size = 3, fontface = "italic") +
    facet_wrap(~metric_label, scales = "free_y", nrow = lay$nrow, ncol = lay$ncol) +
    scale_fill_manual(values = cfg$palette$cohort) +
    labs(title = title, x = NULL, y = ylab, fill = "Cohort") +
    theme_violin_panel()
  attr(p, "size") <- lay$size
  p
}

# One metric, one panel, with Dunn pairwise brackets (Bonferroni). The y-axis is cut at
# Q3 + 3*IQR (as in the steps plot) and the number of people above it is printed.
plot_single_metric_violin <- function(d, yvar, cfg, title, ylab, ref_line = NULL) {
  d <- copy(d)[!is.na(get(yvar))]
  d[, cohort := factor(cohort, levels = cfg$cohorts)]
  d[, y := get(yvar)]
  upper <- max(d[, .(u = stats::quantile(y, 0.75) + 3 * stats::IQR(y)), by = cohort]$u)
  out_n <- d[, .(n_out = sum(y > upper)), by = cohort][n_out > 0L][
    , `:=`(label = paste0("n=", n_out, " above"), y_pos = upper * 0.98)]
  dunn <- rstatix::dunn_test(as.data.frame(d), y ~ cohort, p.adjust.method = "bonferroni")
  dunn <- rstatix::add_xy_position(dunn, x = "cohort")            # x positions of the brackets
  dunn <- dunn[dunn$p.adj < 0.05, ]                               # only significant pairs are drawn
  top <- max(d[y <= upper, y]); step <- 0.08 * (top - min(d$y))
  dunn$y.position <- top + step * seq_len(nrow(dunn))             # heights set from the data, not typed by hand
  # aesthetics live on the layers (not on ggplot()) so the bracket layer inherits none
  p <- ggplot() +
    geom_violin(data = d[y <= upper], aes(x = cohort, y = y, fill = cohort), color = NA, alpha = 0.8) +
    geom_boxplot(data = d, aes(x = cohort, y = y), width = 0.1, outlier.shape = NA, fill = "white", color = "gray30") +
    coord_cartesian(ylim = c(min(d$y), top + step * (nrow(dunn) + 1L))) +
    scale_fill_manual(values = cfg$palette$cohort) +
    labs(title = title, x = NULL, y = ylab, fill = "Cohort") +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "bottom")
  if (nrow(dunn) > 0L) p <- p + ggpubr::stat_pvalue_manual(dunn, label = "p.adj.signif", tip.length = 0.01, size = 3)
  if (!is.null(ref_line)) p <- p + geom_hline(yintercept = ref_line, linetype = "dashed", color = "gray40", linewidth = 0.6)
  if (nrow(out_n) > 0L) {
    p <- p + geom_text(data = out_n, aes(x = cohort, y = y_pos, label = label), inherit.aes = FALSE,
                       size = 3, color = "#C0392B", vjust = 1, fontface = "italic")
  }
  p
}

# Person-level means by cohort within strata (sex, age group, race or CCI group):
# dodged violins coloured by stratum, one panel per metric; KW p-value between
# cohorts printed for each stratum level.
plot_metric_by_stratum <- function(long, stratum, cfg, title, ylab = "Mean per day", ncol = 2L) {
  pal <- switch(stratum, sex = cfg$palette$sex, age_group = cfg$palette$age_group,
                cci_group = cfg$palette$cci, race_plot = cfg$palette$race,
                stop("unknown stratum: ", stratum, call. = FALSE))
  d <- long[!is.na(value) & !is.na(get(stratum))]
  d[, `:=`(cohort = factor(cohort, levels = cfg$cohorts), metric_label = nice_label(metric),
           stratum_f = factor(get(stratum), levels = names(pal)))]
  kw <- d[, .(n_c = uniqueN(cohort), p = kw_p(value, cohort)), by = .(metric_label, stratum_f)][
    n_c > 1L & !is.na(p)]
  kw[, `:=`(label = paste0("KW p", ifelse(p < 0.001, "<0.001", paste0("=", round(p, 3))), " (", stratum_f, ")"),
            vjust = 1.5 + 1.4 * (as.integer(stratum_f) - 1L))]
  dodge <- position_dodge(width = 0.9)
  lay <- violin_layout(d[, uniqueN(paste(cohort, stratum_f)), by = metric_label][, max(V1)], uniqueN(d$metric_label), ncol)
  p <- ggplot(d, aes(x = cohort, y = value, fill = stratum_f)) +
    geom_violin(color = NA, alpha = 0.8, position = dodge) +
    geom_boxplot(width = 0.08, outlier.size = 0.3, outlier.alpha = 0.6, color = "gray30", position = dodge) +
    geom_text(data = kw, aes(x = Inf, y = Inf, label = label, vjust = vjust), hjust = 1.05,
              inherit.aes = FALSE, size = 2.8, fontface = "italic") +
    facet_wrap(~metric_label, scales = "free_y", nrow = lay$nrow, ncol = lay$ncol) +
    scale_fill_manual(values = pal, drop = TRUE) +
    labs(title = title, subtitle = "Kruskal-Wallis test between cohorts within each group",
         x = NULL, y = ylab, fill = nice_label(stratum)) +
    theme_violin_panel()
  attr(p, "size") <- lay$size
  p
}

# ---- missingness ------------------------------------------------------------------------------------

# ECDF of the longest gap per person ("share of people with longest gap <= x days").
plot_gap_ecdf <- function(pg, cfg, max_x = 14L) {
  d <- copy(pg)[, cohort := factor(cohort, levels = cfg$cohorts)]
  ggplot(d, aes(x = pmin(gap_max, max_x), color = cohort)) +
    stat_ecdf(linewidth = 0.8) +
    facet_wrap(~stream, ncol = 1) +
    scale_color_manual(values = cfg$palette$cohort) +
    scale_x_continuous(breaks = c(0, 3, 7, max_x)) +
    labs(title = "ECDF of Longest Continuous Missing Run",
         subtitle = "Y = proportion of people with longest gap <= x days",
         x = "Longest consecutive missing days", y = "Cumulative proportion", color = "Cohort") +
    theme_bw(base_size = 11) + theme(strip.text = element_text(face = "bold"), legend.position = "bottom")
}

# Distribution of gap lengths (share of all gaps in each length group) by cohort.
plot_gap_length_distribution <- function(gaps, cfg) {
  d <- copy(gaps)[, `:=`(cohort = factor(cohort, levels = cfg$cohorts),
                         len = gap_length_bin(gap_length, cfg))]
  s <- d[, .N, by = .(stream, cohort, len)][, share := 100 * N / sum(N), by = .(stream, cohort)]
  ggplot(s, aes(x = len, y = share, fill = cohort)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.65) +
    facet_wrap(~stream) +
    scale_fill_manual(values = cfg$palette$cohort) +
    scale_y_continuous(labels = function(x) paste0(x, "%"), expand = expansion(mult = c(0, 0.08))) +
    labs(title = "Gap Length Distribution", x = "Gap length (days)", y = "% of gaps", fill = "Cohort") +
    theme_minimal_top()
}

# Onset vs continuation probabilities (with bootstrap intervals) by cohort.
plot_transition <- function(tr, cfg) {
  long <- rbindlist(list(
    tr[, .(stream, cohort, measure = "P(missing tomorrow | observed today)", est = p_onset, lo = p_onset_low, hi = p_onset_high)],
    tr[, .(stream, cohort, measure = "P(missing tomorrow | missing today)", est = p_continue, lo = p_continue_low, hi = p_continue_high)]))
  long[, cohort := factor(cohort, levels = cfg$cohorts)]
  style <- cfg$plot$interval_style %||% "box"
  add_interval_columns(long, "est", "lo", "hi", style = style)
  ggplot(long, aes(x = cohort, y = est_, color = cohort, fill = cohort)) +
    interval_layer(style, horizontal = FALSE, color_var = "cohort", show.legend = FALSE) +
    facet_grid(stream ~ measure, scales = "free_y", labeller = labeller(measure = label_wrap_gen(30))) +
    scale_color_manual(values = cfg$palette$cohort) + scale_fill_manual(values = cfg$palette$cohort) +
    labs(title = "Do Missing Days Come in Runs?", x = NULL, y = "Probability", color = "Cohort", fill = "Cohort",
         caption = interval_caption(style)) +
    theme_violin_panel()
}

# One configuration for every "missing share by ..." line figure (study week, calendar year, weekday, age):
# a line per cohort with a shaded band for the 95% bootstrap interval, one panel per stream (panels listed in ROWS,
# one above the other, like every line figure), the cohort palette,
# theme_bw(base_size = 11), bold strips, legend at the bottom. `discrete = TRUE` for a categorical x (weekday).
plot_missing_rate <- function(rates, xvar, cfg, xlab, title, discrete = FALSE, subtitle = NULL) {
  if (nrow(rates) == 0L) {                       # e.g. no age with >= 30 people in any cohort: say so instead of failing
    return(ggplot() + annotate("text", x = 0, y = 0, label = "No group has enough people to be shown\n(at least 30 per point)", size = 4) +
             labs(title = title) + theme_void(base_size = 11))
  }
  d <- copy(rates)[, cohort := factor(cohort, levels = cfg$cohorts)]
  if (discrete) d[, (xvar) := factor(get(xvar), levels = DOW_LABELS)]
  ggplot(d, aes(x = .data[[xvar]], y = 100 * estimate, color = cohort, fill = cohort, group = cohort)) +
    geom_ribbon(aes(ymin = 100 * lo, ymax = 100 * hi), alpha = 0.15, color = NA) +
    geom_line(linewidth = 0.8) +
    facet_wrap(~stream, ncol = 1, scales = "free_y") +
    scale_color_manual(values = cfg$palette$cohort) + scale_fill_manual(values = cfg$palette$cohort) +
    labs(title = title, subtitle = subtitle, x = xlab, y = "% of days missing", color = "Cohort", fill = "Cohort") +
    theme_bw(base_size = 11) + theme(strip.text = element_text(face = "bold"), legend.position = "bottom")
}

# Missing share over a time axis (study week or calendar year).
plot_missing_over_time <- function(rates, xvar, cfg, xlab, title)
  plot_missing_rate(rates, xvar, cfg, xlab, title)

# Missing share by weekday and cohort.
plot_missing_by_weekday <- function(rates, cfg)
  plot_missing_rate(rates, "dow", cfg, NULL, "Missing Days by Day of Week", discrete = TRUE)

# Missing share by attained age (years, continuous), per cohort.
plot_missing_by_age <- function(rates, cfg)
  plot_missing_rate(rates, "age", cfg, "Age (years)", "Missing Days by Age",
                    subtitle = "5-year moving window; bands: 95% person-level bootstrap interval")

# Plain names of the strata for titles and legends.
stratum_label <- function(stratum) c(sex = "Sex", age_group = "Age Group", race_plot = "Race", cci_group = "CCI Group")[[stratum]]

# Person-level % missing by cohort within the levels of one stratum (sex, age_group, race_plot, cci_group): the dodged
# violin + box layout of plot_metric_by_stratum(), one panel per stream. Above each cohort the test between the
# levels of the stratum within that cohort (Kruskal-Wallis for >= 3 levels, Wilcoxon for 2) with its Holm-adjusted p value
# (family "groups_within_cohort" of missing_pct_stratum_tests()).
plot_missing_pct_by_stratum <- function(long, stratum, tests, cfg, title) {
  pal <- switch(stratum, sex = cfg$palette$sex, age_group = cfg$palette$age_group,
                cci_group = cfg$palette$cci, race_plot = cfg$palette$race,
                stop("unknown stratum: ", stratum, call. = FALSE))
  d <- long[!is.na(value) & !is.na(get(stratum))]
  d[, `:=`(cohort = factor(cohort, levels = cfg$cohorts), metric_label = nice_label(stream),
           stratum_f = factor(get(stratum), levels = names(pal)))]
  t <- tests[stratum_var == stratum & family == "groups_within_cohort" & test != "none"]
  t[, `:=`(cohort = factor(group, levels = cfg$cohorts), metric_label = nice_label(stream),
           label = paste0(fifelse(test == "Kruskal-Wallis", "KW", "W"), " p",
                          fifelse(p_holm < 0.001, "<0.001", paste0("=", signif(p_holm, 2)))))]
  dodge <- position_dodge(width = 0.9)
  lay <- violin_layout(d[, uniqueN(paste(cohort, stratum_f)), by = metric_label][, max(V1)], uniqueN(d$metric_label), 3L)
  p <- ggplot(d, aes(x = cohort, y = value, fill = stratum_f)) +
    geom_violin(color = NA, alpha = 0.8, position = dodge) +
    geom_boxplot(width = 0.08, outlier.size = 0.3, outlier.alpha = 0.6, color = "gray30", position = dodge) +
    geom_text(data = t, aes(x = cohort, y = Inf, label = label), inherit.aes = FALSE, vjust = 1.4, size = 2.8, fontface = "italic") +
    facet_wrap(~metric_label, scales = "free_y", nrow = lay$nrow, ncol = lay$ncol) +
    scale_fill_manual(values = pal, drop = TRUE) +
    scale_x_discrete(drop = FALSE) +
    labs(title = title,
         subtitle = "Test between the groups within each cohort (KW = Kruskal-Wallis; W = Wilcoxon rank-sum, two groups)\np Holm-adjusted over cohorts",
         x = NULL, y = "% of days missing (one value per person)", fill = stratum_label(stratum)) +
    theme_violin_panel()
  attr(p, "size") <- lay$size
  p
}

# Observed (black) vs missing (white) days for a seeded sample of people - the only
# place a calendar grid is built, and only for the sample.
plot_missing_heatmap <- function(days, cfg, n_people = 300L, title = "Observed vs Missing Days") {
  ids <- sample_person_ids(days$person_id, n_people, cfg$seed)
  d <- days[person_id %in% ids, .(person_id, date)]
  w <- d[, .(first_date = min(date), last_date = max(date)), by = person_id]
  w[, window_days := as.integer(last_date - first_date) + 1L]
  grid <- w[, .(study_day = seq_len(window_days)), by = .(person_id, first_date, window_days)]
  grid[, date := first_date + study_day - 1L]
  grid[, observed := d[grid, on = c("person_id", "date"), .N, by = .EACHI]$N > 0L]
  grid[, person_id := factor(person_id, levels = w[order(-window_days), person_id])]
  ggplot(grid, aes(x = study_day, y = person_id, fill = observed)) +
    geom_tile(linewidth = 0) +
    scale_fill_manual(values = c(`FALSE` = "white", `TRUE` = "black"),
                      labels = c(`FALSE` = "Missing", `TRUE` = "Observed"), name = NULL) +
    labs(title = title, x = "Study day", y = "Person") +
    theme_minimal(base_size = 11) +
    theme(panel.grid = element_blank(), axis.text.y = element_blank(), legend.position = "bottom",
          legend.background = element_rect(color = "black", linewidth = 0.5))
}

# ---- model results ---------------------------------------------------------------------------------------

# Forest plot of odds ratios. `tab` needs: label, or, or_low, or_high, and optionally group.
plot_forest <- function(tab, title, xlab = "Odds ratio (log scale)", group_col = NULL, style = "box") {
  d <- copy(tab)[!is.na(or)]
  d[, label := factor(label, levels = rev(unique(label)))]
  add_interval_columns(d, "or", "or_low", "or_high", log = TRUE, style = style)
  p <- ggplot(d, aes(x = est_, y = label)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
    interval_layer(style, horizontal = TRUE) +
    scale_x_log10(breaks = c(0.5, 0.8, 1, 1.25, 1.5, 2, 3, 5)) + expand_limits(x = 1) +
    labs(title = title, x = xlab, y = NULL, caption = interval_caption(style)) +
    theme_bw(base_size = 11) + theme(strip.text = element_text(face = "bold"))
  if (!is.null(group_col)) p <- p + facet_wrap(stats::as.formula(paste("~", group_col)), scales = "free_y")
  p
}

# Predicted next-day gap risk as the prior day gets less/more active than usual. Population-average
# risk (integrated over the person random intercept) from the unadjusted joint model with all other
# predictors at their average; the bars are the 95% interval from the uncertainty of the intercept
# and the within-person slope only (the person variance is held fixed).
plot_risk_curve <- function(risk, cfg, title = "Gap Risk the Next Day by Yesterday's Activity") {
  d <- copy(risk)[!is.na(risk)]
  p <- ggplot(d, aes(x = within_z, y = 100 * risk)) +
    geom_line(color = "#2C3E50", linewidth = 0.8) + geom_point(color = "#2C3E50", size = 2) +
    facet_wrap(~stream, ncol = 1, scales = "free_y") +
    scale_x_continuous(breaks = c(-1, 0, 1), labels = c("-1 SD", "usual", "+1 SD")) +
    labs(title = title, x = "Yesterday vs the person's usual (SD)", y = "% chance a gap starts tomorrow") +
    theme_bw(base_size = 11) + theme(strip.text = element_text(face = "bold"))
  if (all(c("risk_low", "risk_high") %in% names(d)) && any(!is.na(d$risk_low)))
    p <- p + geom_errorbar(aes(ymin = 100 * risk_low, ymax = 100 * risk_high), width = 0.12, color = "#2C3E50")
  p
}

# Within-person odds ratios of the primary (joint, adjusted) onset model: one panel per stream, one
# row per prior-day predictor. OR > 1: a lower-than-usual day goes with higher odds that a gap starts tomorrow.
plot_onset_forest <- function(coefs, cfg, title = "Prior-Day Activity and the Odds That a Gap Starts Tomorrow") {
  d <- coefs[role == "primary" & grepl("^w_", term) & !is.na(or)]
  d[, label := paste0(nice_label(sub("^w_", "", term)))]
  d[, panel := paste0(nice_label(stream), "  (", format(n_events[1L], big.mark = ","), " gaps, ", format(n_persons[1L], big.mark = ","), " people)"), by = stream]
  style <- cfg$plot$interval_style %||% "box"
  add_interval_columns(d, "or", "or_low", "or_high", log = TRUE, style = style)
  ggplot(d, aes(x = est_, y = label)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
    interval_layer(style, horizontal = TRUE) +
    facet_wrap(~panel, scales = "free_y", ncol = 1) +
    scale_x_log10(breaks = c(0.5, 0.8, 1, 1.25, 1.5, 2, 3)) + expand_limits(x = 1) +
    labs(title = title, x = "Odds ratio per 1 SD below the person's usual (log scale)", y = NULL, caption = interval_caption(style)) +
    theme_bw(base_size = 11) + theme(strip.text = element_text(face = "bold"))
}

# The primary within-person effect by cohort (separate fits) with the cohort-interaction test.
plot_cohort_forest <- function(cohort_specific, cohort_interaction, cfg,
                               title = "Primary Effect by Cohort") {
  d <- cohort_specific[!is.na(or)]
  d[, cohort := factor(cohort_stratum, levels = rev(cfg$cohorts))]
  all_pairs <- CJ(stream = unique(cohort_specific$stream), cohort_stratum = cfg$cohorts)
  miss <- all_pairs[!d, on = c("stream", "cohort_stratum")]
  sub <- if (nrow(miss)) paste0("Not estimable (too few events): ", paste0(miss$cohort_stratum, " [", miss$stream, "]", collapse = ", "))
  lab <- if (nrow(cohort_interaction) > 0L)
    cohort_interaction[!is.na(p), .(stream, label = paste0("Do cohorts differ? LRT chisq(", df, ")=", round(chisq, 1), ", p", ifelse(p < 0.001, "<0.001", paste0("=", signif(p, 2)))))] else data.table(stream = character(), label = character())
  style <- cfg$plot$interval_style %||% "box"
  add_interval_columns(d, "or", "or_low", "or_high", log = TRUE, style = style)
  ggplot(d, aes(x = est_, y = cohort, color = cohort, fill = cohort)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
    interval_layer(style, horizontal = TRUE, color_var = "cohort", show.legend = FALSE) +
    geom_text(data = lab, aes(x = Inf, y = Inf, label = label), inherit.aes = FALSE, hjust = 1.05, vjust = 1.5,
              size = 3, fontface = "italic") +
    facet_wrap(~stream) +
    scale_color_manual(values = cfg$palette$cohort) + scale_fill_manual(values = cfg$palette$cohort) +
    scale_x_log10(breaks = c(0.5, 0.8, 1, 1.25, 1.5, 2, 3)) + expand_limits(x = 1) +
    labs(title = title, subtitle = sub, x = "Odds ratio per 1 SD below usual steps (log scale)", y = NULL,
         caption = interval_caption(style)) +
    theme_bw(base_size = 11) + theme(strip.text = element_text(face = "bold"))
}

# Model B: estimated mean activity on the first full day back by length of the gap (lines, 95% CI)
# against the observed group means (points, 95% CI that respects clustering within people).
plot_gap_dose_response <- function(curve, observed, cfg,
                                   title = "Activity After a Gap, by How Long the Gap Was") {
  lab <- c(B_back2_unadj = "Unadjusted", B_back2_adj = "Adjusted (cohort, age, sex)")
  cv <- curve[model_id %in% names(lab)][, model := factor(lab[model_id], levels = lab)]
  ob <- observed[model_id == "B_back2_unadj"][, n_lab := paste0("n=", format(n_gaps, big.mark = ","))]
  ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
    geom_ribbon(data = cv, aes(x = gap_length, ymin = lo, ymax = hi, fill = model), alpha = 0.18) +
    geom_line(data = cv, aes(x = gap_length, y = est, color = model), linewidth = 0.9) +
    geom_pointrange(data = ob, aes(x = x, y = mean, ymin = lo, ymax = hi), color = "gray20", size = 0.35) +
    geom_text(data = ob, aes(x = x, y = -Inf, label = n_lab), vjust = -0.6, size = 2.6, color = "gray30") +
    facet_wrap(~stream, ncol = 1, scales = "free_y") +
    scale_x_continuous(trans = "log2", breaks = c(1, 2, 4, 8, 14)) +
    scale_color_manual(values = c("#2C3E50", "#E74C3C")) + scale_fill_manual(values = c("#2C3E50", "#E74C3C")) +
    labs(title = title,
         x = "Gap length (days, log scale)",
         y = "Activity on first full day back\n(SD from the person's usual)", color = NULL, fill = NULL) +
    theme_bw(base_size = 11) + theme(strip.text = element_text(face = "bold"), legend.position = "bottom")
}

# Model C: predicted probability of each gap-length group when the day before the gap was 1 SD below
# usual / usual / 1 SD above usual (population-average, 95% CI).
plot_ordinal_probs <- function(pred, cfg, title = "Gap Length Given Yesterday's Activity") {
  d <- copy(pred)[!is.na(prob)]
  d <- d[, .SD[adjusted == max(adjusted)], by = stream]          # the adjusted fit where there is one
  d[, z := factor(within_z, levels = c(-1, 0, 1), labels = c("1 SD below usual", "usual", "1 SD above usual"))]
  style <- cfg$plot$interval_style %||% "box"
  add_interval_columns(d, "prob", "prob_low", "prob_high", scale = 100, style = style)
  ggplot(d, aes(x = length_group, y = est_, color = z, fill = z)) +
    interval_layer(style, horizontal = FALSE, color_var = "z", width = 0.7, position = position_dodge(width = 0.8)) +
    facet_wrap(~stream, scales = "free_y") +
    scale_color_manual(values = c("#E74C3C", "#7F8C8D", "#3498DB")) + scale_fill_manual(values = c("#E74C3C", "#7F8C8D", "#3498DB")) +
    labs(title = title, x = "Gap length (days)", y = "% of gaps in this length group", color = "Day before the gap",
         fill = "Day before the gap", caption = interval_caption(style)) +
    theme_bw(base_size = 11) + theme(strip.text = element_text(face = "bold"), legend.position = "bottom")
}

# ---- what each figure shows (written next to the figures) ----------------------------------------------

# Captions by figure name; names ending in "_*" match every figure that starts with that prefix.
figure_captions <- function() c(
  missing_pct_by_cohort = "Share of days without a record, one value per person (violin + box). KW p = Kruskal-Wallis test between cohorts. Descriptive.",
  longest_gap_by_cohort = "Longest run of consecutive missing days per person. KW p = Kruskal-Wallis between cohorts. Descriptive.",
  gap_ecdf = "Cumulative share of people whose longest gap is at most x days (x cut at 14). Curves further to the right = longer gaps.",
  gap_length_distribution = "Share of each cohort's gaps in each length group (1, 2, 3, 4-7, 8-14, 15+ days). Gaps are strongly right-skewed, hence grouped.",
  "heatmap_*" = "Observed (black) vs missing (white) days for a seeded sample of 300 people, longest follow-up first.",
  transition = "P(missing tomorrow | observed today) and P(missing tomorrow | missing today) by cohort with person-level bootstrap 95% intervals. A large gap between the two means missing days come in runs.",
  missing_by_study_week = "Share of days missing by week since each person's first record; line = pooled share, shaded band = 95% person-level bootstrap interval. Descriptive; not adjusted.",
  missing_by_year = "Share of days missing by calendar year; line = pooled share, shaded band = 95% person-level bootstrap interval. Descriptive; not adjusted.",
  missing_by_weekday = "Share of days missing by weekday and cohort; line = pooled share, shaded band = 95% person-level bootstrap interval. Formal tests are in weekday_tests.tsv.",
  missing_by_age = "Share of days missing by attained age (years, continuous) and cohort; line = pooled share, shaded band = 95% person-level bootstrap interval; each age pools the days at age -2 to +2 years (a 5-year moving window); age at each day = age at first record + elapsed time (about +/- 0.5 year). Descriptive; not adjusted. Only ages with >= 30 people are drawn.",
  missing_pct_by_sex = "Person-level % of days missing by cohort and sex (violin + box). Label: test between the groups within each cohort (Kruskal-Wallis, Wilcoxon for two groups), Holm-adjusted over cohorts; effect sizes (epsilon-squared / Cliff's delta) are in missing_pct_by_stratum_tests.tsv. Descriptive.",
  missing_pct_by_age_group = "Person-level % of days missing by cohort and age group (violin + box). Label: Kruskal-Wallis between age groups within each cohort, Holm-adjusted over cohorts; effect sizes in missing_pct_by_stratum_tests.tsv. Descriptive; age is at the first record.",
  missing_pct_by_race_plot = "Person-level % of days missing by cohort and race group (violin + box). Label: Kruskal-Wallis between race groups within each cohort, Holm-adjusted over cohorts; groups with fewer than 5 people in a cohort are left out of the test. Descriptive.",
  missing_pct_by_cci_group = "Person-level % of days missing by cohort and Charlson comorbidity group (violin + box). Label: Kruskal-Wallis between CCI groups within each cohort, Holm-adjusted over cohorts. Descriptive; CCI is not used as an adjuster elsewhere.",
  "metrics_*" = "Person-level mean of each metric by cohort (or by cohort within a subgroup). Descriptive; KW p = Kruskal-Wallis between cohorts.",
  onset_primary_forest = "Primary gap-onset model (gap >= 1 day, joint, adjusted). Odds ratio that a gap starts tomorrow per 1 SD BELOW the person's own usual value yesterday (OR > 1: lower-than-usual days precede gaps). Odds ratios, not risk ratios; conditional on the person's random intercept. Association, not causation.",
  forest_gap_definition = "The primary within-person steps effect when a gap is defined as >= 1, 2, 3 or 7 missing days (model D). Same odds-ratio reading as above.",
  risk_curve = "Population-average chance that a gap starts tomorrow when yesterday was 1 SD below / at / 1 SD above the person's usual steps (unadjusted joint model; bars = 95% interval from the intercept and slope only). Absolute risks for illustration; they depend on the random-intercept variance.",
  cohort_forest = "The primary effect fitted separately in each cohort. The label is the likelihood-ratio test (LRT) of the whole model with vs without a cohort x effect interaction: chisq with df = cohorts - 1, small p = the effect of prior-day steps differs between cohorts. It is one overall test of difference, not a comparison of single pairs, and boxes that overlap or do not overlap are not a substitute for it. Separate-cohort estimates are conditional on that cohort's own random intercept.",
  forest_sensitivity = "The primary effect under different assumptions (trailing 14-day baseline, no partial-wear days, first gap only, main-sleep-only, conditional logit, cluster-robust SE). Spline and random-slope analyses are likelihood-ratio tests and are listed in sensitivity.tsv, not drawn here.",
  gap_dose_response = "Model B (exploratory). Mean activity on the first full day after a gap (SD from the person's usual) by gap length, with the fitted straight line in log2(gap length) and its 95% CI; points are observed group means with person-clustered 95% CI. Conditional on a gap having occurred.",
  gap_length_probs = "Model C (exploratory, shown only where proportional odds held). Population-average probability of each gap-length group when the day before the gap was 1 SD below / at / 1 SD above usual steps, with 95% intervals. Not shown for streams where proportional odds was rejected (see model_C_gap_length.tsv)."
)

# figures drawn as box-and-whisker intervals of an estimate (see interval_layer)
INTERVAL_FIGURES <- c("transition", "onset_primary_forest", "forest_gap_definition", "cohort_forest", "forest_sensitivity", "gap_length_probs")
INTERVAL_FIGURE_NOTE <- " One box per estimate: line = estimate, box = central 50% interval, whiskers = 95% interval - the uncertainty of the estimate, not the spread of the data (the 50% interval is derived from the 95% interval, exact for a normal sampling distribution). cfg$plot$interval_style = \"se\" (dot +/- 1 SD) or \"ci95\" (dot + 95% interval) draws the alternatives."

caption_for <- function(name) {
  cp <- figure_captions()
  if (name %in% INTERVAL_FIGURES && name %in% names(cp)) return(paste0(cp[[name]], INTERVAL_FIGURE_NOTE))
  if (name %in% names(cp)) return(cp[[name]])
  pre <- names(cp)[endsWith(names(cp), "_*")]
  hit <- pre[startsWith(name, sub("\\*$", "", pre))]
  if (length(hit)) cp[[hit[1L]]] else NA_character_
}
