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
# Every figure that shows an estimate with its uncertainty can be drawn in three styles (cfg$plot$interval_style):
#   "se"   (default for the MODEL figures) one dot and +/- 1 SD of the estimate (its standard error)
#   "ci95" one dot and the 95% interval
#   "box"  ONE BOX PER ESTIMATE:
#     line in the box = the estimate     box = the central 50% interval     whiskers = the 95% interval
# The box is NOT the distribution of the data: it is the uncertainty of the same estimate. It is obtained from the reported
# 95% interval by scaling each side of it with qnorm(0.75) / qnorm(0.975) = 0.344 (exact for a normal sampling distribution,
# on the log scale for odds ratios; an approximation for the percentile intervals of probabilities and bootstrap estimates).
# The cohort transition probabilities have their own setting, cfg$plot$transition_style (default "box").
INTERVAL_K <- stats::qnorm(0.75) / stats::qnorm(0.975)

# Adds the columns est_, lo_, hi_, q25_, q75_ to `d`, in place: the estimate, the ends of the whiskers / bar and the box.
# Rows whose interval can be drawn: estimate and both limits finite (and positive on a log axis), and a standard error of the
# log odds ratio below `max_se` (an interval of 1e17-fold width is a term the data cannot estimate, e.g. a year with almost no
# events, not a result). Those rows are left out of the figure; the tables keep them.
drawable_interval <- function(d, est, lo, hi, log = FALSE, max_se = 10) {
  e <- d[[est]]; l <- d[[lo]]; h <- d[[hi]]
  ok <- is.finite(e) & is.finite(l) & is.finite(h)
  if (log) ok <- ok & e > 0 & l > 0 & h > 0 & (base::log(pmax(h, 1e-300)) - base::log(pmax(l, 1e-300))) / (2 * stats::qnorm(0.975)) < max_se
  else ok <- ok & (h - l) / (2 * stats::qnorm(0.975)) < max_se
  ok
}

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
  style <- cfg$plot$transition_style %||% "box"
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

# Every effect of a "within-person" predictor is expressed in the NORMAL direction: per +1 SD HIGHER than the person's usual
# (a higher value is a positive step; OR = exp(coefficient)). OR < 1: higher values go with lower odds of the outcome.
# Tick marks of every odds-ratio axis (log scale); the ones outside the data range are not drawn.
OR_BREAKS <- c(0.1, 0.25, 0.5, 0.67, 0.8, 1, 1.25, 1.5, 2, 4, 10)

# Forest plot of odds ratios. `tab` needs: label, or, or_low, or_high, and optionally group.
plot_forest <- function(tab, title, xlab = "Odds ratio per +1 SD (log scale)", group_col = NULL, style = "se") {
  d <- copy(tab)[!is.na(or)]
  d <- d[drawable_interval(d, "or", "or_low", "or_high", log = TRUE)]
  d[, label := factor(label, levels = rev(unique(label)))]
  add_interval_columns(d, "or", "or_low", "or_high", log = TRUE, style = style)
  p <- ggplot(d, aes(x = est_, y = label)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
    interval_layer(style, horizontal = TRUE) +
    scale_x_log10(breaks = OR_BREAKS) + expand_limits(x = 1) +
    labs(title = title, x = xlab, y = NULL, caption = interval_caption(style)) +
    theme_bw(base_size = 11) + theme(strip.text = element_text(face = "bold"))
  if (!is.null(group_col)) p <- p + facet_wrap(stats::as.formula(paste("~", group_col)), scales = "free_y")
  p
}

# Predicted next-day gap risk as the prior day gets less (-1 SD) or more (+1 SD) active than usual. Population-average
# risk (integrated over the person random intercept) from the unadjusted joint model with all other
# predictors at their average; the bars are the 95% interval from the uncertainty of the intercept
# and the within-person slope only (the person variance is held fixed).
plot_risk_curve <- function(risk, cfg, title = "Gap Risk the Next Day by Yesterday's Activity") {
  d <- copy(risk)[!is.na(risk)]
  p <- ggplot(d, aes(x = within_z, y = 100 * risk)) +
    geom_line(color = "#2C3E50", linewidth = 0.8) + geom_point(color = "#2C3E50", size = 2) +
    facet_wrap(~stream, ncol = 1, scales = "free_y") +
    scale_x_continuous(breaks = c(-1, 0, 1), labels = c("-1 SD", "0 (usual)", "+1 SD")) +
    labs(title = title, x = "Yesterday's steps vs usual (SD)", y = "% chance a gap starts tomorrow") +
    theme_bw(base_size = 11) + theme(strip.text = element_text(face = "bold"))
  if (all(c("risk_low", "risk_high") %in% names(d)) && any(!is.na(d$risk_low)))
    p <- p + geom_errorbar(aes(ymin = 100 * risk_low, ymax = 100 * risk_high), width = 0.12, color = "#2C3E50")
  p
}

# Within-person odds ratios of the primary (joint, adjusted) onset model: one panel per stream, one
# row per prior-day predictor. OR < 1: a higher-than-usual day goes with lower odds that a gap starts tomorrow.
plot_onset_forest <- function(coefs, cfg, title = "Prior-Day Activity and the Odds That a Gap Starts Tomorrow") {
  d <- coefs[role == "primary" & grepl("^w_", term) & !is.na(or)]
  d <- d[drawable_interval(d, "or", "or_low", "or_high", log = TRUE)]
  d[, label := paste0(nice_label(sub("^w_", "", term)))]
  d[, panel := paste0(nice_label(stream), "  (", format(n_events[1L], big.mark = ","), " gaps, ", format(n_persons[1L], big.mark = ","), " people)"), by = stream]
  style <- cfg$plot$interval_style %||% "se"
  add_interval_columns(d, "or", "or_low", "or_high", log = TRUE, style = style)
  ggplot(d, aes(x = est_, y = label)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
    interval_layer(style, horizontal = TRUE) +
    facet_wrap(~panel, scales = "free_y", ncol = 1) +
    scale_x_log10(breaks = OR_BREAKS) + expand_limits(x = 1) +
    labs(title = title, x = "Odds ratio per +1 SD (log scale)", y = NULL, caption = interval_caption(style)) +
    theme_bw(base_size = 11) + theme(strip.text = element_text(face = "bold"))
}

# The primary within-person effect by cohort (separate fits) with the cohort-interaction test.
plot_cohort_forest <- function(cohort_specific, cohort_interaction, cfg,
                               title = "Primary Effect by Cohort") {
  d <- cohort_specific[!is.na(or)]
  d <- d[drawable_interval(d, "or", "or_low", "or_high", log = TRUE)]
  d[, cohort := factor(cohort_stratum, levels = rev(cfg$cohorts))]
  all_pairs <- CJ(stream = unique(cohort_specific$stream), cohort_stratum = cfg$cohorts)
  miss <- all_pairs[!d, on = c("stream", "cohort_stratum")]
  sub <- if (nrow(miss)) paste0("Not estimable (too few events): ", paste0(miss$cohort_stratum, " [", miss$stream, "]", collapse = ", "))
  lab <- if (nrow(cohort_interaction) > 0L)
    cohort_interaction[!is.na(p), .(stream, label = paste0("Do cohorts differ? LRT chisq(", df, ")=", round(chisq, 1), ", p", ifelse(p < 0.001, "<0.001", paste0("=", signif(p, 2)))))] else data.table(stream = character(), label = character())
  style <- cfg$plot$interval_style %||% "se"
  add_interval_columns(d, "or", "or_low", "or_high", log = TRUE, style = style)
  ggplot(d, aes(x = est_, y = cohort, color = cohort, fill = cohort)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
    interval_layer(style, horizontal = TRUE, color_var = "cohort", show.legend = FALSE) +
    geom_text(data = lab, aes(x = Inf, y = Inf, label = label), inherit.aes = FALSE, hjust = 1.05, vjust = 1.5,
              size = 3, fontface = "italic") +
    facet_wrap(~stream) +
    scale_color_manual(values = cfg$palette$cohort) + scale_fill_manual(values = cfg$palette$cohort) +
    scale_x_log10(breaks = OR_BREAKS) + expand_limits(x = 1) +
    labs(title = title, subtitle = sub, x = "Odds ratio per +1 SD steps (log scale)", y = NULL,
         caption = interval_caption(style)) +
    theme_bw(base_size = 11) + theme(strip.text = element_text(face = "bold"))
}

# Model B: estimated mean activity on the first full day back at a few gap lengths (cfg$plot$gap_days: 1, 2, 3, 5 and 13 days),
# one forest row per gap length with two colours - the unadjusted and the adjusted fit - and the interval of the estimate
# (cfg$plot$interval_style). The values are the fitted straight line in log2(gap length) of each specification, so they are
# the model's estimates at those lengths (0 = the person's usual activity).
plot_gap_dose_response <- function(curve, observed = NULL, cfg,
                                   title = "Activity After a Gap, by How Long the Gap Was") {
  lab <- c(B_back2_unadj = "Unadjusted", B_back2_adj = "Adjusted")
  days <- cfg$plot$gap_days %||% c(1L, 2L, 3L, 5L, 13L)
  cv <- curve[model_id %in% names(lab) & gap_length %in% days][, model := factor(lab[model_id], levels = lab)]
  cv[, day := factor(gap_length, levels = rev(sort(unique(gap_length))))]
  style <- cfg$plot$interval_style %||% "se"
  add_interval_columns(cv, "est", "lo", "hi", style = style)
  cols <- c("#2C3E50", "#E74C3C")
  ggplot(cv, aes(x = est_, y = day, color = model, fill = model)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
    interval_layer(style, horizontal = TRUE, color_var = "model", width = 0.6, position = position_dodge(width = 0.7)) +
    facet_wrap(~stream, ncol = 1, scales = "free_x") +
    scale_color_manual(values = cols) + scale_fill_manual(values = cols) +
    guides(color = guide_legend(reverse = TRUE), fill = guide_legend(reverse = TRUE)) +
    labs(title = title,
         x = "Activity on first full day back (SD); 0 = usual", y = "Gap length (days)",
         color = NULL, fill = NULL, caption = interval_caption(style)) +
    theme_bw(base_size = 11) + theme(strip.text = element_text(face = "bold"), legend.position = "bottom")
}

# Model C: predicted probability of each gap-length group when the day before the gap was at -1 SD /
# 0 (usual) / +1 SD of the person's usual steps (population-average, with the interval of the estimate).
plot_ordinal_probs <- function(pred, cfg, title = "Gap Length Given Yesterday's Activity") {
  d <- copy(pred)[!is.na(prob)]
  d <- d[, .SD[adjusted == max(adjusted)], by = stream]          # the adjusted fit where there is one
  d[, z := factor(within_z, levels = c(-1, 0, 1), labels = c("-1 SD", "0 (usual)", "+1 SD"))]
  style <- cfg$plot$interval_style %||% "se"
  add_interval_columns(d, "prob", "prob_low", "prob_high", scale = 100, style = style)
  ggplot(d, aes(x = length_group, y = est_, color = z, fill = z)) +
    interval_layer(style, horizontal = FALSE, color_var = "z", width = 0.7, position = position_dodge(width = 0.8)) +
    facet_wrap(~stream, scales = "free_y") +
    scale_color_manual(values = c("#E74C3C", "#7F8C8D", "#3498DB")) + scale_fill_manual(values = c("#E74C3C", "#7F8C8D", "#3498DB")) +
    labs(title = title, x = "Gap length (days)", y = "% of gaps in this length group", color = "Steps before the gap vs usual",
         fill = "Steps before the gap vs usual", caption = interval_caption(style)) +
    theme_bw(base_size = 11) + theme(strip.text = element_text(face = "bold"), legend.position = "bottom")
}

# ---- one forest plot with ALL predictors per model and stream ---------------------------------------------------------------

# Plain-language label of a predictor (the name after w_ / b_), with the transform that was applied.
PREDICTOR_LABELS <- c(act_steps = "Steps", act_sedentary = "Sedentary minutes", act_mvpa = "MVPA minutes",
                      sleep_inbed = "Time in bed", sleep_eff = "Sleep efficiency",
                      hr_cover = "Heart-rate wear minutes", hr_elev = "Elevated heart-rate minutes")
predictor_label <- function(x) ifelse(x %in% names(PREDICTOR_LABELS), PREDICTOR_LABELS[x], nice_label(x))

# Short label of a model term for the y axis of the forest plots (units and reference levels are explained in the subtitle / caption).
term_label <- function(term) {
  fixed <- c(log2_gap = "Gap length", w_pre_day_steps = "Steps before the gap", w_pre_day_mvpa = "MVPA before the gap", age10 = "Age",
             prev_gap_len_l = "Previous gap length", any_prior_gap = "Any earlier gap", days_since_return_l = "Days since return",
             cci = "CCI score", cci_missing = "CCI missing")
  out <- vapply(term, function(t) {
    if (t %in% names(fixed)) fixed[[t]]
    else if (grepl("^t?w_", t)) paste0("Within: ", predictor_label(sub("^t?w_", "", t)))
    else if (grepl("^b_", t)) paste0("Between: ", predictor_label(sub("^b_", "", t)))
    else if (grepl("^cohort", t)) paste0("Cohort: ", sub("^cohort", "", t))
    else if (grepl("^sex2", t)) paste0("Sex: ", sub("^sex2", "", t))
    else if (grepl("^dow_t", t)) paste0("Weekday: ", sub("^dow_t", "", t))
    else if (grepl("^device", t)) paste0("Device: ", sub("^device", "", t))
    else if (grepl("^cci_grp", t)) paste0("CCI group: ", sub("^cci_grp", "", t))
    else if (grepl("^any_prior_gap", t)) "Any earlier gap"
    else t
  }, character(1))
  unname(out)
}

# Which block of the forest plot a term belongs to.
term_block <- function(term) {
  lv <- c("Within-person (per +1 SD above usual)", "Between-person (per +1 SD)", "Covariates and gap history")
  factor(ifelse(grepl("^t?w_", term), lv[1L], ifelse(grepl("^b_", term), lv[2L], lv[3L])), levels = lv)
}

forest_theme <- function() theme_bw(base_size = 11) + theme(strip.text = element_text(face = "bold"), legend.position = "bottom")

# Height (inches) of a forest figure with `n_rows` estimates and `n_panels` panels.
forest_height <- function(n_rows, n_panels = 1L) max(3.5, 1.4 + 0.27 * n_rows + 0.45 * n_panels)

# Gap-onset model of ONE stream (the primary, joint, adjusted model): every fixed effect: within-person terms, between-person terms, covariates
# and gap history (in that order). The intercept and the spline terms of study time are left out (their coefficients have no reading one by one).
plot_onset_forest_all <- function(coefs, s, cfg) {
  d <- coefs[role == "primary" & stream == s & term != "(Intercept)" & !grepl("^splines::ns", term) & !is.na(or)]
  n_undrawn <- sum(!drawable_interval(d, "or", "or_low", "or_high", log = TRUE))
  d <- d[drawable_interval(d, "or", "or_low", "or_high", log = TRUE)]
  if (nrow(d) == 0L) return(NULL)
  d[, label := term_label(term)]
  d[, label := factor(label, levels = rev(unique(label)))]
  style <- cfg$plot$interval_style %||% "se"
  add_interval_columns(d, "or", "or_low", "or_high", log = TRUE, style = style)
  n_txt <- paste0(format(d$n_events[1L], big.mark = ","), " gaps, ", format(d$n_persons[1L], big.mark = ","), " people; gap >= 1 day, joint adjusted model.",
                  if (n_undrawn > 0L) paste0(" ", n_undrawn, " term(s) with too few events to estimate are not drawn (see the table)."))
  p <- ggplot(d, aes(x = est_, y = label)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
    interval_layer(style, horizontal = TRUE) +
    scale_x_log10(breaks = OR_BREAKS) + expand_limits(x = 1) +
    labs(title = paste0("All Predictors of a Gap Starting Tomorrow (", nice_label(s), " Stream)"),
         subtitle = paste0(n_txt, "\n", covariate_note(cfg)),
         x = "Odds ratio (log scale)", y = NULL, caption = interval_caption(style)) +
    forest_theme()
  attr(p, "size") <- c(10, forest_height(nrow(d), 1L) + 0.4)
  p
}

# Sensitivity analyses of ONE stream: every within-person predictor under every analysis (the pre-specified one is marked).
# `all_terms`: the sensitivity_all_terms table. The spline and random-slope analyses are tests, not estimates, and are not drawn.
SENSITIVITY_ORDER <- c("main", "full_period_standardization", "trailing_baseline", "exclude_partial_wear", "first_gap_only", "cci_removed", "cci_as_groups", "cci_within_cohort",
                       "main_sleep_only", "conditional_logit", "cluster_robust_logit")
plot_sensitivity_forest_stream <- function(all_terms, s, cfg) {
  d <- all_terms[stream == s & !is.na(or) & analysis %in% SENSITIVITY_ORDER & grepl("^t?w_", term)]   # the within-person predictors (all terms: sensitivity_all_terms_* figures)
  d <- d[drawable_interval(d, "or", "or_low", "or_high", log = TRUE)]
  if (nrow(d) == 0L) return(NULL)
  prim <- sub("^w_", "", cfg$primary_term)
  d[, pred := sub("^t?w_", "", term)]
  d[, panel := predictor_label(pred)]
  d[pred == prim, panel := paste0(panel, "  (pre-specified test)")]
  d[, panel := factor(panel, levels = unique(panel[order(pred != prim, match(pred, unique(pred)))]))]
  d[, label := factor(nice_label(analysis), levels = rev(nice_label(SENSITIVITY_ORDER)))]
  style <- cfg$plot$interval_style %||% "se"
  add_interval_columns(d, "or", "or_low", "or_high", log = TRUE, style = style)
  p <- ggplot(d, aes(x = est_, y = label)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
    interval_layer(style, horizontal = TRUE) +
    facet_wrap(~panel, ncol = 1) +
    scale_x_log10(breaks = OR_BREAKS) + expand_limits(x = 1) +
    labs(title = paste0("Sensitivity Analyses, All Within-Person Predictors (", nice_label(s), " Stream)"),
         x = "Odds ratio per +1 SD (log scale)", y = NULL, caption = interval_caption(style)) +
    forest_theme()
  attr(p, "size") <- c(9, forest_height(nrow(d), uniqueN(d$panel)))
  p
}

# Model B of ONE stream: every fixed effect of the four specifications (change in activity on the day back, in SD units).
GAPB_LABELS <- c(B_back2_unadj = "Unadjusted: first full day back", B_back2_adj = "Adjusted: first full day back",
                 B_back2_adj_pre = "Adjusted + steps before the gap: first full day back",
                 B_back1_adj = "Adjusted: first day back",
                 B_back2_adj_fullperiod = "Adjusted, full-period standardization: first full day back")
plot_gapB_forest <- function(all_terms, s, cfg) {
  d <- all_terms[stream == s & !is.na(estimate) & !is.na(se) & model_id %in% names(GAPB_LABELS)]
  if (nrow(d) == 0L) return(NULL)
  d[, `:=`(panel = factor(GAPB_LABELS[model_id], levels = GAPB_LABELS), label = term_label(term))]
  d[, label := factor(label, levels = rev(unique(term_label(c("log2_gap", "w_pre_day_steps", d$term)))))]
  style <- cfg$plot$interval_style %||% "se"
  add_interval_columns(d, "estimate", "conf_low", "conf_high", style = style)
  p <- ggplot(d, aes(x = est_, y = label)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
    interval_layer(style, horizontal = TRUE) +
    facet_wrap(~panel, ncol = 1, scales = "free_y") +
    labs(title = paste0("Model B, All Predictors (", nice_label(s), " Stream): Activity After a Gap"),
         subtitle = ov_wrap(paste0("Gap length: per doubling. Steps before the gap: per +1 SD above usual. ", covariate_note(cfg)), 150),
         x = "Change in activity on the day back (SD); 0 = no change",
         y = NULL, caption = interval_caption(style)) +
    forest_theme()
  attr(p, "size") <- c(10, forest_height(nrow(d), uniqueN(d$panel)))
  p
}

# Model C of ONE stream: every coefficient of the ordinal model as an odds ratio of a LONGER gap (only where proportional odds held).
GAPC_LABELS <- c(C_unadj = "Unadjusted", C_adj = "Adjusted", C_adj_fullperiod = "Adjusted, full-period standardization")
plot_gapC_forest <- function(all_terms, s, cfg) {
  d <- all_terms[stream == s & !is.na(or) & !is.na(se) & model_id %in% names(GAPC_LABELS)]
  d <- d[drawable_interval(d, "or", "or_low", "or_high", log = TRUE)]
  if (nrow(d) == 0L) return(NULL)
  d[, `:=`(panel = factor(GAPC_LABELS[model_id], levels = GAPC_LABELS), label = term_label(term))]
  d[, label := factor(label, levels = rev(unique(term_label(c("w_pre_day_steps", d$term)))))]
  style <- cfg$plot$interval_style %||% "se"
  add_interval_columns(d, "or", "or_low", "or_high", log = TRUE, style = style)
  p <- ggplot(d, aes(x = est_, y = label)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
    interval_layer(style, horizontal = TRUE) +
    facet_wrap(~panel, ncol = 1, scales = "free_y") +
    scale_x_log10(breaks = OR_BREAKS) + expand_limits(x = 1) +
    labs(title = paste0("Model C, All Predictors (", nice_label(s), " Stream): Odds of a Longer Gap"),
         x = "Odds ratio of a longer gap (log scale)",
         subtitle = ov_wrap(paste0("Steps before the gap: per +1 SD above usual. ", covariate_note(cfg)), 150),
         y = NULL, caption = interval_caption(style)) +
    forest_theme()
  attr(p, "size") <- c(10, forest_height(nrow(d), uniqueN(d$panel)))
  p
}

# The one extra figure per stream: the observed (model-free) % of days followed by a gap, by decile of yesterday's
# within-person steps, one line per cohort. It shows the shape of the association that the odds ratio summarises.
plot_onset_deciles <- function(dec, s, cfg) {
  d <- dec[stream == s & !is.na(rate)]
  if (nrow(d) == 0L) return(NULL)
  d[, `:=`(cohort = factor(cohort, levels = cfg$cohorts), dec = factor(decile))]
  style <- cfg$plot$interval_style %||% "se"
  add_interval_columns(d, "rate", "rate_low", "rate_high", scale = 100, style = style)
  dg <- position_dodge(width = 0.6)
  ggplot(d, aes(x = dec, y = est_, color = cohort, fill = cohort, group = cohort)) +
    geom_line(linewidth = 0.5, alpha = 0.6, position = dg) +
    interval_layer(style, horizontal = FALSE, color_var = "cohort", width = 0.5, position = dg) +
    scale_color_manual(values = cfg$palette$cohort) + scale_fill_manual(values = cfg$palette$cohort) +
    labs(title = paste0("Observed Gap Rate by Yesterday's Steps (", nice_label(s), " Stream)"),
         x = "Decile of yesterday's steps vs usual (1 = lowest, 10 = highest)",
         y = "% of days followed by a gap", color = "Cohort", fill = "Cohort", caption = interval_caption(style)) +
    theme_bw(base_size = 11) + theme(strip.text = element_text(face = "bold"), legend.position = "bottom")
}

# ---- what each figure shows (written next to the figures) ----------------------------------------------

# Captions by figure name; names ending in "_*" match every figure that starts with that prefix.
figure_captions <- function() c(
  device_version_missing_and_composition = "Watch version of each person, by cohort. One bar per cohort (100%): the grey segment at the left (Missing) is the share of people with NO usable device: no device record, scale only (a scale records weight only, no activity, sleep or heart rate) or rows without a version. The brown segment is people whose only device is a clip-on tracker or phone tracking: they have activity data but no watch. The other colours are the main watch family (the wrist-worn version with most records per person). Subtitle: Pearson chi-square of the missing share across cohorts; Monte-Carlo chi-square of the composition (missing included, and among people with a device only), each with the bias-corrected Cramer's V and its bootstrap interval. Wilson 95% intervals, pairwise comparisons with the reference cohort (Fisher exact, Holm), the logistic model adjusted for age and sex and the cell residuals are in the device_* tables. Descriptive; association only.",
  device_version_within_cohort = "The device versions within each cohort (people with a device). One bar per cohort and version (no intervals; Wilson intervals are in device_version_by_cohort.tsv): the share of the cohort's people with a device whose main device is that version; versions are grouped smartwatch / band / other / clip-on or phone. A * after a version means that its share differs between cohorts (Monte-Carlo chi-square of the cohort x has-this-version table, Holm over the versions). Subtitle: overall Monte-Carlo chi-square of the cohort x version table with the bias-corrected Cramer's V and its bootstrap interval. Descriptive; association only.",
  "outcome_by_variable_*" = "Person-level % of days missing (one value per person and stream) against ONE model variable, the streams side by side. Continuous variable: one line per cohort showing the mean % missing over 10 quantile bins of the variable (bins pooled over the cohorts), band = 95% interval of the mean, so the cohorts can be compared directly. Categorical variable: one violin per cohort within each level (the distribution of the person-level % missing, white dot = median). The statistics are written below the figure: the overall test and effect size (continuous: Spearman rho with its Bonett-Wright interval; categorical: Kruskal-Wallis or Wilcoxon with epsilon-squared or Cliff's delta and a bootstrap interval), Holm over the variables of the stream; for continuous variables the rank correlation within each cohort; and whether the association differs between cohorts (quasi-binomial variable x cohort or level x cohort F test, Holm). Between-person levels are each person's mean over observed days (geometric mean for log-transformed predictors). Cells with fewer than 10 people are not drawn. Descriptive; association only.",
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
  missing_pct_by_cci_group = "Person-level % of days missing by cohort and Charlson comorbidity group (violin + box). Label: Kruskal-Wallis between CCI groups within each cohort, Holm-adjusted over cohorts. Descriptive. The CCI score is a covariate of the models; its collinearity with cohort is reported (GVIF) and tested in sensitivity analyses.",
  "metrics_*" = "Person-level mean of each metric by cohort (or by cohort within a subgroup). Descriptive; KW p = Kruskal-Wallis between cohorts.",
  onset_primary_forest = "Primary gap-onset model (gap >= 1 day, joint, adjusted). Odds ratio that a gap starts tomorrow per +1 SD HIGHER than the person's own usual value yesterday (OR < 1: higher-than-usual days are followed by fewer gaps; OR > 1: by more gaps). Odds ratios, not risk ratios; conditional on the person's random intercept. Association, not causation.",
  forest_gap_definition = "The primary within-person steps effect when a gap is defined as >= 1, 2, 3 or 7 missing days (model D). Same odds-ratio reading as above.",
  risk_curve = "Population-average chance that a gap starts tomorrow when yesterday's steps were -1 SD / 0 (usual) / +1 SD from the person's usual (unadjusted joint model; bars = 95% interval from the intercept and slope only). Absolute risks for illustration; they depend on the random-intercept variance.",
  cohort_forest = "The primary effect fitted separately in each cohort. The label is the likelihood-ratio test (LRT) of the whole model with vs without a cohort x effect interaction: chisq with df = cohorts - 1, small p = the effect of prior-day steps differs between cohorts. It is one overall test of difference, not a comparison of single pairs, and intervals that overlap or do not overlap are not a substitute for it. Separate-cohort estimates are conditional on that cohort's own random intercept.",
  "wear_strata_onset_all_terms_other_*" = "Gap-onset model (primary, joint, adjusted) fitted separately in people with CONSISTENT wear (dark) and INCONSISTENT wear (red): every term except the intercept and the study-time spline (device family and CCI score included). Odds ratios; within and between terms per +1 SD, cohorts against Healthy. The strata contain different people and are compared with a likelihood-ratio test of (steps x stratum) in wear_strata_interaction.tsv - a difference is not read from which stratum's bar crosses 1. Strata are defined by observed-day coverage, which is built from the same days as the outcome (see wear_consistency in STATISTICS.md).",
  "wear_strata_gapB_all_terms_*" = "Model B (adjusted, first full day back) by wear consistency: every coefficient, in SD of the person's usual steps (per doubling of the gap, age per 10 years, cohorts against Healthy). Exploratory.",
  "wear_strata_gapC_all_terms_*" = "Model C (ordinal, adjusted) by wear consistency: odds ratios of a longer gap for every term (steps before the gap per +1 SD). Exploratory; only where proportional odds held in that stratum.",
  wear_strata_primary_effect = "The primary effect (odds ratio per +1 SD of prior-day steps) for people with consistent vs inconsistent wear, for the onset model and model C, one panel per stream. The label is the formal test that the strata differ (likelihood-ratio test of the interaction for the onset model; Wald test of the difference of two independent estimates for model C). Model B is compared in wear_strata_comparison.tsv.",
  wear_consistency_distribution = "Observed-day coverage (observed days / days between the person's first and last record) per stream, cumulative share of people by cohort. Dashed line: the classification cut-off (consistent >= cut-off); the number of people per stratum is written next to it. Descriptive.",
  wear_consistency_composition = "Standardized mean differences (consistent - inconsistent wear) in age, sex, CCI, follow-up length and cohort share, per stream. Dashed lines: |SMD| = 0.1. Shows how the strata differ in who they contain; it is not adjusted for anything.",
  "post_gap_activity_observed_*" = "Observed (model-free) steps after a gap, by day since the return (day 0 = first observed day after the gap, window of 8 days) and gap length. Value = steps (log) as a z-score against the person's usual BEFORE the gap, in within-person SD; 0 = usual, + = more steps. Dots: observed group means with the interval of the mean, clustered by person; n gaps and people are in each panel title (the number contributing to later days is smaller: a day without a record adds nothing). Red diamond: the Model B estimate (unadjusted) for day 1 at the group's median gap length; it is model-based and is not the same quantity as the observed group mean.",
  "collinearity_gvif_*" = "Multicollinearity of every term of the models of this stream: VIF-equivalent = GVIF^(1/df) from the covariance of the fixed effects (car::vif approach; terms with several coefficients - cohort, weekday, spline, year - are one term). Dashed lines: 5 (high) and 10 (severe). A high value widens the confidence interval of that term and makes it unstable; it does not bias the other estimates. No term was removed because of it. Tables: model_collinearity.tsv, model_term_tests.tsv.",
  "gap_definition_all_terms_other_*" = "Gap-onset model (joint, adjusted) with the gap defined as >= 1 (primary, dark), 2, 3 or 7 days: every term except the intercept and the study-time spline (device family and CCI score included). Odds ratios; within and between terms per +1 SD. Compares what the model estimates when only longer gaps count as the outcome.",
  "sensitivity_all_terms_other_*" = "ALL terms of the gap-onset model of this stream under the primary analysis (main, dark) and every sensitivity analysis that estimates them: full-period standardization, trailing 14-day baseline, no partial-wear days, first gap only, main-sleep-only (sleep), conditional logit (only the time-varying terms exist) and cluster-robust SE. Odds ratios; device family and CCI score included, calendar year is not a covariate. Spline and random-slope analyses are tests (sensitivity.tsv).",
  "exploratory_grid_effects_*" = "Exploratory specifications of this stream: the within- and between-person odds ratio of each predictor when it is the only predictor, when its group is together, and in the joint model, unadjusted and adjusted. Exploratory; no multiplicity adjustment. Covariate coefficients are in onset_coefficients.tsv.",
  "onset_forest_other_*" = "Gap-onset model of this stream (gap >= 1 day, joint, adjusted), every predictor, in this order: within-person terms (odds ratio that a gap starts tomorrow per +1 SD HIGHER than the person's own usual; OR < 1: higher-than-usual days are followed by fewer gaps), between-person terms (a person's usual level, per +1 SD), then covariates (cohort against Healthy, age per 10 years, sex against Female, weekday against Monday, device family against the largest family, CCI score per +1 point) and gap history (previous gap length and days since the last return per doubling of 1 + days, any earlier gap yes/no). The natural-spline terms of study time are in onset_coefficients.tsv. Odds ratios, not risk ratios; association, not causation.",
  "sensitivity_forest_*" = "Sensitivity analyses of this stream: every within-person predictor (one panel each; the pre-specified test is marked) under each analysis - the main fit, a trailing 14-day baseline, no partial-wear days, first gap only, main-sleep-only (sleep), conditional logit and cluster-robust SE. Odds ratio per +1 SD higher than the person's usual (OR < 1: higher-than-usual days are followed by fewer gaps). Estimates that agree across analyses are robust to those assumptions. The spline and random-slope analyses are likelihood-ratio tests and are in sensitivity.tsv.",
  "gapB_forest_*" = "Model B of this stream (exploratory), every fixed effect of the four specifications: change in activity on the first (full) day back, in SD of the person's usual, per doubling of the gap, per +1 SD of steps the day before the gap, per 10 years of age, and cohort (against Healthy) / sex (against Female) differences. 0 = no change. Conditional on a gap having occurred.",
  "gapB_mvpa_effect_*" = "Model B (exploratory) with MVPA before the gap as an additional covariate: the change in activity on the first full day back (SD of the person's usual steps) per +1 SD of MVPA minutes (log) on the last observed day BEFORE the gap (against the person's usual MVPA before that day), per doubling of the gap and per +1 SD of steps before the gap, in the specifications without MVPA (dark) and with MVPA (coloured). 0 = no change. MVPA is measured on the same day and standardized like the pre-gap steps; it is never the MVPA of the days after the gap (that is part of the outcome). Steps and MVPA of one day are correlated: see vif_equiv in gap_mvpa_effects.tsv.",
  "gapC_mvpa_effect_*" = "Model C (ordinal, exploratory; only where proportional odds held) with MVPA before the gap as an additional covariate: odds ratios of a LONGER gap per +1 SD of steps and of MVPA minutes (log) on the last observed day before the gap (against the person's usual before that day), without MVPA (dark) and with MVPA (red). Steps and MVPA of one day are correlated, so the steps odds ratio is the association beyond MVPA; see vif_equiv in gap_mvpa_effects.tsv.",
  "gapC_forest_*" = "Model C of this stream (exploratory; only where proportional odds held): every coefficient of the ordinal mixed model as an odds ratio of a LONGER gap (steps the day before the gap per +1 SD above usual, OR > 1 = longer gaps after more active days; age per 10 years; cohorts against Healthy, sex against Female).",
  "onset_deciles_*" = "Observed, model-free % of days followed by a gap, by decile of yesterday's within-person steps (1 = lowest relative to the person's usual), one line per cohort; the deciles are cut over all days at risk. Unadjusted; it shows the shape of the association that the odds ratio summarises (a falling line = higher-than-usual days are followed by fewer gaps). Intervals respect several days per person (person-clustered).",
  forest_sensitivity = "The primary effect under different assumptions (trailing 14-day baseline, no partial-wear days, first gap only, main-sleep-only, conditional logit, cluster-robust SE). Spline and random-slope analyses are likelihood-ratio tests and are listed in sensitivity.tsv, not drawn here.",
  gap_dose_response = "Model B (exploratory). Estimated mean activity on the first full day after a gap (SD from the person's usual) at gap lengths of 1, 2, 3, 5 and 13 days, one row per length, dark = unadjusted fit, red = adjusted fit (cohort, age, sex), with the interval of the estimate. The values come from the fitted straight line in log2(gap length), so they are model estimates at those lengths, not group means (observed group means are in model_B_observed.tsv). Conditional on a gap having occurred.",
  gap_length_probs = "Model C (exploratory, shown only where proportional odds held). Population-average probability of each gap-length group when the day before the gap was -1 SD / 0 (usual) / +1 SD from the person's usual steps, with the interval of the estimate. Not shown for streams where proportional odds was rejected (see model_C_gap_length.tsv)."
)

# figures that show an estimate with an interval (see interval_layer); `*_` entries match every figure with that prefix
INTERVAL_FIGURES <- c("transition", "onset_primary_forest", "forest_gap_definition", "gap_dose_response", "cohort_forest", "forest_sensitivity", "gap_length_probs",
                      "onset_forest_other_", "sensitivity_forest_", "gapB_forest_", "gapC_forest_", "gapB_mvpa_effect_", "gapC_mvpa_effect_", "onset_deciles_",
                      "gap_definition_all_terms_", "sensitivity_all_terms_", "wear_strata_", "post_gap_activity_observed_", "exploratory_grid_effects_")
interval_note <- function(style = "se") {
  if (identical(style, "box")) " One box per estimate: line = estimate, box = central 50% interval, whiskers = 95% interval - the uncertainty of the estimate, not the spread of the data (the 50% interval is derived from the 95% interval, exact for a normal sampling distribution). cfg$plot$interval_style = \"se\" (dot +/- 1 SD) or \"ci95\" (dot + 95% interval) draws the alternatives."
  else if (identical(style, "ci95")) " One dot per estimate with its 95% interval (the uncertainty of the estimate, not the spread of the data). cfg$plot$interval_style = \"se\" (dot +/- 1 SD) or \"box\" draws the alternatives."
  else " One dot per estimate with a bar of +/- 1 SD of the estimate (its standard error; about a 68% interval) - the uncertainty of the estimate, not the spread of the data. cfg$plot$interval_style = \"ci95\" (dot + 95% interval) or \"box\" (box + whiskers) draws the alternatives."
}

# Caption of a figure; `cfg` decides which interval style the note describes (transition has its own setting).
caption_for <- function(name, cfg = NULL) {
  cp <- figure_captions()
  style <- if (identical(name, "transition")) cfg$plot$transition_style %||% "box" else cfg$plot$interval_style %||% "se"
  is_interval <- name %in% INTERVAL_FIGURES || any(endsWith(INTERVAL_FIGURES, "_") & startsWith(name, INTERVAL_FIGURES))
  note <- if (is_interval) interval_note(style) else ""
  if (name %in% names(cp)) return(paste0(cp[[name]], note))
  pre <- names(cp)[endsWith(names(cp), "_*")]
  hit <- pre[startsWith(name, sub("\\*$", "", pre))]
  if (length(hit)) return(paste0(cp[[hit[1L]]], note))
  NA_character_
}
