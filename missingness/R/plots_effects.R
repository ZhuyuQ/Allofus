# ==============================================================================
# plots_effects.R
# ------------------------------------------------------------------------------
# Figures that compare estimates ACROSS analyses of the same model (primary vs sensitivity / gap definition / run-in
# variant), the multicollinearity figure, the observed post-gap activity figure and the run-in figures.
# They reuse the helpers and the look of plots.R (add_interval_columns, interval_layer, term_label, forest_theme,
# OR_BREAKS, cohort palette, theme_bw(11), bold strips, legend at the bottom).
# ==============================================================================

# Colours of the groups compared in one forest: the FIRST group is the primary analysis (dark navy), the others follow.
EFFECT_COLORS <- c("#2C3E50", "#E74C3C", "#3498DB", "#27AE60", "#F39C12", "#8E44AD", "#16A085", "#95A5A6")

# Label of a model TERM GROUP (a factor / spline / year term that needs several coefficients is one group).
term_group_label <- function(term, df = 1L) {
  fixed <- c(cohort = "Cohort", age10 = "Age", sex2 = "Sex", dow_t = "Weekday", `splines::ns(study_day, 3)` = "Study time (spline)",
             device = "Device family", cci = "CCI score", cci_missing = "CCI missing", cci_grp = "CCI group")
  lab <- ifelse(term %in% names(fixed), fixed[term], term_label(term))
  ifelse(df > 1L, paste0(lab, " [", df, " df]"), lab)
}

# One forest of ALL terms of one stream's model, comparing the groups of `group_col` (analyses, gap definitions, strata ...).
#   tab         one row per (group, term) with stream, term, group_col and the effect columns
#   group_labels named character vector: value of group_col -> legend label; its ORDER sets the order, the first is the primary
#   part        "other" (the only one): predictors first, then covariates and gap history
#   measure     "or": odds ratios (log axis, reference 1)   "estimate": estimates of a linear model (reference 0)
# The spline terms of study time and the intercept are never drawn: their coefficients have no reading one by one (they are
# tested jointly, see the term_tests tables).
plot_effects_by_group <- function(tab, s, cfg, group_col, group_labels, part = "other", title, subtitle = NULL,
                                  measure = c("or", "estimate"), xlab = NULL) {
  part <- match.arg(part); measure <- match.arg(measure)
  cols <- if (measure == "or") c("or", "or_low", "or_high") else c("estimate", "conf_low", "conf_high")
  d <- copy(tab)[stream == s]
  d[, term_key := sub("^tw_", "w_", term)]
  d <- d[is_reportable_term(term_key) & !is.na(get(cols[1L])) & as.character(get(group_col)) %in% names(group_labels)]
  d <- d[drawable_interval(d, cols[1L], cols[2L], cols[3L], log = measure == "or")]
  if (nrow(d) == 0L) return(NULL)
  group_labels <- group_labels[names(group_labels) %in% as.character(d[[group_col]])]          # only the groups this model has (no empty legend entries)
  d[, grp := factor(group_labels[as.character(get(group_col))], levels = unname(group_labels))]
  ord <- d[order(match(as.character(get(group_col)), names(group_labels)), match(term_key, unique(term_key)))]
  d[, label := term_label(term_key)]
  d[, label := factor(label, levels = rev(unique(term_label(ord$term_key))))]
  style <- cfg$plot$interval_style %||% "se"
  add_interval_columns(d, cols[1L], cols[2L], cols[3L], log = measure == "or", style = style)
  n_g <- uniqueN(d$grp); pal <- stats::setNames(EFFECT_COLORS[seq_len(n_g)], levels(d$grp)[seq_len(n_g)])
  ref <- if (measure == "or") 1 else 0
  p <- ggplot(d, aes(x = est_, y = label, color = grp, fill = grp)) +
    geom_vline(xintercept = ref, linetype = "dashed", color = "gray40") +
    interval_layer(style, horizontal = TRUE, color_var = "grp", width = 0.6, position = position_dodge(width = 0.75)) +
    scale_color_manual(values = pal, drop = FALSE) + scale_fill_manual(values = pal, drop = FALSE) +
    guides(color = guide_legend(reverse = TRUE, nrow = if (n_g > 3L) 2L else 1L), fill = guide_legend(reverse = TRUE, nrow = if (n_g > 3L) 2L else 1L)) +
    labs(title = title, subtitle = subtitle, color = NULL, fill = NULL, y = NULL, caption = interval_caption(style),
         x = xlab %||% if (measure == "or") "Odds ratio (log scale)" else "Estimate (SD units)") +
    forest_theme()
  if (measure == "or") p <- p + scale_x_log10(breaks = OR_BREAKS) + expand_limits(x = 1) else p <- p + expand_limits(x = 0)
  attr(p, "size") <- c(10, min(30, max(4, 1.3 + 0.3 * uniqueN(d$label) * max(1, 0.4 * n_g))))
  p
}

# The reading of every covariate, for the subtitles of the odds-ratio figures (the device reference is the one used in the models).
covariate_note <- function(cfg) {
  ref <- cfg$device$reference_used %||% cfg$device$reference
  paste0("Within and between: per +1 SD. Age: per 10 years. CCI: per +1 point. Previous gap length and days since return: per doubling of (1 + days). Cohorts against Healthy; sex against Female; weekdays against Monday",
         if (isTRUE(cfg$model$adjust_device)) paste0("; device family against ", ref %||% "the largest family") else "", ".")
}
subtitle_or <- function(cfg) paste0(ov_wrap(covariate_note(cfg), 100), "\nThe first colour (dark) is the primary analysis.")

# ---- multicollinearity ---------------------------------------------------------------------------------------------

# GVIF-based diagnostic (collinearity.R) of the models of ONE stream: one panel per model, one dot per term (a factor /
# spline / year set is one term with its df), x = VIF-equivalent GVIF^(1/df) on a log axis, dashed lines at 5 and 10.
#   tab : stream, analysis (panel label), term, df, vif_equiv, flag
plot_collinearity <- function(tab, s, cfg) {
  if (is.null(tab) || nrow(tab) == 0L) return(NULL)
  d <- copy(tab)[stream == s & !is.na(vif_equiv)]
  if (nrow(d) == 0L) return(NULL)
  d[, label := term_group_label(sub("^tw_", "w_", term), df)]
  d[, panel := factor(analysis, levels = unique(analysis))]
  d[, flag := factor(flag, levels = c("ok", "high", "severe"))]
  d[, label := factor(label, levels = rev(unique(label)))]
  ggplot(d, aes(x = pmax(vif_equiv, 1), y = label, color = flag)) +
    geom_vline(xintercept = VIF_HIGH, linetype = "dashed", color = "#E67E22") +
    geom_vline(xintercept = VIF_SEVERE, linetype = "dashed", color = "#E74C3C") +
    geom_point(size = 2.2) +
    facet_wrap(~panel, ncol = 1, scales = "free_y") +
    scale_x_log10(breaks = c(1, 2, 5, 10, 20, 50, 100)) +
    scale_color_manual(values = c(ok = "#2C3E50", high = "#E67E22", severe = "#E74C3C"), drop = FALSE) +
    labs(title = paste0("Multicollinearity of the Model Terms (", nice_label(s), " Stream)"),
         subtitle = ov_wrap("VIF-equivalent = GVIF^(1/df) from the covariance of the fixed effects. Dashed lines: 5 (high) and 10 (severe). A high value widens that term's interval; it does not bias the others.", 100),
         x = "VIF-equivalent (log scale)", y = NULL, color = "Flag") +
    forest_theme()
}

# ---- observed activity after gaps ------------------------------------------------------------------------------------

# Observed (model-free) steps after a gap by day since the return, one panel per gap-length group (observed_post_gap_trajectory).
#   Dots and bars: observed group mean of the z-score of steps against the person's usual BEFORE the gap, and +/- 1 SD of that
#   mean (person-clustered; the interval style follows cfg$plot$interval_style). Red diamonds: the Model B (unadjusted) estimate
#   for day 1, the day the model describes, at the group's median gap length - a model-based value, drawn apart from the data.
plot_post_gap_observed <- function(traj, curve, s, cfg) {
  d <- copy(traj)[stream == s]
  if (nrow(d) == 0L) return(NULL)
  style <- cfg$plot$interval_style %||% "se"
  add_interval_columns(d, "mean", "lo", "hi", style = style)
  n0 <- d[days_since_return == min(days_since_return), .(gap_group, n_gaps, n_persons)]
  lv <- levels(d$gap_group)
  lab <- stats::setNames(vapply(lv, function(g) {
    r <- n0[gap_group == g]
    if (nrow(r)) sprintf("%s%s: %s gaps, %s people", g, if (g == "all gaps") "" else if (g == "1") " day" else " days", format(r$n_gaps, big.mark = ","), format(r$n_persons, big.mark = ",")) else g
  }, ""), lv)
  d[, panel := factor(lab[as.character(gap_group)], levels = lab)]
  mb <- NULL
  if (!is.null(curve) && nrow(curve) > 0L) {
    cv <- curve[stream == s & model_id == "B_back2_unadj"]
    mdl <- unique(d[days_since_return == 1L, .(gap_group, med = pmax(1L, pmin(round(median_gap_length), max(cv$gap_length))))])
    mdl <- mdl[med <= max(cv$gap_length) & !gap_group %in% "15+"]
    mb <- merge(mdl, cv[, .(med = gap_length, est, lo, hi)], by = "med")
    if (nrow(mb)) { add_interval_columns(mb, "est", "lo", "hi", style = style); mb[, panel := factor(lab[as.character(gap_group)], levels = lab)] } else mb <- NULL
  }
  p <- ggplot(d, aes(x = days_since_return, y = est_)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
    geom_line(color = "#2C3E50", linewidth = 0.5, alpha = 0.6) +
    geom_pointrange(aes(ymin = lo_, ymax = hi_), color = "#2C3E50", size = 0.35)
  if (!is.null(mb)) p <- p + geom_pointrange(data = mb, aes(x = 1, y = est_, ymin = lo_, ymax = hi_), inherit.aes = FALSE,
                                             color = "#E74C3C", shape = 5, size = 0.5, stroke = 1)
  p + facet_wrap(~panel, ncol = 2) +
    scale_x_continuous(breaks = sort(unique(d$days_since_return))) +
    labs(title = paste0("Observed Steps After a Gap (", nice_label(s), " Stream)"),
         subtitle = "Day 0 = first observed day after the gap. Value = steps vs the person's usual before the gap (SD; 0 = usual, + = more).\nDots: observed means, interval of the mean (clustered by person).\nRed diamond: Model B estimate for day 1 (unadjusted, at the group's median gap length).",
         x = "Days since return", y = "Steps vs usual before the gap (SD)", caption = interval_caption(style)) +
    forest_theme() + theme(legend.position = "none")
}

# ---- run-in analysis (runin.R) ----------------------------------------------------------------------------------------------------
# Every figure states its hypothesis (H1 wear trait -> later gaps, H2 the steps effect depends on the run-in, H3 robustness of the primary
# effect) in the subtitle or caption; descriptive figures are unadjusted and use only days AFTER the run-in.

RUNIN_CAT_FILLS <- c("#DEEBF7", "#9ECAE1", "#4292C6", "#2171B5", "#08306B")      # sequential: no missing days (light) -> most (dark)
runin_cat_palette <- function(lv) stats::setNames(RUNIN_CAT_FILLS[round(seq(1, length(RUNIN_CAT_FILLS), length.out = length(lv)))], lv)
fmt_p <- function(p) ifelse(is.na(p), "NA", ifelse(p < 0.001, "<0.001", sprintf("%.2g", p)))
stream_levels <- function(x) factor(x, levels = intersect(c("activity", "sleep", "hr"), unique(x)))

subtitle_runin <- function(cfg) paste0(ov_wrap(paste0("Only days after the first ", cfg$runin$days, " days of follow-up are modelled (the run-in itself never is). Run-in term: per doubling of (1 + missing days in the run-in), centred. ",
                                                        covariate_note(cfg)), 100), "\nThe first colour (dark) is the model without a run-in term.")

# The tests of H2 as caption lines (statistics are written below the figure, outside the panels).
runin_h2_caption <- function(tests, which = "primary") {
  if (is.null(tests) || nrow(tests) == 0L) return(NULL)
  t2 <- tests[model %in% c("onset", "model C") & analysis == which & startsWith(hypothesis, "H2") & !is.na(p)]
  if (nrow(t2) == 0L) return(NULL)
  ov_wrap(paste0("H2, steps x run-in interaction", if (which == "sensitivity") " (sensitivity)" else "", ": ",
                 paste(t2[, sprintf("%s %s%s: %s p = %s, Holm p = %s", model, nice_label(stream), if (which == "sensitivity") paste0(" [", sub(",.*$", "", form), "]") else "", ifelse(!is.na(p_lrt), "LRT", "Wald"), fmt_p(p), fmt_p(p_holm))], collapse = " | ")), 150)
}

# Missing days in the run-in by cohort and stream (100% bars). `dist`: runin_distribution.
plot_runin_distribution <- function(dist, cfg) {
  if (is.null(dist) || nrow(dist) == 0L) return(NULL)
  d <- copy(dist)
  d[, `:=`(cohort = factor(cohort, levels = c(cfg$cohorts, "All cohorts")), category = factor(as.character(category), levels = runin_category_labels(cfg)), stream_f = stream_levels(stream))]
  d[, stream_f := factor(nice_label(as.character(stream_f)), levels = nice_label(levels(stream_f)))]
  d[, txt := fifelse(pct >= 4, sprintf("%.0f%%", pct), "")]
  d[, txtcol := fifelse(as.integer(category) >= 3L, "white", "black")]
  n_lab <- unique(d[, .(stream_f, cohort, total_people)])
  ggplot(d, aes(x = cohort, y = pct, fill = category)) +
    geom_col(position = position_stack(reverse = TRUE), width = 0.75) +
    geom_text(aes(label = txt, color = txtcol), position = position_stack(vjust = 0.5, reverse = TRUE), size = 3, show.legend = FALSE) +
    geom_text(data = n_lab, aes(x = cohort, y = 103, label = paste0("n = ", format(total_people, big.mark = ","))), inherit.aes = FALSE, size = 3, vjust = 0) +
    scale_color_identity() + scale_fill_manual(values = runin_cat_palette(levels(d$category)), drop = FALSE) +
    scale_y_continuous(limits = c(0, 109), breaks = seq(0, 100, 25)) +
    facet_wrap(~stream_f, nrow = 1) +
    labs(title = sprintf("Run-in: Missing Days in the First %d Days of Follow-up", cfg$runin$days),
         subtitle = ov_wrap("Description, no model. 0 = a record on every run-in day. Gaps analysed later start AFTER the run-in, so this measure always comes first. Question: how consistently do people wear the device at the start, and does it differ by cohort?", 100),
         x = NULL, y = "Share of people (%)", fill = "Missing days\nin the run-in") +
    forest_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
}

# H1, descriptive: the daily rate of new gaps after the run-in by run-in group. `rate`: runin_rate; `assoc`: runin_association (statistics in the caption).
plot_runin_rate <- function(rate, assoc, cfg) {
  if (is.null(rate) || nrow(rate) == 0L) return(NULL)
  d <- copy(rate)
  d[, `:=`(category = factor(category, levels = runin_category_labels(cfg)), cohort = factor(cohort, levels = c(cfg$cohorts, "All cohorts")), stream_f = stream_levels(stream))]
  d[, stream_f := factor(nice_label(as.character(stream_f)), levels = nice_label(levels(stream_f)))]
  d <- d[drawable_interval(d, "rate", "rate_low", "rate_high", log = FALSE)]
  if (nrow(d) == 0L) return(NULL)
  style <- cfg$plot$interval_style %||% "se"
  add_interval_columns(d, "rate", "rate_low", "rate_high", log = FALSE, style = style)
  pal <- c(cfg$palette$cohort, `All cohorts` = "#2C3E50")
  stat <- if (!is.null(assoc) && nrow(assoc)) vapply(unique(assoc$stream), function(st) {
    a <- assoc[stream == st]
    paste0(nice_label(st), ": ", paste(sprintf("%s rho = %.2f (%.2f to %.2f), Holm p %s", a$cohort, a$rho, a$rho_low, a$rho_high, fmt_p(a$p_holm)), collapse = "; "))
  }, "") else character()
  ggplot(d, aes(x = category, y = est_, color = cohort, fill = cohort, group = cohort)) +
    geom_line(data = function(x) x[stats::ave(seq_along(est_), stream_f, cohort, FUN = length) > 1L], position = position_dodge(width = 0.7), linewidth = 0.4, alpha = 0.5) +
    interval_layer(style, horizontal = FALSE, color_var = "cohort", width = 0.6, position = position_dodge(width = 0.7)) +
    scale_color_manual(values = pal) + scale_fill_manual(values = pal) +
    facet_wrap(~stream_f, nrow = 1, scales = "free_y") +
    expand_limits(y = 0) +
    labs(title = sprintf("H1: Missing Days in the First %d Days and Gaps Afterwards", cfg$runin$days),
         subtitle = ov_wrap("H1: people with more missing days in the run-in start new gaps more often later (rates rise from left to right). Unadjusted, descriptive; only days after the run-in. Interval: people are the clusters. The adjusted test of H1 is the run-in main effect (runin_tests.tsv).", 100),
         x = "Missing days in the run-in", y = "New gaps per 100 observed days (after the run-in)", color = NULL, fill = NULL,
         caption = paste(c(interval_caption(style), if (length(stat)) ov_wrap(paste0("Person-level Spearman correlation, run-in missing days vs the person's later rate: ", paste(stat, collapse = " | ")), 150)), collapse = "\n")) +
    forest_theme()
}

# H2, descriptive: gap rate by quintile of prior-day steps, one line per run-in group, one panel per cohort (one stream).
plot_runin_steps_by_group <- function(steps, s, cfg) {
  if (is.null(steps) || nrow(steps) == 0L) return(NULL)
  d <- copy(steps)[stream == s & n_events >= 20L]
  if (nrow(d) == 0L) return(NULL)
  d[, `:=`(category = factor(category, levels = runin_category_labels(cfg)), cohort = factor(cohort, levels = cfg$cohorts))]
  d <- d[drawable_interval(d, "rate", "rate_low", "rate_high", log = FALSE)]
  if (nrow(d) == 0L) return(NULL)
  style <- cfg$plot$interval_style %||% "se"
  add_interval_columns(d, "rate", "rate_low", "rate_high", log = FALSE, style = style)
  pal <- stats::setNames(c("#6BAED6", "#2171B5", "#E67E22", "#C0392B", "#7F2704")[seq_along(levels(d$category))], levels(d$category))
  ggplot(d, aes(x = mean_z, y = est_, color = category, fill = category, group = category)) +
    geom_line(data = function(x) x[stats::ave(seq_along(est_), cohort, category, FUN = length) > 1L], linewidth = 0.5, alpha = 0.7) +
    interval_layer(style, horizontal = FALSE, color_var = "category", width = 0.12) +
    scale_color_manual(values = pal, drop = FALSE) + scale_fill_manual(values = pal, drop = FALSE) +
    facet_wrap(~cohort, ncol = 2, scales = "free_y") + expand_limits(y = 0) +
    labs(title = paste0("H2: Gap Rate by Prior-Day Steps and Run-in (", nice_label(s), " Stream)"),
         subtitle = ov_wrap("H2: if the steps-gap association depends on the run-in, the lines differ in SLOPE between groups (not only in height). Unadjusted, descriptive; days after the run-in only. x: yesterday's steps vs the person's usual (mean z-score of each quintile; higher = more steps). Bins with fewer than 20 gaps are not drawn.", 100),
         x = "Yesterday's steps vs the person's usual (SD)", y = "New gaps per 100 observed days", color = "Missing days\nin the run-in", fill = "Missing days\nin the run-in",
         caption = interval_caption(style)) +
    forest_theme()
}

# H3 (robustness): the primary effect without and with each coding of the run-in. `pe`: runin_primary_effect; `tests`: runin_tests (H2 label); `variants` / `sensitivity`: the primary figure shows the reference and the run-in + interaction model, the sensitivity figure all forms.
plot_runin_primary_effect <- function(pe, tests, cfg, variants = RUNIN_PRIMARY_VARIANTS, sensitivity = FALSE) {
  if (is.null(pe) || nrow(pe) == 0L) return(NULL)
  d <- copy(pe)[grepl("odds", measure) & as.character(variant) %in% variants]
  d <- d[drawable_interval(d, "effect", "low", "high", log = TRUE)]
  if (nrow(d) == 0L) return(NULL)
  d[, `:=`(vlab = factor(RUNIN_LABELS[as.character(variant)], levels = rev(unname(RUNIN_LABELS[intersect(RUNIN_ONSET_VARIANTS, variants)]))), stream_f = stream_levels(stream))]
  d[, `:=`(stream_lab = factor(nice_label(as.character(stream_f)), levels = nice_label(levels(stream_f))), model = factor(model, levels = c("onset", "model C", "model B")))]
  d[, panel := factor(paste0(model, ": ", nice_label(as.character(stream_f))), levels = unique(paste0(model, ": ", nice_label(as.character(stream_f)))[order(model != "onset", stream_f)]))]
  cap_h2 <- runin_h2_caption(tests, if (sensitivity) "sensitivity" else "primary")
  style <- cfg$plot$interval_style %||% "se"
  add_interval_columns(d, "effect", "low", "high", log = TRUE, style = style)
  pal <- stats::setNames(EFFECT_COLORS[match(levels(d$vlab), unname(RUNIN_LABELS[RUNIN_ONSET_VARIANTS]))], levels(d$vlab))
  p <- ggplot(d, aes(x = est_, y = vlab, color = vlab, fill = vlab)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
    interval_layer(style, horizontal = TRUE, color_var = "vlab", show.legend = FALSE) +
    scale_color_manual(values = pal) + scale_fill_manual(values = pal) +
    scale_x_log10(breaks = OR_BREAKS) + expand_limits(x = 1) + facet_grid(model ~ stream_lab, scales = "free_y", space = "free_y")
  p + labs(title = if (sensitivity) "Sensitivity: Primary Effect With Other Forms of the Run-in" else "H3: Is the Primary Effect Explained by Wear in the First Days?", y = NULL, x = "Odds ratio per +1 SD of steps (log scale)",
           subtitle = ov_wrap(paste0("Same model, same days after the run-in; only the run-in terms differ. Onset: a gap starts tomorrow. Model C: a longer gap. 1 = no association. Similar estimates in every row = the steps effect is not explained by early wear. Run-in terms are centred, so the steps effect is at the average run-in.",
                                     if (sensitivity) " Here: the main-effect-only model, run-in as categories and run-in linear (onset)." else ""), 100),
           caption = paste(c(interval_caption(style), cap_h2), collapse = "\n")) + forest_theme()
}

# H2 (model): the steps effect at 0, 1, 3, 7 missing days in the run-in. `eff`: runin_effect_by_missing; `tests`: runin_tests.
plot_runin_steps_by_missing <- function(eff, tests, cfg) {
  if (is.null(eff) || nrow(eff) == 0L) return(NULL)
  d <- copy(eff)[model %in% c("onset", "model C")]
  d <- d[drawable_interval(d, "effect", "low", "high", log = TRUE)]
  if (nrow(d) == 0L) return(NULL)
  d[, `:=`(stream_f = stream_levels(stream), mf = factor(missing_days))]
  lv <- unique(d[order(model != "onset", stream_f), paste0(model, ": ", nice_label(as.character(stream_f)))])
  d[, panel := factor(paste0(model, ": ", nice_label(as.character(stream_f))), levels = lv)]
  style <- cfg$plot$interval_style %||% "se"
  add_interval_columns(d, "effect", "low", "high", log = TRUE, style = style)
  cap_h2 <- runin_h2_caption(tests)
  p <- ggplot(d, aes(x = mf, y = est_)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
    geom_line(data = function(x) x[stats::ave(seq_along(est_), panel, FUN = length) > 1L], aes(group = 1), color = "#2C3E50", alpha = 0.5) +
    interval_layer(style, horizontal = FALSE, color = "#2C3E50", width = 0.4) +
    scale_y_log10(breaks = OR_BREAKS) + facet_wrap(~panel, ncol = 3, scales = "free_y")
  p + labs(title = "H2: Does the Steps Effect Depend on Wear in the First Days?", x = "Missing days in the run-in", y = "Odds ratio per +1 SD of steps (log scale)",
           subtitle = ov_wrap("From the model with the steps x run-in interaction (run-in on the log2(1 + missing days) scale). Flat points = no modification; 1 = no association.", 100),
           caption = paste(c(interval_caption(style), cap_h2), collapse = "\n")) + forest_theme()
}

# Shape check of the run-in coding (onset model): categories vs the log2 coding at each group's median.
plot_runin_shape <- function(coef_on, models, dist, shape, cfg) {
  if (is.null(coef_on) || nrow(coef_on) == 0L || is.null(dist) || nrow(dist) == 0L) return(NULL)
  cat <- coef_on[runin_variant == "categorical" & grepl("^cov0_cat", term) & !is.na(or), .(stream, category = sub("^cov0_cat", "", term), effect = or, low = or_low, high = or_high)]
  main <- coef_on[runin_variant == "main" & term == "cov0", .(stream, b = estimate, se)]
  med <- dist[cohort == "All cohorts", .(stream, category = as.character(category), m = median_missing_days_in_category)][!is.na(m)]
  fit <- merge(med, main, by = "stream")
  if (nrow(cat) == 0L || nrow(fit) == 0L) return(NULL)
  fit[, L := log2(1 + m)]
  fit <- fit[, .(stream, category, effect = exp(b * L), low = exp((b - 1.96 * se) * L), high = exp((b + 1.96 * se) * L))]
  ref <- fit[category == runin_category_labels(cfg)[1L], .(stream, category, effect = 1, low = 1, high = 1)]
  d <- rbindlist(list(cbind(cat, coding = "Categories (no assumed shape)"), cbind(fit[category != runin_category_labels(cfg)[1L]], coding = "log2(1 + missing days) at the group median"),
                      cbind(ref, coding = "Categories (no assumed shape)")), use.names = TRUE)
  d <- d[drawable_interval(d, "effect", "low", "high", log = TRUE)]
  d[, `:=`(category = factor(category, levels = runin_category_labels(cfg)), stream_f = stream_levels(stream))]
  d[, stream_f := factor(nice_label(as.character(stream_f)), levels = nice_label(levels(stream_f)))]
  style <- cfg$plot$interval_style %||% "se"
  add_interval_columns(d, "effect", "low", "high", log = TRUE, style = style)
  aic <- if (!is.null(shape) && nrow(shape)) shape[variant == "categorical" & !is.na(delta_aic_vs_main), sprintf("%s: AIC(categories) - AIC(log2) = %.1f", nice_label(stream), delta_aic_vs_main)]
  pal <- c("Categories (no assumed shape)" = "#2C3E50", "log2(1 + missing days) at the group median" = "#E74C3C")
  ggplot(d, aes(x = category, y = est_, color = coding, fill = coding)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
    interval_layer(style, horizontal = FALSE, color_var = "coding", width = 0.5, position = position_dodge(width = 0.6)) +
    scale_color_manual(values = pal) + scale_fill_manual(values = pal) +
    scale_y_log10(breaks = OR_BREAKS) + facet_wrap(~stream_f, nrow = 1) +
    labs(title = "H1 (shape): Odds of a Later Gap by Missing Days in the Run-in", x = "Missing days in the run-in", y = "Odds ratio of a gap starting (log scale), against 0 missing days",
         color = NULL, fill = NULL,
         subtitle = ov_wrap("Onset model after the run-in. Points: each group against 0 missing days, no shape assumed. Second colour: what the log2(1 + missing days) coding of the main model implies at the group's median. Close agreement supports the log scale; a gap between them means the shape is not log-linear.", 100),
         caption = paste(c(interval_caption(style), if (length(aic)) paste(aic, collapse = " | ")), collapse = "\n")) +
    forest_theme()
}

# ---- the exploratory grid ---------------------------------------------------------------------------------------------------

# Label of an onset specification set.
spec_set_label <- function(set_name) {
  vapply(set_name, function(x) {
    if (x == "joint") "All predictors together"
    else if (x == "own_multi") "Own-stream predictors together"
    else if (x == "act_multi") "Activity predictors together"
    else if (grepl("^(own_uni|act_uni|uni)_", x)) paste(predictor_label(sub("^(own_uni|act_uni|uni)_", "", x)), "alone")
    else x
  }, character(1), USE.NAMES = FALSE)
}

# Exploratory specifications of ONE stream: the within- and between-person odds ratio of every predictor in every specification
# (alone / grouped / all together), unadjusted and adjusted. The adjustment covariates are in onset_coefficients.tsv.
plot_exploratory_grid <- function(coefs, s, cfg) {
  d <- coefs[stream == s & k == cfg$onset$primary_k & role %in% c("primary", "exploratory") & grepl("^[wb]_", term) & !is.na(or)]
  d <- d[drawable_interval(d, "or", "or_low", "or_high", log = TRUE)]              # intervals that cannot be estimated cannot be drawn on a log axis
  if (nrow(d) == 0L || !any(d$role == "exploratory")) return(NULL)
  d[, `:=`(spec = spec_set_label(set_name), adj = factor(ifelse(adjusted, "Adjusted", "Unadjusted"), levels = c("Unadjusted", "Adjusted")),
           panel = term_label(term))]
  d[, panel := factor(panel, levels = unique(panel))]
  d[, spec := factor(spec, levels = rev(unique(spec)))]
  style <- cfg$plot$interval_style %||% "se"
  add_interval_columns(d, "or", "or_low", "or_high", log = TRUE, style = style)
  p <- ggplot(d, aes(x = est_, y = spec, color = adj, fill = adj)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
    interval_layer(style, horizontal = TRUE, color_var = "adj", width = 0.6, position = position_dodge(width = 0.7)) +
    facet_wrap(~panel, ncol = 2) +
    scale_color_manual(values = c(Unadjusted = "#95A5A6", Adjusted = "#2C3E50")) + scale_fill_manual(values = c(Unadjusted = "#95A5A6", Adjusted = "#2C3E50")) +
    guides(color = guide_legend(reverse = TRUE), fill = guide_legend(reverse = TRUE)) +
    scale_x_log10(breaks = OR_BREAKS) + expand_limits(x = 1) +
    labs(title = paste0("Exploratory Specifications: Odds Ratios of Every Predictor (", nice_label(s), " Stream)"),
         subtitle = "Each row is a model that contains the predictors named in the row; within and between: per +1 SD. Exploratory, not multiplicity-adjusted.",
         x = "Odds ratio (log scale)", y = NULL, color = NULL, fill = NULL, caption = interval_caption(style)) +
    forest_theme()
  attr(p, "size") <- c(10, max(4, 1.5 + 0.35 * uniqueN(d$spec) * ceiling(uniqueN(d$panel) / 2)))
  p
}

# ---- MVPA before the gap in models B and C ---------------------------------------------------------------------------------

GAPB_MVPA_LABELS <- c(B_back2_adj = "Adjusted (cohort, age, sex)", B_back2_adj_mvpa = "Adjusted + MVPA before the gap",
                      B_back2_adj_pre = "Adjusted + steps before the gap", B_back2_adj_pre_mvpa = "Adjusted + steps and MVPA before the gap")
GAPC_MVPA_LABELS <- c(C_adj = "Adjusted (cohort, age, sex)", C_adj_mvpa = "Adjusted + MVPA before the gap")

# The MVPA covariate of ONE stream, with and without it: the gap-length term, the steps term and the MVPA term of each
# specification (the other covariates are in gapB_forest_ / gapC_forest_). Dark = without MVPA, coloured = with MVPA.
#   all_terms : model_B_all_terms or model_C_all_terms
plot_gap_mvpa <- function(all_terms, s, cfg, model = c("B", "C")) {
  model <- match.arg(model)
  lab <- if (model == "B") GAPB_MVPA_LABELS else GAPC_MVPA_LABELS
  if (is.null(all_terms) || nrow(all_terms) == 0L) return(NULL)
  d <- all_terms[stream == s & model_id %in% names(lab) & term %in% c("log2_gap", "w_pre_day_steps", "w_pre_day_mvpa") & !is.na(estimate)]
  if (!any(d$term == "w_pre_day_mvpa")) return(NULL)                               # no MVPA fit for this stream
  # keep a specification without MVPA only next to its MVPA version
  has <- unique(d$model_id[d$term == "w_pre_day_mvpa"])
  keep <- c(has, GAP_MVPA_PAIRS[has]); d <- d[model_id %in% keep]
  p <- if (model == "B") {
    plot_effects_by_group(d, s, cfg, "model_id", lab, "other", paste0("Model B: MVPA Before the Gap (", nice_label(s), " Stream)"),
                          "Change in activity on the first full day back, per +1 SD of MVPA / steps\nbefore the gap and per doubling of the gap. Dark: without MVPA; coloured: with MVPA.",
                          measure = "estimate", xlab = "Change in activity on the day back (SD); 0 = no change")
  } else {
    plot_effects_by_group(d, s, cfg, "model_id", lab, "other", paste0("Model C: MVPA Before the Gap (", nice_label(s), " Stream)"),
                          "Odds of a longer gap per +1 SD of steps / MVPA before the gap (ordinal model).\nDark: without MVPA; red: with MVPA.",
                          measure = "or", xlab = "Odds ratio of a longer gap (log scale)")
  }
  sz <- attr(p, "size")
  p <- p + guides(color = guide_legend(reverse = TRUE, nrow = 2), fill = guide_legend(reverse = TRUE, nrow = 2))
  attr(p, "size") <- sz
  p
}

# ---- device version: missing share and composition by cohort, in ONE figure -----------------------------------------------

# One horizontal 100% bar per cohort: the grey segment at the left is the share of people WITHOUT a watch version
# (no device record, scale only, no version), the brown segment is people with only a clip-on or phone (activity data, no watch), the other colours are the main watch family. The row
# labels carry n and the missing share; the subtitle carries the tests (device_tests.tsv); the Wilson intervals are in the tables.
plot_device_composition <- function(ds, cfg) {
  comp <- ds$device_composition_by_cohort; mb <- ds$device_missing_by_cohort; tt <- ds$device_tests
  if (is.null(comp) || nrow(comp) == 0L) return(NULL)
  miss_lab <- cfg$device$missing_label
  mb <- mb[cohort != "All cohorts"]
  lv <- cfg$cohorts[cfg$cohorts %in% comp$cohort]
  lab <- mb[, .(cohort, lab = sprintf("%s\nn = %s\nmissing %.1f%%", cohort, format(n_people, big.mark = ",", trim = TRUE), pct_missing))]
  d <- copy(comp)
  d[, cohort := factor(cohort, levels = rev(lv), labels = rev(lab$lab[match(lv, lab$cohort)]))]
  cats <- levels(comp$category)
  d[, category := factor(as.character(category), levels = rev(cats))]       # stacking puts the first level at the far end: Missing is drawn first (left)
  pal <- cfg$palette$device[cats]; pal[is.na(pal)] <- "#7F7F7F"
  dark <- vapply(pal, function(h) { rgb <- grDevices::col2rgb(h); sum(rgb * c(0.299, 0.587, 0.114)) < 140 }, logical(1))
  d[, txt_col := ifelse(dark[as.character(category)], "white", "black")]
  d[, lab_txt := ifelse(pct >= 4, sprintf("%.0f%%", pct), "")]
  f <- function(r, col, fmt = "%.2f") if (nrow(r) && is.finite(r[[col]])) sprintf(fmt, r[[col]]) else "NA"
  pf <- function(x) if (is.na(x)) "NA" else if (x < 0.001) "<0.001" else sprintf("%.3f", x)
  t1 <- tt[grepl("^Missing", test)]; t2 <- tt[grepl("missing included", test)]; t3 <- tt[grepl("known device", test)]
  sub <- paste0("Missing differs between cohorts: chi-square(", t1$df, ") = ", f(t1, "chisq", "%.1f"), ", p ", ifelse(t1$p_asymptotic < 0.001, "<0.001", paste0("= ", pf(t1$p_asymptotic))),
                ", Cramer's V = ", f(t1, "cramers_v"), " (95% CI ", f(t1, "v_low"), "-", f(t1, "v_high"), ").\n",
                "Composition differs, missing included: Monte-Carlo p ", ifelse(t2$p_monte_carlo < 0.001, "<0.001", paste0("= ", pf(t2$p_monte_carlo))),
                ", V = ", f(t2, "cramers_v"), " (", f(t2, "v_low"), "-", f(t2, "v_high"), ").\n",
                "Composition differs, known devices only: Monte-Carlo p ", ifelse(t3$p_monte_carlo < 0.001, "<0.001", paste0("= ", pf(t3$p_monte_carlo))),
                ", V = ", f(t3, "cramers_v"), " (", f(t3, "v_low"), "-", f(t3, "v_high"), ").")
  ggplot(d, aes(x = pct, y = cohort, fill = category, group = category)) +
    geom_col(width = 0.72, color = "white", linewidth = 0.3) +
    geom_text(aes(label = lab_txt, color = txt_col, group = category), position = position_stack(vjust = 0.5), size = 3, show.legend = FALSE) +
    scale_color_identity() +
    scale_fill_manual(values = pal, breaks = cats, drop = FALSE) +
    scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20), labels = function(x) paste0(x, "%"), expand = expansion(mult = c(0, 0.01))) +
    labs(title = "Watch Version: Missing and Composition by Cohort", subtitle = sub, x = "Share of people in the cohort", y = NULL, fill = "Main device",
         caption = ov_wrap(paste0("Grey (Missing) = no device record, scale only (a scale records weight only) or no version. Brown = clip-on or phone only: activity data, no watch. Other colours = the watch version with most records per person, grouped into families. ",
                                  "Labels below 4% are not printed; all values and their Wilson 95% intervals are in device_composition_by_cohort.tsv."), 150)) +
    guides(fill = guide_legend(nrow = 2)) +
    forest_theme()
}

# ---- person-level % of days missing against every model variable --------------------------------------------------

ov_wrap <- function(x, width) paste(strwrap(x, width), collapse = "\n")
ov_fmt_p <- function(p) ifelse(is.na(p), "NA", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))

# ONE figure per model variable, the streams that have it side by side, ONE row:
#   continuous variable  one line per cohort: mean % of days missing over 10 quantile bins of the variable (bins pooled over the cohorts), band = 95% interval of the mean;
#   categorical variable one violin per cohort within each level (the distribution of the person-level % missing; white dot = median).
# All statistics are written BELOW the figure (outside the panels): the overall test and effect size, the correlation within each cohort and the test of whether the
# association differs between cohorts, per stream.
#   ovs: list (by stream) of outcome_variable_analysis() results;  var: the variable name (e.g. "cohort", "age_first", "act_steps")
plot_outcome_variable <- function(ovs, var, cfg) {
  oc <- cfg$outcome_vars
  streams <- names(ovs)[vapply(ovs, function(a) !is.null(a$tests) && var %in% a$tests$variable, logical(1))]
  if (!length(streams)) return(NULL)
  is_cat <- var %in% names(OV_CATEGORICAL)
  cells <- rbindlist(lapply(streams, function(s) ovs[[s]]$cells[variable == var]))
  cells <- cells[n >= oc$min_cell]
  if (nrow(cells) == 0L) return(NULL)
  tt <- rbindlist(lapply(streams, function(s) ovs[[s]]$tests[variable == var]), fill = TRUE)
  bc <- rbindlist(lapply(streams, function(s) ovs[[s]]$by_cohort[variable == var]), fill = TRUE)
  var_label <- tt$variable_label[1L]
  sl <- stats::setNames(nice_label(streams), streams); lv_stream <- unname(sl[streams])
  ord_cat <- unique(c(cfg$cohorts, c("Female", "Male", "Other"), cfg$strata$cci_labels, "Unknown", levels(ovs[[streams[1L]]]$data$watch)))
  cells[, `:=`(xid = paste(stream, key, sep = "|"), cohort = factor(cohort, levels = cfg$cohorts))]
  cells[, ord := if (is_cat) match(x_label, ord_cat) else as.integer(key)]
  setorder(cells, stream, ord)
  xlev <- unique(cells$xid); xlab <- stats::setNames(cells$x_label[match(xlev, cells$xid)], xlev)
  cells[, `:=`(xid = factor(xid, levels = xlev), stream_f = factor(sl[stream], levels = lv_stream))]

  # ---- the statistics, written below the figure ----
  stat_lines <- vapply(streams, function(s) {
    t1 <- tt[stream == s]; r <- bc[stream == s]
    first <- if (is_cat) paste0(ifelse(t1$test == "Kruskal-Wallis", sprintf("Kruskal-Wallis H(%d) = %.1f", as.integer(t1$df), t1$statistic), "Wilcoxon rank-sum"), ", p(Holm) ", ov_fmt_p(t1$p_holm), "; ",
                                ifelse(t1$effect_name == "epsilon-squared", "epsilon-squared = ", "Cliff's delta = "), sprintf("%.3f (95%% CI %.3f to %.3f)", t1$effect, t1$effect_low, t1$effect_high))
             else paste0("Spearman rho = ", sprintf("%.3f (95%% CI %.3f to %.3f)", t1$effect, t1$effect_low, t1$effect_high), ", p(Holm) ", ov_fmt_p(t1$p_holm))
    byc <- if (!is_cat && nrow(r)) paste0("; rho within cohort: ", paste0(r$cohort, " ", sprintf("%.2f", r$rho), collapse = ", ")) else ""
    dif <- if ("p_cohort_difference" %in% names(t1) && !is.na(t1$p_cohort_difference)) paste0("; cohorts differ (", ifelse(is_cat, "pattern across levels", "slope"), "): p(Holm) ", ov_fmt_p(t1$p_cohort_difference_holm)) else ""
    paste(strwrap(paste0(sl[[s]], ": ", first, byc, dif), 150, exdent = 4), collapse = "\n")
  }, character(1))
  note <- ov_wrap(paste0(if (is_cat) "Violin = distribution of the person-level % of days missing, white dot = median; cohorts side by side within each level. "
                         else "Lines = mean % missing per cohort within 10 bins of the variable (bins from everybody pooled), band = 95% interval of the mean. ",
                         "Tests: overall rank test (Holm over the variables of the stream); the cohort-difference test is the ", if (is_cat) "level x cohort" else "variable x cohort",
                         " F test of the quasi-binomial model (Holm). Cells with fewer than ", oc$min_cell, " people are not drawn. Details: outcome_variable_tests / _by_cohort / _contrasts / _models / _cells (.tsv)."), 150)
  cap <- paste0("Statistics\n", paste(stat_lines, collapse = "\n"), "\n\n", note)
  n_lines <- length(strsplit(cap, "\n", fixed = TRUE)[[1L]])

  dodge <- position_dodge(width = 0.85)
  p <- ggplot() +
    facet_grid(. ~ stream_f, scales = "free_x") +
    scale_x_discrete(labels = xlab) +
    scale_color_manual(values = cfg$palette$cohort, breaks = cfg$cohorts) + scale_fill_manual(values = cfg$palette$cohort, breaks = cfg$cohorts) +
    labs(title = paste0("% of Days Missing by ", var_label), x = var_label, y = "% of days missing (one value per person)", color = "Cohort", fill = "Cohort", caption = cap) +
    theme_bw(base_size = 11) + theme(strip.text = element_text(face = "bold"), legend.position = "bottom", panel.grid.minor = element_blank(),
                                      axis.text.x = element_text(angle = if (is_cat) 35 else 45, hjust = 1, size = 8),
                                      plot.caption = element_text(size = 8.5, hjust = 0, lineheight = 1.05, margin = margin(t = 10)))
  if (is_cat) {
    pr <- rbindlist(lapply(streams, function(s) {
      d <- ovs[[s]]$data[!is.na(get(var)) & !is.na(pct_missing), .(stream = s, key = as.character(get(var)), cohort = factor(as.character(cohort), levels = cfg$cohorts), y = pct_missing)]
      d[, xid := paste(stream, key, sep = "|")]
      d[paste(xid, cohort) %in% paste(cells$xid, cells$cohort)]
    }))
    pr[, `:=`(xid = factor(xid, levels = xlev), stream_f = factor(sl[stream], levels = lv_stream))]
    p <- p + geom_violin(data = pr, aes(x = xid, y = y, fill = cohort), position = dodge, scale = "width", width = 0.8, color = NA, alpha = 0.85) +
      stat_summary(data = pr, aes(x = xid, y = y, group = cohort), fun = stats::median, geom = "point", position = dodge, shape = 21, fill = "white", color = "gray25", size = 1.3)
  } else {
    p <- p + geom_ribbon(data = cells, aes(x = xid, ymin = low, ymax = high, fill = cohort, group = cohort), alpha = 0.18) +
      geom_line(data = cells, aes(x = xid, y = mean_pct, color = cohort, group = cohort), linewidth = 0.8) +
      geom_point(data = cells, aes(x = xid, y = mean_pct, color = cohort), size = 1.5)
  }
  attr(p, "size") <- c(max(7, 3.4 * length(streams) + 2.4), 5.4 + 0.2 * n_lines)
  p
}

# The device VERSIONS within each cohort (people with a device): one row per version, grouped smartwatch / band / other / clip-on or phone, one BAR per
# cohort = % of the cohort's people with a device whose main device is that version (no intervals; they are in device_version_by_cohort.tsv). A * after the
# version means the share differs between the cohorts (Monte-Carlo chi-square of the cohort x "has this version" table, Holm over the versions); the subtitle
# has the overall test.
plot_device_versions <- function(ds, cfg) {
  ver <- ds$device_version_by_cohort; vt <- ds$device_version_tests; tt <- ds$device_tests
  if (is.null(ver) || nrow(ver) == 0L) return(NULL)
  grp_lv <- c("Smartwatch", "Band", "Other", cfg$device$nonwrist_label)
  vt <- copy(vt)[, `:=`(rare = grepl("^Other .* versions$", device_version))]
  vt[, `:=`(group = factor(group, levels = grp_lv))]
  vt <- vt[order(group, rare, -n_people)]
  vt[, lab := sprintf("%s  (n = %s)%s", device_version, format(n_people, big.mark = ",", trim = TRUE), ifelse(!is.na(p_holm) & p_holm < 0.05, "  *", ""))]
  d <- merge(copy(ver)[, device_version := as.character(device_version)], vt[, .(device_version, lab, group)], by = "device_version", suffixes = c("_old", ""))
  d[, `:=`(lab = factor(lab, levels = rev(vt$lab)), group = factor(group, levels = grp_lv), cohort = factor(as.character(cohort), levels = rev(cfg$cohorts)))]
  d[, txt := ifelse(pct_of_known >= 1, sprintf("%.1f%%", pct_of_known), "")]
  nk <- unique(d[, .(cohort, n_known)])
  clab <- stats::setNames(sprintf("%s (n = %s)", nk$cohort, format(nk$n_known, big.mark = ",", trim = TRUE)), as.character(nk$cohort))
  t4 <- tt[grepl("VERSION", test)]
  sub <- paste0("Share of each cohort's people with a device whose main device is this version (n = people with a device; each cohort's bars add up to 100%).\n",
                "Versions differ between cohorts overall: Monte-Carlo p ", ifelse(t4$p_monte_carlo < 0.001, "<0.001", sprintf("= %.3f", t4$p_monte_carlo)),
                ", Cramer's V = ", sprintf("%.2f", t4$cramers_v), " (95% CI ", sprintf("%.2f", t4$v_low), "-", sprintf("%.2f", t4$v_high), "). * = this version's share differs between cohorts, p(Holm) < 0.05.")
  dodge <- position_dodge(width = 0.82)
  p <- ggplot(d, aes(x = pct_of_known, y = lab, fill = cohort)) +
    geom_col(position = dodge, width = 0.78) +
    geom_text(aes(label = txt), position = dodge, hjust = -0.15, size = 2.3, color = "gray25") +
    facet_grid(group ~ ., scales = "free_y", space = "free_y") +
    scale_fill_manual(values = cfg$palette$cohort, breaks = cfg$cohorts, labels = clab[cfg$cohorts]) +
    scale_x_continuous(labels = function(x) paste0(x, "%"), expand = expansion(mult = c(0, 0.12))) +
    labs(title = "Device Versions Within Each Cohort", subtitle = ov_wrap(sub, 115), x = "% of the cohort's people with a device", y = NULL, fill = "Cohort",
         caption = ov_wrap(paste0("Each person counts once, under the version with most records (wrist-worn devices first). Versions used by fewer than ", cfg$device$min_version_people,
                                  " people are pooled as 'Other ... versions'. Labels below 1% are not printed. Counts, Wilson intervals and the per-version tests: device_version_by_cohort.tsv, device_version_tests.tsv."), 125)) +
    guides(fill = guide_legend(nrow = 1)) +
    forest_theme() + theme(plot.subtitle = element_text(size = 8.5), plot.caption = element_text(size = 8, hjust = 0), strip.text.y = element_text(angle = 0, face = "bold"),
                           panel.grid.major.y = element_blank())
  attr(p, "size") <- c(11, max(5, 2.2 + 0.62 * nrow(vt)))
  p
}
