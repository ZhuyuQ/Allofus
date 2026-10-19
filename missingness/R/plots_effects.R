# ==============================================================================
# plots_effects.R
# ------------------------------------------------------------------------------
# Figures that compare estimates ACROSS analyses of the same model (primary vs sensitivity / gap definition / wear
# stratum), the multicollinearity figure, the observed post-gap activity figure and the wear-consistency figures.
# They reuse the helpers and the look of plots.R (add_interval_columns, interval_layer, term_label, forest_theme,
# OR_BREAKS, cohort palette, theme_bw(11), bold strips, legend at the bottom).
# ==============================================================================

# Colours of the groups compared in one forest: the FIRST group is the primary analysis (dark navy), the others follow.
EFFECT_COLORS <- c("#2C3E50", "#E74C3C", "#3498DB", "#27AE60", "#F39C12", "#8E44AD", "#16A085", "#95A5A6")

# Label of a model TERM GROUP (a factor / spline / year term that needs several coefficients is one group).
term_group_label <- function(term, df = 1L) {
  fixed <- c(cohort = "Cohort", age10 = "Age", sex2 = "Sex", dow_t = "Weekday", `splines::ns(study_day, 3)` = "Study time (spline)",
             `factor(year_t)` = "Calendar year")
  lab <- ifelse(term %in% names(fixed), fixed[term], term_label(term))
  ifelse(df > 1L, paste0(lab, " [", df, " df]"), lab)
}

# One forest of ALL terms of one stream's model, comparing the groups of `group_col` (analyses, gap definitions, strata ...).
#   tab         one row per (group, term) with stream, term, group_col and the effect columns
#   group_labels named character vector: value of group_col -> legend label; its ORDER sets the order, the first is the primary
#   part        "other": every term except the calendar years (predictors first, then covariates and gap history)
#               "year" : the calendar-year terms only
#   measure     "or": odds ratios (log axis, reference 1)   "estimate": estimates of a linear model (reference 0)
# The spline terms of study time and the intercept are never drawn: their coefficients have no reading one by one (they are
# tested jointly, see the term_tests tables).
plot_effects_by_group <- function(tab, s, cfg, group_col, group_labels, part = c("other", "year"), title, subtitle = NULL,
                                  measure = c("or", "estimate"), xlab = NULL) {
  part <- match.arg(part); measure <- match.arg(measure)
  cols <- if (measure == "or") c("or", "or_low", "or_high") else c("estimate", "conf_low", "conf_high")
  d <- copy(tab)[stream == s]
  d[, term_key := sub("^tw_", "w_", term)]
  d <- d[is_reportable_term(term_key) & !is.na(get(cols[1L])) & as.character(get(group_col)) %in% names(group_labels)]
  is_year <- grepl("^factor\\(year_t\\)", d$term_key)
  d <- if (part == "year") d[is_year] else d[!is_year]
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
    guides(color = guide_legend(reverse = TRUE), fill = guide_legend(reverse = TRUE)) +
    labs(title = title, subtitle = subtitle, color = NULL, fill = NULL, y = NULL, caption = interval_caption(style),
         x = xlab %||% if (measure == "or") "Odds ratio (log scale)" else "Estimate (SD units)") +
    forest_theme()
  if (measure == "or") p <- p + scale_x_log10(breaks = OR_BREAKS) + expand_limits(x = 1) else p <- p + expand_limits(x = 0)
  attr(p, "size") <- c(10, min(30, max(4, 1.3 + 0.3 * uniqueN(d$label) * max(1, 0.4 * n_g))))
  p
}

SUBTITLE_OR <- paste0("Within and between: per +1 SD. Age: per 10 years. Cohorts are compared with Healthy; sex with Female; weekdays with Sunday;",
                      "\nother terms per +1 unit. The first colour (dark) is the primary analysis.")

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
         subtitle = "VIF-equivalent = GVIF^(1/df) from the covariance of the fixed effects.\nDashed lines: 5 (high) and 10 (severe). A high value widens that term's interval; it does not bias the others.",
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

# ---- wear consistency ---------------------------------------------------------------------------------------------------

# Distribution of observed-day coverage per stream (ECDF by cohort) with the classification cut-off.
plot_wear_distribution <- function(by_person, criteria, cfg) {
  d <- copy(by_person); d[, cohort := factor(cohort, levels = cfg$cohorts)]
  thr <- copy(criteria)
  ggplot(d, aes(x = coverage, color = cohort)) +
    stat_ecdf(linewidth = 0.8) +
    geom_vline(data = thr, aes(xintercept = threshold), linetype = "dashed", color = "gray30", inherit.aes = FALSE) +
    geom_text(data = thr, aes(x = threshold, y = 0.12, label = sprintf("cut-off %.3f\n%s consistent / %s inconsistent", threshold,
                                                                       format(n_consistent, big.mark = ","), format(n_inconsistent, big.mark = ","))),
              hjust = 1.03, size = 3, inherit.aes = FALSE) +
    facet_wrap(~stream, ncol = 1) +
    scale_color_manual(values = cfg$palette$cohort) +
    labs(title = "Wear Consistency: Share of Days With a Record", x = "Observed-day coverage (observed days / days between first and last record)",
         y = "Cumulative share of people", color = "Cohort") +
    forest_theme()
}

# Differences between the two wear-consistency strata in the people's characteristics (standardized mean differences).
plot_wear_composition <- function(comp, cfg) {
  d <- copy(comp)[!is.na(smd)]
  d[, variable := factor(variable, levels = rev(unique(variable)))]
  ggplot(d, aes(x = smd, y = variable, color = stream)) +
    geom_vline(xintercept = 0, color = "gray50") +
    geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", color = "gray60") +
    geom_point(size = 2.4, position = position_dodge(width = 0.5)) +
    scale_color_manual(values = c(activity = "#2C3E50", sleep = "#3498DB", hr = "#E74C3C")) +
    labs(title = "Who Wears the Device Consistently?", x = "Standardized mean difference (Consistent - Inconsistent)", y = NULL, color = "Stream",
         subtitle = "Dashed lines: |SMD| = 0.1. Positive = higher / more common among consistent wearers.") +
    forest_theme()
}

# The primary effect in the two strata, side by side, for the odds-ratio models (onset and model C) and every stream, with
# the formal test of the difference. `cmp`: wear_strata_comparison(). Model B (a linear model, SD units) is drawn in its own
# forest (wear_strata_gapB_all_terms_<stream>) and compared in the table.
plot_wear_primary_effect <- function(cmp, cfg) {
  if (is.null(cmp) || nrow(cmp) == 0L) return(NULL)
  cmp <- cmp[grepl("odds", effect_measure)]
  if (nrow(cmp) == 0L) return(NULL)
  long <- rbindlist(lapply(WEAR_LEVELS, function(lv) {
    cmp[, .(model, stream, stratum = lv, est = get(paste0("effect_", tolower(lv))), lo = get(paste0("low_", tolower(lv))),
            hi = get(paste0("high_", tolower(lv))))]
  }))
  long[, stratum := factor(stratum, levels = WEAR_LEVELS)]
  style <- cfg$plot$interval_style %||% "se"
  add_interval_columns(long, "est", "lo", "hi", log = TRUE, style = style)
  long[, panel := factor(paste0(model, ": ", stream))]
  lab <- cmp[, .(panel = factor(paste0(model, ": ", stream), levels = levels(long$panel)),
                 txt = sprintf("%s: p = %s", ifelse("p_interaction_lrt" %in% names(cmp) & model == "onset (primary)", "LRT, strata differ", "Wald, strata differ"),
                               ifelse(p_difference < 0.001, "<0.001", format(signif(p_difference, 2)))))]
  ggplot(long, aes(x = est_, y = stratum, color = stratum, fill = stratum)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
    interval_layer(style, horizontal = TRUE, color_var = "stratum", show.legend = FALSE) +
    geom_text(data = lab, aes(x = Inf, y = Inf, label = txt), hjust = 1.05, vjust = 1.5, size = 3, fontface = "italic", inherit.aes = FALSE) +
    facet_wrap(~panel, ncol = 3) +
    scale_color_manual(values = c(Consistent = "#2C3E50", Inconsistent = "#E74C3C")) + scale_fill_manual(values = c(Consistent = "#2C3E50", Inconsistent = "#E74C3C")) +
    scale_x_log10(breaks = OR_BREAKS) + expand_limits(x = 1) +
    labs(title = "Primary Effect by Wear Consistency", y = NULL, x = "Odds ratio per +1 SD of steps (log scale)",
         subtitle = "Onset model: odds that a gap starts tomorrow. Model C: odds of a longer gap. 1 = no association.",
         caption = interval_caption(style)) +
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
# (no device record, scale only, phone / clip-on only, ...), the coloured segments are the main watch family of the others. The row
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
         caption = paste0("Grey = no watch version (no device record, scale only, phone or clip-on only, or no version). Colours = the watch version with most records per person, grouped into families.\n",
                          "Labels below 4% are not printed; all values and their Wilson 95% intervals are in device_composition_by_cohort.tsv.")) +
    guides(fill = guide_legend(nrow = 2)) +
    forest_theme()
}

# The same cohort x watch-version table read the other way: one 100% bar per watch version (and per "Missing"), coloured by COHORT. The top
# bar is the cohort mix of all people: if cohort and watch version were unrelated, every bar would look like it.
plot_device_cohort_mix <- function(ds, cfg) {
  mx <- ds$device_cohort_mix; tt <- ds$device_tests
  if (is.null(mx) || nrow(mx) == 0L) return(NULL)
  d <- copy(mx)
  cats <- levels(d$category)
  lab <- unique(d[, .(category, n_category)])[, lab := sprintf("%s\nn = %s", category, format(n_category, big.mark = ",", trim = TRUE))]
  d[, category := factor(category, levels = rev(cats), labels = rev(lab$lab[match(cats, lab$category)]))]
  d[, cohort := factor(cohort, levels = rev(cfg$cohorts))]
  pal <- cfg$palette$cohort
  dark <- vapply(pal, function(h) { rgb <- grDevices::col2rgb(h); sum(rgb * c(0.299, 0.587, 0.114)) < 140 }, logical(1))
  d[, `:=`(lab_txt = ifelse(pct_of_category >= 4, sprintf("%.0f%%", pct_of_category), ""), txt_col = ifelse(dark[as.character(cohort)], "white", "black"))]
  t2 <- tt[grepl("missing included", test)]
  pf <- function(x) if (is.na(x)) "NA" else if (x < 0.001) "<0.001" else sprintf("%.3f", x)
  sub <- paste0("Cohort and watch version are associated: Monte-Carlo p ", ifelse(t2$p_monte_carlo < 0.001, "<0.001", paste0("= ", pf(t2$p_monte_carlo))),
                ", Cramer's V = ", sprintf("%.2f", t2$cramers_v), " (95% CI ", sprintf("%.2f", t2$v_low), "-", sprintf("%.2f", t2$v_high), ").\n",
                ov_wrap("The top bar is the cohort mix of everybody; a bar that differs from it shows where a cohort is over- or under-represented (cell residuals: device_composition_residuals.tsv).", 130))
  ggplot(d, aes(x = pct_of_category, y = category, fill = cohort, group = cohort)) +
    geom_col(width = 0.72, color = "white", linewidth = 0.3) +
    geom_text(aes(label = lab_txt, color = txt_col, group = cohort), position = position_stack(vjust = 0.5), size = 3, show.legend = FALSE) +
    scale_color_identity() +
    scale_fill_manual(values = pal, breaks = cfg$cohorts) +
    scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20), labels = function(x) paste0(x, "%"), expand = expansion(mult = c(0, 0.01))) +
    labs(title = "Cohort Mix of Each Watch Version", subtitle = sub, x = "Share of the people with this watch version", y = NULL, fill = "Cohort",
         caption = "Each bar = 100% of the people with that watch version (Missing = no watch version). Labels below 4% are not printed; Wilson 95% intervals are in device_cohort_mix.tsv.") +
    forest_theme()
}

# ---- person-level % of days missing against every model variable --------------------------------------------------

ov_wrap <- function(x, width) paste(strwrap(x, width), collapse = "\n")
ov_fmt_p <- function(p) ifelse(is.na(p), "NA", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))

# Categorical variables (cohort, sex, watch version) in ONE figure: a column per variable, a row per stream; violin + box of the person-level
# % of days missing by level, n per level under the violin, the overall test and its effect size in the corner of each panel.
#   ovs: list (by stream) of outcome_variable_analysis() results
plot_outcome_categorical <- function(ovs, cfg) {
  min_n <- cfg$outcome_vars$min_n
  rows <- rbindlist(lapply(names(ovs), function(s) {
    d <- ovs[[s]]$data
    rbindlist(lapply(intersect(names(OV_CATEGORICAL), names(d)), function(v)
      d[!is.na(get(v)) & !is.na(pct_missing), .(stream = s, variable = v, level = as.character(get(v)), value = pct_missing)]))
  }))
  if (nrow(rows) == 0L) return(NULL)
  rows[, n := .N, by = .(stream, variable, level)]
  rows <- rows[n >= min_n]                                                       # levels with fewer people are not drawn (nor tested)
  lv <- unique(c(cfg$cohorts, c("Female", "Male", "Other"), if ("watch" %in% names(ovs[[1L]]$data)) levels(ovs[[1L]]$data$watch)))
  rows[, `:=`(level = factor(level, levels = lv), variable_f = factor(OV_CATEGORICAL[variable], levels = OV_CATEGORICAL), stream_f = factor(nice_label(stream), levels = nice_label(names(ovs))))]
  pal <- c(cfg$palette$cohort, cfg$palette$sex, cfg$palette$device); pal <- pal[!duplicated(names(pal))]
  cnt <- unique(rows[, .(stream_f, variable_f, level, n)])
  tt <- rbindlist(lapply(ovs, `[[`, "tests"))[type == "categorical"]
  tt[, `:=`(variable_f = factor(OV_CATEGORICAL[variable], levels = OV_CATEGORICAL), stream_f = factor(nice_label(stream), levels = nice_label(names(ovs))))]
  tt[, label := paste0(ifelse(test == "Kruskal-Wallis", sprintf("Kruskal-Wallis H(%d) = %.1f", as.integer(df), statistic), "Wilcoxon rank-sum"),
                       "\np(Holm) ", ov_fmt_p(p_holm), "\n", ifelse(effect_name == "epsilon-squared", "epsilon-sq = ", "Cliff's delta = "), sprintf("%.3f", effect),
                       "\n(95% CI ", sprintf("%.3f", effect_low), " to ", sprintf("%.3f", effect_high), ")")]
  dodge <- position_dodge(width = 0.9)
  ggplot(rows, aes(x = level, y = value, fill = level)) +
    geom_violin(color = NA, alpha = 0.8, scale = "width", width = 0.85) +
    geom_boxplot(width = 0.12, outlier.size = 0.3, outlier.alpha = 0.5, color = "gray25", fill = "white", alpha = 0.7) +
    geom_text(data = cnt, aes(x = level, y = -Inf, label = paste0("n = ", format(n, big.mark = ","))), inherit.aes = FALSE, vjust = -0.5, size = 2.6, color = "gray25") +
    geom_text(data = tt, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE, hjust = -0.03, vjust = 1.15, size = 2.7, fontface = "italic", lineheight = 0.95) +
    facet_grid(stream_f ~ variable_f, scales = "free", space = "free_x") +
    scale_fill_manual(values = pal, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0.07, 0.38))) +
    labs(title = "Person-Level % of Days Missing by Model Variable (Categorical)", x = NULL, y = "% of days missing (one value per person)",
         subtitle = ov_wrap("Violin = distribution, white box = median and quartiles. Corner: overall rank test across the levels (Holm-adjusted over all variables of the stream) and its effect size with a bootstrap 95% interval.", 150),
         caption = ov_wrap(paste0("Levels with fewer than ", min_n, " people are not drawn or tested. Watch version: Missing = no watch version (see device_version_missing_and_composition). Post-hoc contrasts and adjusted odds ratios: outcome_variable_contrasts.tsv, outcome_variable_models.tsv."), 160)) +
    theme_bw(base_size = 11) + theme(strip.text = element_text(face = "bold"), axis.text.x = element_text(angle = 35, hjust = 1), plot.subtitle = element_text(size = 8.5),
                                      plot.caption = element_text(size = 8, hjust = 0), panel.grid.minor = element_blank())
}

# Continuous variables of ONE stream: age and the between-person level of every model predictor. Light points = people (coloured by
# cohort, a seeded sample), black = mean % missing in quantile bins with the 95% interval of the mean; the corner gives the Spearman
# correlation (Bonett-Wright interval, Holm over the variables of the stream) and the adjusted odds ratio per +1 SD.
plot_outcome_continuous <- function(a, stream, cont_spec, cfg) {
  d <- a$data; bins <- a$bins
  if (is.null(bins) || nrow(bins) == 0L) return(NULL)
  cv <- c(list(age_first = list(label = "Age at first record (years)", display = "age_first", log_x = FALSE)), cont_spec)
  cv <- cv[names(cv) %in% unique(bins$variable)]
  lab_levels <- vapply(cv, `[[`, "", "label")
  set.seed(cfg$seed)
  ids <- sample_person_ids(d$person_id, cfg$outcome_vars$points_max, cfg$seed)
  ymax <- stats::quantile(d$pct_missing, 0.99, na.rm = TRUE)
  pts <- rbindlist(lapply(names(cv), function(v) {
    x <- d[person_id %in% ids & !is.na(get(cv[[v]]$display)), .(x = get(cv[[v]]$display), y = pct_missing, cohort)]
    q <- stats::quantile(x$x, c(0.005, 0.995), na.rm = TRUE)
    x[x >= q[1L] & x <= q[2L] & y <= ymax][, variable_f := cv[[v]]$label]
  }))
  pts[, `:=`(variable_f = factor(variable_f, levels = lab_levels), cohort = factor(cohort, levels = cfg$cohorts))]
  bb <- copy(bins)[, variable_f := factor(variable_label, levels = lab_levels)]
  tt <- a$tests[type == "continuous"]; md <- a$models[adjustment != "unadjusted" & term == "per +1 SD"]
  tt <- merge(tt, md[, .(variable, or = odds_ratio, or_low, or_high)], by = "variable", all.x = TRUE)
  tt[, variable_f := factor(variable_label, levels = lab_levels)]
  tt[, label := paste0("Spearman rho = ", sprintf("%.3f", effect), " (", sprintf("%.3f", effect_low), " to ", sprintf("%.3f", effect_high), "), p(Holm) ", ov_fmt_p(p_holm),
                       ifelse(is.na(or), "", paste0("\nadjusted OR per +1 SD = ", sprintf("%.2f", or), " (", sprintf("%.2f", or_low), "-", sprintf("%.2f", or_high), ")")))]
  ncol <- min(3L, length(cv)); nr <- ceiling(length(cv) / ncol)
  p <- ggplot() +
    geom_point(data = pts, aes(x = x, y = y, color = cohort), alpha = 0.18, size = 0.7) +
    geom_line(data = bb, aes(x = x_display, y = mean_pct), color = "black", linewidth = 0.5) +
    geom_pointrange(data = bb, aes(x = x_display, y = mean_pct, ymin = low, ymax = high), color = "black", size = 0.35) +
    geom_text(data = tt, aes(x = -Inf, y = Inf, label = label), hjust = -0.03, vjust = 1.2, size = 2.7, fontface = "italic", lineheight = 0.95) +
    facet_wrap(~variable_f, scales = "free_x", ncol = ncol) +
    scale_color_manual(values = cfg$palette$cohort, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.18))) +
    coord_cartesian(ylim = c(0, ymax * 1.18)) +
    labs(title = paste0("Person-Level % of Days Missing by Model Variable (Continuous): ", nice_label(stream), " Stream"),
         subtitle = ov_wrap("Black: mean % missing within quantile bins of the variable with its 95% interval. Points: people (seeded sample), coloured by cohort. Corner: rank correlation across all people and the odds ratio of a day being missing per +1 SD of the variable, adjusted for cohort, sex and age (quasi-binomial model); Holm adjustment over the variables of the stream.", 150),
         x = NULL, y = "% of days missing (one value per person)", color = "Cohort",
         caption = ov_wrap("Between-person level = the person's mean over observed days (geometric mean for log-transformed predictors). Points outside the 0.5-99.5% range of x or above the 99th percentile of y are not drawn; bins and tests use everybody.", 160)) +
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 2))) +
    theme_bw(base_size = 11) + theme(strip.text = element_text(face = "bold"), legend.position = "bottom", plot.subtitle = element_text(size = 8.5),
                                      plot.caption = element_text(size = 8, hjust = 0), panel.grid.minor = element_blank())
  attr(p, "size") <- c(13, 3.6 * nr + 2.2)
  p
}
