# ==============================================================================
# plots_stream_model.R - the three views of the stream model (stream_model.R) and the writing of its tables.
#   1 model-free: the share of risk-set rows followed by each category, by decile of w (one panel per category, own y-axis, person-clustered 95% intervals)
#   2 adjusted association: forest of the odds ratios per +1 SD (common log axis, the 0.97-1.03 equivalence region shaded, the negative control marked, the probably-worn
#     split of category 3 as secondary rows) and the stacked-model contrasts 1 vs 2 and 1 vs 3 with Holm-adjusted intervals
#   3 absolute risk: the population-average probability of each category across w from -2 to +2 (one panel per category, own y-axis; models A and B side by side) and the
#     stacked mix of all five statuses (model B)
#   plus the events per category and cohort against the thresholds, and the sensitivity forest (primary estimate first)
# ==============================================================================

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

STM_MODEL_LABELS <- c(c1 = "1 All missing (main outcome)", c2 = "2 Activity only", c3 = "3 Sleep missing", c4 = "4 Inconsistent (negative control)",
                      c3_probably_worn_20 = "3a Sleep missing, probably worn (20 h)", c3_probably_not_worn_20 = "3b Sleep missing, probably not worn (20 h)",
                      c3_probably_worn_22 = "3a Sleep missing, probably worn (22 h)", c3_probably_not_worn_22 = "3b Sleep missing, probably not worn (22 h)")
STM_ROLE_COLORS <- c(`main outcome` = "#C0392B", comparison = "#2C3E50", `negative control` = "#8E44AD", secondary = "#7F8C8D")

stm_role <- function(model) fcase(model == "c1", "main outcome", model %in% c("c2", "c3"), "comparison", model == "c4", "negative control", default = "secondary")

# View 1: observed shares by decile of w.
plot_stm_deciles <- function(dec) {
  d <- copy(dec)[, label := factor(label, levels = unique(label))]
  ggplot(d, aes(decile, share_pct)) + geom_pointrange(aes(ymin = low, ymax = high), size = 0.4, color = "#2C3E50") + geom_line(color = "grey60") +
    facet_wrap(~label, nrow = 1, scales = "free_y") + scale_x_continuous(breaks = 1:10) +
    labs(title = "Observed Share of Days Followed by Each Status, by Decile of Yesterday's Steps", x = "Decile of w (within-person deviation of steps on d - 1; 10 = most above usual)",
         y = "% of risk-set rows followed by the category", caption = "Model-free. 95% intervals with people as clusters. Each panel has its own y-axis (baseline rates differ by orders of magnitude).") +
    theme_bw(base_size = 10) + theme(strip.text = element_text(face = "bold", size = 8))
}

# View 2: forest of adjusted odds ratios per +1 SD and the stacked-model contrasts.
plot_stm_forest <- function(tests, contrasts, scfg) {
  d <- copy(tests)[variant == "main"][, `:=`(label = factor(STM_MODEL_LABELS[model], levels = rev(STM_MODEL_LABELS)), role = stm_role(model))]
  A <- ggplot(d, aes(or, label, color = role)) + annotate("rect", xmin = scfg$equiv_margin[1L], xmax = scfg$equiv_margin[2L], ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.6) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
    geom_linerange(aes(xmin = or_low, xmax = or_high), linewidth = 0.5) + geom_linerange(aes(xmin = or90_low, xmax = or90_high), linewidth = 1.6) + geom_point(aes(shape = role), size = 2.5) +
    scale_x_log10() + scale_color_manual(values = STM_ROLE_COLORS) + scale_shape_manual(values = c(`main outcome` = 16, comparison = 16, `negative control` = 17, secondary = 1)) +
    labs(title = "Adjusted Association of Yesterday's Steps With Tomorrow's Status", x = "Odds ratio per +1 SD of w (log scale)", y = NULL, color = NULL, shape = NULL,
         caption = "Thin line: 95% interval; thick line: 90% interval (equivalence). Grey band: 0.97-1.03. Subject-specific odds ratios (mixed models).") + theme_bw(base_size = 10) + theme(legend.position = "bottom")
  c2 <- copy(contrasts)[variant == "main"]
  B <- ggplot(c2, aes(ratio_of_or, contrast)) + geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") + geom_linerange(aes(xmin = ratio_low_holm, xmax = ratio_high_holm), linewidth = 0.8, color = "#2C3E50") +
    geom_point(size = 2.5, color = "#2C3E50") + scale_x_log10() + labs(title = "Contrasts (Multinomial Model)", x = "Ratio of odds ratios (Holm-adjusted 95% interval)", y = NULL,
         caption = "Population-averaged scale; one multinomial model with person-clustered errors.") + theme_bw(base_size = 10)
  out <- gridExtra::arrangeGrob(A, B, ncol = 2, widths = c(2, 1.2)); attr(out, "size") <- c(13, 5.5); out
}

# View 3: absolute risk across w, model A (mixed, random intercept integrated out) and model B (multinomial).
plot_stm_absolute <- function(curve) {
  d <- curve[variant == "main" & model %in% c("c1", "c2", "c3", "c4")][, label := factor(STM_MODEL_LABELS[model], levels = STM_MODEL_LABELS[c("c1", "c2", "c3", "c4")])]
  d[, source := factor(source, c("multinomial", "mixed"), c("B multinomial", "A mixed"))]
  ggplot(d, aes(w, prob_pct, color = source, fill = source)) + geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.15, color = NA) + geom_line(linewidth = 0.8) +
    facet_wrap(~label, nrow = 1, scales = "free_y") + scale_x_continuous(breaks = -2:2) + scale_color_manual(values = c("#2C3E50", "#C0392B")) + scale_fill_manual(values = c("#2C3E50", "#C0392B")) +
    labs(title = "Population-Average Probability of Each Status on the Next Day", x = "w (SD of steps on day d above the person's usual)", y = "% of days (population average)", color = NULL, fill = NULL,
         caption = "Averaged over the observed covariates (model A also over the random intercept); 95% bands from 2,000 coefficient draws. Each panel has its own y-axis.") +
    theme_bw(base_size = 10) + theme(strip.text = element_text(face = "bold", size = 8), legend.position = "bottom")
}

# The mix of the five statuses across w (model B probabilities add to 100%).
plot_stm_stacked <- function(curve) {
  d <- curve[variant == "main" & source == "multinomial"][, label := factor(STM_LABELS[sub("c", "", model)], levels = rev(STM_LABELS))]
  ggplot(d, aes(w, prob_pct, fill = label)) + geom_area(position = "stack", alpha = 0.9) + scale_x_continuous(breaks = -2:2, expand = c(0, 0)) + scale_y_sqrt(expand = c(0, 0)) +
    scale_fill_manual(values = setNames(c("#BDC3C7", "#C0392B", "#2C3E50", "#E67E22", "#8E44AD"), STM_LABELS)) +
    labs(title = "How the Mix of Next-Day Statuses Changes With Steps", x = "w (SD of steps on day d above the person's usual)", y = "% of days (square-root scale)", fill = NULL,
         caption = "Multinomial model, covariates as observed; the segments add to 100%. Square-root axis so the rare statuses are visible.") + theme_bw(base_size = 10)
}

# Events per category and cohort against the thresholds; the events lost by the wear rule (relative to the 16 h set) are drawn on top when that variant is present.
plot_stm_events <- function(ec, scfg) {
  m <- ec[category > 0L & group != "all"]; mn <- m[variant == "main"][, label := factor(label, levels = STM_LABELS[-1L])]
  lost <- NULL
  if ("wear_ge16h" %in% ec$variant) { w16 <- ec[variant == "wear_ge16h" & category > 0L & group != "all", .(group, category, e16 = events)]
    lost <- merge(mn[, .(group, category, label, events)], w16, by = c("group", "category"))[, lost := pmax(0, e16 - events)] }
  p <- ggplot(mn, aes(group, events)) + geom_col(fill = "#2C3E50")
  if (!is.null(lost)) p <- p + geom_col(data = lost, aes(y = lost), fill = "grey70", position = position_stack())
  p + geom_hline(yintercept = scfg$min_events_effect, linetype = "dashed", color = "#C0392B") + geom_hline(yintercept = scfg$min_events_equiv, linetype = "dotted", color = "#8E44AD") +
    facet_wrap(~label, nrow = 1, scales = "free_y") + labs(title = "Events per Next-Day Status and Cohort", x = NULL, y = "events (rows followed by the status)",
      caption = "Red dashed: 1,500 events (effect claims); purple dotted: 6,000 (equivalence claims). Grey: events lost by the 20 h wear rule (against the 16 h set).") +
    theme_bw(base_size = 10) + theme(axis.text.x = element_text(angle = 30, hjust = 1))
}

# Sensitivity forest: the odds ratio of every category under every variant (the main estimate first).
plot_stm_sensitivity <- function(tests, scfg) {
  d <- tests[model %in% c("c1", "c2", "c3", "c4")][, `:=`(label = factor(STM_MODEL_LABELS[model], levels = STM_MODEL_LABELS[c("c1", "c2", "c3", "c4")]),
                                                       variant = factor(variant, levels = rev(unique(c("main", variant)))))]
  ggplot(d, aes(or, variant)) + annotate("rect", xmin = scfg$equiv_margin[1L], xmax = scfg$equiv_margin[2L], ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.6) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") + geom_linerange(aes(xmin = or_low, xmax = or_high), linewidth = 0.6, color = "#2C3E50") + geom_point(size = 2, color = "#2C3E50") +
    facet_wrap(~label, nrow = 1, scales = "free_x") + scale_x_log10() +
    labs(title = "Sensitivity of the Odds Ratios per +1 SD of w", x = "Odds ratio (log scale, 95% interval)", y = NULL, caption = "Grey band: 0.97-1.03. Each panel has its own x-axis.") + theme_bw(base_size = 10)
}

stm_figures <- function(res, dec, scfg, event_counts = NULL) {
  f <- list(stream_model_observed_deciles = plot_stm_deciles(dec), stream_model_forest = plot_stm_forest(res$tests, res$contrasts, scfg), stream_model_absolute_risk = plot_stm_absolute(res$risk_curve))
  if (nrow(res$risk_curve[source == "multinomial"])) f$stream_model_mix_by_steps <- plot_stm_stacked(res$risk_curve)
  if (!is.null(event_counts)) f$stream_model_events <- plot_stm_events(event_counts, scfg)
  if (uniqueN(res$tests$variant) > 1L) f$stream_model_sensitivity <- plot_stm_sensitivity(res$tests, scfg)
  f
}

# Tables to <outdir>/tables/stream_model_<name>.tsv and figures to <outdir>/figures/.
stm_write <- function(tabs, figs, outdir, dpi = 150) {
  dir.create(file.path(outdir, "tables"), recursive = TRUE, showWarnings = FALSE); dir.create(file.path(outdir, "figures"), recursive = TRUE, showWarnings = FALSE)
  for (nm in names(tabs)) if (is.data.frame(tabs[[nm]]) && nrow(tabs[[nm]])) fwrite(tabs[[nm]], file.path(outdir, "tables", paste0("stream_model_", nm, ".tsv")), sep = "\t")
  for (nm in names(figs)) { sz <- attr(figs[[nm]], "size") %||% c(12, 4.5); ggplot2::ggsave(file.path(outdir, "figures", paste0(nm, ".png")), figs[[nm]], width = sz[1L], height = sz[2L], dpi = dpi, bg = "white") }
  invisible(names(tabs))
}
