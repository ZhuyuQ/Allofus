# ==============================================================================
# compare_bands.R
# ------------------------------------------------------------------------------
# Side-by-side comparison of two FINISHED runs of the missingness pipeline on different people (for example compliance >= 70%
# and 50-60%). It only reads the saved result tables of each run (<results>/tables/*.tsv) and, for the note on the "+1 SD"
# scale, the day tables of each run; no model is refitted. Figures are drawn in the style of the existing ones; the statistics
# are in tables outside the plots.
#
#   out_dirs : named list, band label -> results folder (the first band is the reference), e.g.
#              list(">=70%" = ".../results", "50-60%" = ".../results/5060")
#   stores   : named list, band label -> day store (open_day_store(cfg, L1 folder)); only for compare_within_sd
#
# Statistics (tables, never printed on the figures)
#   * The two groups are different people (the bands cannot overlap), so the two estimates of a term are independent:
#     log ratio of odds ratios = beta_cmp - beta_ref, SE = sqrt(se_ref^2 + se_cmp^2), ratio = exp(log ratio) with a 95% interval.
#     No p-value is given (descriptive comparison); a ratio whose interval contains 1 is not evidence that the groups are the same.
#   * "+1 SD" is each model's own: every run standardizes with its own pooled within-person SD of log(1 + steps), so the two odds
#     ratios are not per the same raw change. The SDs are recomputed from the day tables and written to compare_within_sd.tsv and
#     the figure caption; the estimates are plotted as the models gave them.
#   * Missingness figures are descriptive (the % of days missing differs by design: compliance defines the groups).
# ==============================================================================

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

CMP_BAND_COLORS   <- c("#2C3E50", "#E67E22")                     # first (reference) band, second band
CMP_STREAM_LABELS <- c(activity = "Activity", sleep = "Sleep", hr = "Heart rate")
cmp_colors <- function(bands) stats::setNames(CMP_BAND_COLORS[seq_along(bands)], bands)

# Table `name` of every band, stacked, with the column `band` (a factor in the order of out_dirs).
cmp_read <- function(out_dirs, name) {
  parts <- lapply(names(out_dirs), function(b) {
    f <- file.path(out_dirs[[b]], "tables", paste0(name, ".tsv"))
    if (!file.exists(f)) stop(sprintf("table '%s' of band '%s' not found: %s (run the stage that writes it first)", name, b, f), call. = FALSE)
    cbind(band = b, data.table::fread(f))
  })
  out <- rbindlist(parts, fill = TRUE)
  out[, band := factor(band, levels = names(out_dirs))]
  out[]
}

# ---- the primary effect -------------------------------------------------------------------------------------------------

# The primary term of the primary model (joint, adjusted, gap >= 1 day) per stream and band, with the size of the model.
cmp_primary <- function(out_dirs, cfg) {
  oc <- cmp_read(out_dirs, "onset_coefficients"); om <- cmp_read(out_dirs, "onset_models")
  p <- oc[role == "primary" & cohort_stratum == "all" & term == cfg$primary_term, .(band, stream, estimate, se, or, or_low, or_high, status)]
  n <- om[role == "primary" & cohort_stratum == "all", .(band, stream, n_rows, n_persons, n_events, n_persons_event)]
  out <- merge(p, n, by = c("band", "stream"), all.x = TRUE)
  out[, stream := factor(stream, levels = names(CMP_STREAM_LABELS))]
  setorder(out, stream, band)
  out[, stream := as.character(stream)]
  out[]
}

# The pooled within-person SD of log(1 + steps) that "+1 SD" means in each run (the activity predictor as the models of each
# stream see it: the sleep and heart-rate models use the activity days that also have a record in that stream).
# store: a day store of one run. Returns stream, sd_within.
compare_within_sd <- function(store, cfg) {
  act <- store$get("activity", cfg$predictors$act_steps$col)
  rbindlist(lapply(names(CMP_STREAM_LABELS), function(s) {
    a <- if (s == "activity") act else restrict_to_stream_days(act, store$get(s, character()))
    sc <- attr(predictor_table(a, "act_steps", cfg, keep_x = FALSE), "scales")
    data.table(stream = s, sd_within = sc$sd_within[1L])
  }))
}

# Ratio of the two odds ratios (second band relative to the first) per stream; independent groups, see the header.
cmp_ror <- function(prim, bands, sds = NULL) {
  a <- prim[band == bands[1L]]; b <- prim[band == bands[2L]]
  m <- merge(a, b, by = "stream", suffixes = c("_ref", "_cmp"))
  m[, c("band_ref", "band_cmp") := NULL]
  m[, `:=`(log_ratio = estimate_cmp - estimate_ref, se_log_ratio = sqrt(se_ref^2 + se_cmp^2))]
  m[, `:=`(ratio = exp(log_ratio), ratio_low = exp(log_ratio - 1.96 * se_log_ratio), ratio_high = exp(log_ratio + 1.96 * se_log_ratio),
           reference = bands[1L], compared = bands[2L])]
  if (!is.null(sds)) {
    pick <- function(b) { x <- sds[band == b]; x$sd_within[match(m$stream, x$stream)] }
    m[, `:=`(within_sd_ref = pick(bands[1L]), within_sd_cmp = pick(bands[2L]))]
  }
  m[, stream := factor(stream, levels = names(CMP_STREAM_LABELS))]
  setorder(m, stream); m[, stream := as.character(stream)]
  m[]
}

# One forest: the primary term per stream, the bands in colour. sds (optional): band, stream, sd_within, for the caption.
plot_cmp_primary <- function(prim, cfg, sds = NULL) {
  bands <- levels(prim$band); style <- cfg$plot$interval_style %||% "se"
  d <- copy(prim)[!is.na(or)]; d <- d[drawable_interval(d, "or", "or_low", "or_high", log = TRUE)]
  d[, stream_lab := factor(CMP_STREAM_LABELS[stream], levels = rev(CMP_STREAM_LABELS))]
  add_interval_columns(d, "or", "or_low", "or_high", log = TRUE, style = style)
  cap <- interval_caption(style)
  if (!is.null(sds) && nrow(sds)) {
    one <- function(b) paste0(b, ": ", paste(sprintf("%.2f", sds[band == b]$sd_within[match(names(CMP_STREAM_LABELS), sds[band == b]$stream)]), collapse = " / "))
    cap <- paste0(cap, "\nEach estimate is per +1 SD of its own model: within-person SD of log(1 + steps), activity / sleep / heart rate: ", paste(vapply(bands, one, ""), collapse = ";  "))
  }
  ggplot(d, aes(x = est_, y = stream_lab)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
    interval_layer(style, horizontal = TRUE, color_var = "band", position = position_dodge(width = 0.6)) +
    scale_x_log10(breaks = OR_BREAKS) + expand_limits(x = 1) +
    scale_color_manual(values = cmp_colors(bands), breaks = bands) + scale_fill_manual(values = cmp_colors(bands), breaks = bands) +
    labs(title = "Primary Effect by Compliance Group",
         subtitle = ov_wrap("Odds that a gap starts tomorrow per +1 SD higher steps on the last full day than the person's usual (joint adjusted model; OR < 1: a more active day is followed by fewer gaps)", 120),
         x = "Odds ratio per +1 SD (log scale)", y = NULL, color = "Compliance", fill = "Compliance", caption = cap) +
    theme_bw(base_size = 11) + theme(legend.position = "bottom", plot.subtitle = element_text(size = 8.5))
}

# ---- the sensitivity analyses of both bands ------------------------------------------------------------------------------

# One panel per stream, one row per analysis (the tested term of each), the bands in colour.
plot_cmp_sensitivity <- function(sens, cfg) {
  bands <- levels(sens$band); style <- cfg$plot$interval_style %||% "se"
  d <- copy(sens)[!is.na(or)]; d <- d[drawable_interval(d, "or", "or_low", "or_high", log = TRUE)]
  ord <- c(intersect(SENSITIVITY_ORDER, unique(d$analysis)), sort(setdiff(unique(d$analysis), SENSITIVITY_ORDER)))
  lab <- function(x) { o <- nice_label(x); paste0(toupper(substr(o, 1L, 1L)), substring(o, 2L)) }          # nice_label leaves a leading "both" in lower case
  d[, `:=`(label = factor(lab(analysis), levels = rev(lab(ord))),
           stream_lab = factor(CMP_STREAM_LABELS[stream], levels = CMP_STREAM_LABELS))]
  add_interval_columns(d, "or", "or_low", "or_high", log = TRUE, style = style)
  p <- ggplot(d, aes(x = est_, y = label)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
    interval_layer(style, horizontal = TRUE, color_var = "band", position = position_dodge(width = 0.7)) +
    facet_wrap(~stream_lab, nrow = 1) +
    scale_x_log10(breaks = OR_BREAKS) + expand_limits(x = 1) +
    scale_color_manual(values = cmp_colors(bands), breaks = bands) + scale_fill_manual(values = cmp_colors(bands), breaks = bands) +
    labs(title = "Sensitivity Analyses by Compliance Group",
         subtitle = ov_wrap("The tested within-person steps term of every sensitivity analysis (odds ratio per +1 SD of each model's own SD); the likelihood-ratio analyses (spline, random slope) have no odds ratio and are not drawn", 140),
         x = "Odds ratio per +1 SD (log scale)", y = NULL, color = "Compliance", fill = "Compliance", caption = interval_caption(style)) +
    theme_bw(base_size = 11) + theme(legend.position = "bottom", strip.text = element_text(face = "bold"), plot.subtitle = element_text(size = 8.5))
  attr(p, "size") <- c(13, max(4.5, 0.42 * length(ord) + 2.2))
  p
}

# ---- missingness ------------------------------------------------------------------------------------------------------------

# Person-level variable by cohort, the bands side by side (violin + box), one panel per stream.
plot_cmp_person <- function(pm, var, ylab, title, cfg) {
  bands <- levels(pm$band)
  d <- copy(pm)[!is.na(get(var))]
  d[, `:=`(cohort = factor(cohort, levels = cfg$cohorts), stream = factor(stream, levels = names(CMP_STREAM_LABELS), labels = CMP_STREAM_LABELS))]
  ggplot(d, aes(x = cohort, y = .data[[var]], fill = band)) +
    geom_violin(position = position_dodge(width = 0.85), scale = "width", alpha = 0.55, color = NA) +
    geom_boxplot(aes(group = interaction(cohort, band)), position = position_dodge(width = 0.85), width = 0.14, outlier.shape = NA,
                 fill = "white", color = "grey25", linewidth = 0.3) +
    facet_wrap(~stream, nrow = 1, scales = "free_y") +
    scale_fill_manual(values = cmp_colors(bands), breaks = bands) +
    labs(title = title, x = NULL, y = ylab, fill = "Compliance",
         caption = "One value per person. Descriptive: the % of days missing differs by design (compliance defines the groups).") +
    theme_violin_panel() + theme(legend.position = "bottom")
}

# Share of all gaps in each length group, by cohort; one column of panels per band, one row per stream.
plot_cmp_gap_length <- function(gaps, cfg) {
  d <- copy(gaps)[, `:=`(cohort = factor(cohort, levels = cfg$cohorts), len = gap_length_bin(gap_length, cfg),
                         stream = factor(stream, levels = names(CMP_STREAM_LABELS), labels = CMP_STREAM_LABELS))]
  s <- d[, .N, by = .(band, stream, cohort, len)][, share := 100 * N / sum(N), by = .(band, stream, cohort)]
  ggplot(s, aes(x = len, y = share, fill = cohort)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.65) +
    facet_grid(stream ~ band) +
    scale_fill_manual(values = cfg$palette$cohort) +
    scale_y_continuous(labels = function(x) paste0(x, "%"), expand = expansion(mult = c(0, 0.08))) +
    labs(title = "Gap Length Distribution by Compliance Group", x = "Gap length (days)", y = "% of gaps (within stream, cohort and group)", fill = "Cohort") +
    theme_minimal_top()
}

# Onset and continuation probabilities (the `transition` table) by cohort, the bands in colour.
plot_cmp_transition <- function(tr, cfg) {
  bands <- levels(tr$band); style <- cfg$plot$transition_style %||% "box"
  long <- rbindlist(list(
    tr[, .(band, stream, cohort, measure = "P(missing tomorrow | observed today)", est = p_onset, lo = p_onset_low, hi = p_onset_high)],
    tr[, .(band, stream, cohort, measure = "P(missing tomorrow | missing today)", est = p_continue, lo = p_continue_low, hi = p_continue_high)]))
  long[, `:=`(cohort = factor(cohort, levels = cfg$cohorts), stream = factor(stream, levels = names(CMP_STREAM_LABELS), labels = CMP_STREAM_LABELS))]
  add_interval_columns(long, "est", "lo", "hi", style = style)
  ggplot(long, aes(x = cohort, y = est_, color = band, fill = band)) +
    interval_layer(style, horizontal = FALSE, color_var = "band", position = position_dodge(width = 0.7), show.legend = TRUE) +
    facet_grid(stream ~ measure, scales = "free_y", labeller = labeller(measure = label_wrap_gen(30))) +
    scale_color_manual(values = cmp_colors(bands), breaks = bands) + scale_fill_manual(values = cmp_colors(bands), breaks = bands) +
    labs(title = "Do Missing Days Come in Runs? By Compliance Group", x = NULL, y = "Probability", color = "Compliance", fill = "Compliance",
         caption = interval_caption(style)) +
    theme_violin_panel() + theme(legend.position = "bottom")
}

# ---- everything ---------------------------------------------------------------------------------------------------------------

# Reads the saved tables of both runs, writes the tables to <cmp_dir>/tables and the figures to <cmp_dir>/figures.
# Needs from each run: stages "missingness" and "temporal" (missingness tables), "onset" with the primary model and "sensitivity".
# stores (optional): day stores of the runs for the within-person SD; without them the SD note is left out.
cmp_run <- function(out_dirs, cfg, cmp_dir, stores = NULL) {
  if (length(out_dirs) != 2L) stop("exactly two bands are compared", call. = FALSE)
  bands <- names(out_dirs)
  td <- file.path(cmp_dir, "tables"); fd <- file.path(cmp_dir, "figures")
  dir.create(td, recursive = TRUE, showWarnings = FALSE); dir.create(fd, recursive = TRUE, showWarnings = FALSE)
  w <- function(x, name) { data.table::fwrite(x, file.path(td, paste0(name, ".tsv")), sep = "\t"); x }
  fig <- function(p, name, w_ = 9, h_ = 5) {
    if (!is.null(attr(p, "size"))) { w_ <- attr(p, "size")[1L]; h_ <- attr(p, "size")[2L] }
    save_plot(p, file.path(fd, paste0(name, ".png")), w_, h_, cfg)
  }
  res <- list()

  sds <- NULL
  if (!is.null(stores)) {
    sds <- rbindlist(lapply(bands, function(b) cbind(band = b, compare_within_sd(stores[[b]], cfg))))
    sds[, band := factor(band, levels = bands)]
    res$within_sd <- w(sds, "compare_within_sd")
  }
  res$primary <- w(cmp_primary(out_dirs, cfg), "compare_primary")
  res$ratio   <- w(cmp_ror(res$primary, bands, sds), "compare_primary_ratio")
  fig(plot_cmp_primary(res$primary, cfg, sds), "compare_primary_effect", 9, 4.2)

  res$sensitivity <- w(cmp_read(out_dirs, "sensitivity"), "compare_sensitivity")
  fig(plot_cmp_sensitivity(res$sensitivity, cfg), "compare_sensitivity")

  pm <- cmp_read(out_dirs, "person_missingness")
  fig(plot_cmp_person(pm, "pct_missing", "% of days missing", "Person-Level % of Days Missing by Compliance Group", cfg), "compare_missing_pct", 11, 5)
  fig(plot_cmp_person(pm, "gap_max", "Longest gap (days)", "Person-Level Longest Gap by Compliance Group", cfg), "compare_longest_gap", 11, 5)
  fig(plot_cmp_person(pm, "n_gaps", "Gaps per person", "Gaps per Person by Compliance Group", cfg), "compare_gaps_per_person", 11, 5)
  fig(plot_cmp_gap_length(cmp_read(out_dirs, "gap_runs"), cfg), "compare_gap_length_distribution", 11, 8)
  fig(plot_cmp_transition(cmp_read(out_dirs, "transition"), cfg), "compare_transition", 11, 8)
  cm <- cmp_read(out_dirs, "cohort_missingness")
  cm[, stream := factor(stream, levels = names(CMP_STREAM_LABELS))]; setorder(cm, stream, cohort, band); cm[, stream := as.character(stream)]
  res$missingness <- w(cm, "compare_missingness")
  invisible(res)
}
