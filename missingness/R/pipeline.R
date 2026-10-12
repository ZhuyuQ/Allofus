# ==============================================================================
# pipeline.R
# ------------------------------------------------------------------------------
# The two orchestrating functions. They only CALL the functions in the other files
# (no statistics are written here), so the order of the analysis can be read off
# in one place.
#
#   run_etl()                  raw exports (L0) -> compact day tables + persons (L1), run once
#   run_missingness_pipeline() day tables + persons -> tables, figures, metadata (L2)
# ==============================================================================

# ---- ETL: raw exports -> L1 --------------------------------------------------------------------------

# Read every cohort's raw export for each stream, build the day table, stack the
# cohorts, apply the "mostly NA / constant" column rule ONCE on the stacked table and
# write it. Each stream is released as soon as it is on disk (only its first-record dates and a
# short summary are kept), so the ETL needs the memory of ONE stream. Returns
# list(persons, dropped_columns, summary); read the tables back with open_day_store().
# `lister`/`reader` are injectable for tests.
run_etl <- function(cfg, out_dir, lister = gsutil_ls, reader = gsutil_fread,
                    final_pid_path, cci_path) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # people in more than one cohort: resolve ONCE, so every person counts in exactly one cohort
  ov <- resolve_cohort_overlap(read_final_pid(final_pid_path), cfg$overlap$policy, cfg$overlap$priority)
  log_info("cohorts: %d people in final_pid, %d belong to more than one cohort; policy '%s' keeps %d",
           ov$n_unique_people, length(ov$overlap_ids), ov$policy, ov$n_kept)
  if (nrow(ov$overlap) > 0L) print(ov$overlap)
  fwrite(ov$overlap, file.path(out_dir, "cohort_overlap.tsv"), sep = "\t")
  ids_by_cohort <- split(ov$assignment$person_id, ov$assignment$cohort)
  keep_cohort <- function(raw, co) raw[person_id %in% (ids_by_cohort[[co]] %||% numeric())]

  first_dates <- list(); summary <- list(); dropped <- list()
  for (stream in names(cfg$streams)) {
    spec <- cfg$streams[[stream]]
    per_cohort <- lapply(cfg$cohorts, function(co) {
      raw <- keep_cohort(read_export_files(export_path(cfg, co, stream), spec$col_classes,
                                           lister = lister, reader = reader), co)
      if (nrow(raw) == 0L) return(NULL)
      build_day_table(raw, stream, cfg, apply_drop_rule = FALSE)
    })
    day <- rbindlist(per_cohort, use.names = TRUE, fill = TRUE)
    assert_unique_keys(day, c("person_id", "date"), paste(stream, "stacked day table"))
    day <- drop_uninformative_columns(day, protect = c("person_id", "date", "has_main_sleep", "has_nap", "n_records"),
                                      na_frac = cfg$etl$na_drop_frac)
    setkey(day, person_id, date)
    write_day_table(day, out_dir, stream, cfg$storage_format)
    log_info("ETL %s: %d rows, %d people, dropped columns: %s", stream, nrow(day), uniqueN(day$person_id),
             paste(attr(day, "dropped_columns")$column, collapse = ", "))
    first_dates[[stream]] <- day[, .(date = min(date)), by = person_id]
    summary[[stream]] <- data.table(stream = stream, n_rows = nrow(day), n_persons = uniqueN(day$person_id),
                                    n_columns = ncol(day), columns = paste(names(day), collapse = ","), first_date = as.character(min(day$date)),
                                    last_date = as.character(max(day$date)))
    dropped[[stream]] <- cbind(stream = stream, attr(day, "dropped_columns"))
    rm(day, per_cohort); invisible(gc())
  }

  # persons: meta of all cohorts + cohort assignment + CCI
  meta <- rbindlist(lapply(cfg$cohorts, function(co) {
    read_export_files(export_path(cfg, co, "meta"),
                      c(race = "character", ethnicity = "character", sex_at_birth = "character",
                        self_reported_category = "character"),
                      lister = lister, reader = reader)
  }), use.names = TRUE, fill = TRUE)
  meta <- unique(meta, by = "person_id")          # a person exported for two cohorts appears once
  persons <- build_person_table(meta, ov$assignment, read_cci(cci_path),
                                first_record_dates(first_dates), cfg)
  saveRDS(persons, file.path(out_dir, "persons.rds"))
  summary <- rbindlist(summary); dropped <- rbindlist(dropped, fill = TRUE)
  fwrite(summary, file.path(out_dir, "manifest.tsv"), sep = "\t")
  fwrite(dropped, file.path(out_dir, "dropped_columns.tsv"), sep = "\t")
  invisible(list(persons = persons, dropped_columns = dropped, summary = summary))
}

# Open L1 (written by run_etl()). Nothing big is read here: `day_tables` is a day store that reads
# one stream at a time when the pipeline asks for it, and `persons` is the (small) person table.
load_inputs <- function(cfg, dir) {
  store <- open_day_store(cfg, dir)
  list(day_tables = store, persons = store$persons())
}

# ---- analysis: L1 -> L2 --------------------------------------------------------------------------------------------

# Add cohort (and stream) to a person-level table.
with_cohort <- function(x, persons, stream) {
  out <- merge(x, persons[, .(person_id, cohort)], by = "person_id", all.x = TRUE)
  out[, stream := stream]
  out[]
}

# Day-table columns the onset / sensitivity models need for one stream: the predictors' source
# columns and the variables of the partial-wear rule (nothing else is read from disk).
model_columns <- function(cfg, stream) {
  preds <- if (stream == "activity") cfg$activity_predictors else cfg$own_predictors[[stream]]
  unique(c(vapply(cfg$predictors[preds], function(p) p$col, ""), all.vars(str2lang(cfg$onset$partial_wear[[stream]]))))
}

# The day tables a model stage of `stream` works on: the stream itself and (for sleep / hr) the
# activity table, each restricted to the needed columns.
stream_inputs <- function(store, stream, cfg, extra_cols = character()) {
  out <- setNames(list(store$get(stream, unique(c(model_columns(cfg, stream), extra_cols)))), stream)
  if (stream != "activity") out$activity <- store$get("activity", model_columns(cfg, "activity"))
  out
}

# stages: any of "missingness", "temporal", "descriptive", "onset", "gap_models", "sensitivity", "missing_figures"
#   "missing_figures" = only the figures / tables of the missingness summaries that are built from person-level
#   counts (missing % by race / sex / CCI / age group, missing share by study week / year / weekday / age): no
#   model, no heatmaps, no descriptive metrics - a quick way to redraw them.
# day_tables: the day store from load_inputs() (or a named list of day tables)
# model_streams: which streams the model stages (onset, gap_models, sensitivity) run on. Streams are
#   processed one at a time; every finished model is a checkpoint in <outdir>/checkpoints, so a call
#   can be repeated after a crash and only does the work that is missing.
run_missingness_pipeline <- function(day_tables, persons, cfg, outdir,
                                     stages = c("missingness", "temporal", "descriptive",
                                                "onset", "gap_models", "sensitivity"),
                                     onset_roles = c("primary", "secondary_threshold", "exploratory"),
                                     model_streams = NULL, save_figures = TRUE) {
  set.seed(cfg$seed)
  store   <- as_day_store(day_tables, persons)
  streams <- store$streams
  model_streams <- model_streams %||% streams
  bad <- setdiff(model_streams, streams)
  if (length(bad)) stop("unknown stream(s) in model_streams: ", paste(bad, collapse = ", "), call. = FALSE)
  ckpt <- file.path(outdir, cfg$checkpoint$dirname)
  res <- list(); figs <- character(); units <- list()
  has <- function(s) s %in% stages
  fig <- function(name, p, w = 7, h = 5) {
    if (!save_figures) return(invisible(NULL))
    if (!is.null(attr(p, "size"))) { w <- attr(p, "size")[1L]; h <- attr(p, "size")[2L] }     # violin figures with many groups set their own size
    save_plot(p, file.path(outdir, "figures", paste0(name, ".png")), w, h, cfg)
    figs <<- c(figs, name)
  }

  # ---- 0-1. person-level summaries, ONE stream at a time (only a stream's own table is in memory) ----
  if (any(c("missingness", "temporal", "descriptive", "missing_figures") %in% stages)) {
    metric_sets <- descriptive_metrics()
    sel <- list(); gaps <- list(); pg_l <- list(); long_l <- list()
    for (s in streams) {
      d <- store$get(s, cols = if (has("descriptive")) NULL else character())     # NULL = all columns
      gaps[[s]] <- get_gap_runs(d)
      pg_l[[s]] <- with_cohort(summarise_person_gaps(d, gaps[[s]]), persons, s)
      sel[[s]]  <- as.data.table(c(stream = s, check_selection_rule(d, s, cfg)))
      if (has("descriptive")) long_l[[s]] <- person_means_long(person_metric_stats(d, metric_sets[[s]]), persons, s)
      if (has("missingness")) fig(paste0("heatmap_", s), plot_missing_heatmap(d, cfg, title = paste("Observed vs Missing Days:", s)), 7, 5)
      rm(d); invisible(gc())
    }
    pg <- rbindlist(pg_l)
    gaps_all <- rbindlist(lapply(streams, function(s) with_cohort(gaps[[s]], persons, s)))
    if (any(c("missingness", "temporal", "descriptive") %in% stages)) {
      res$selection_check <- rbindlist(sel)
      res$person_missingness <- pg
      res$gap_runs <- gaps_all
    }
  }

  # ---- 2. cohort summary and comparisons ----
  if (has("missingness")) {
    res$cohort_missingness <- cohort_missingness_table(pg)
    vars <- c("pct_missing", "gap_max", "n_gaps")
    cmp <- unlist(lapply(streams, function(s) lapply(vars, function(v) compare_cohorts(pg, v, s))), recursive = FALSE)
    res$cohort_tests_overall  <- rbindlist(lapply(cmp, `[[`, "overall"))
    res$cohort_tests_pairwise <- rbindlist(lapply(cmp, `[[`, "pairwise"))
    res$missing_fraction_models <- rbindlist(lapply(streams, function(s) {
      pgs <- pg[stream == s]
      rbind(cbind(stream = s, fit_missing_fraction(pgs, persons, FALSE)),
            cbind(stream = s, fit_missing_fraction(pgs, persons, TRUE)))
    }), fill = TRUE)
    long_pct <- pg[, .(person_id, cohort, metric = stream, value = pct_missing)]
    long_gap <- pg[, .(person_id, cohort, metric = stream, value = gap_max)]
    fig("missing_pct_by_cohort", plot_metric_violin(long_pct, cfg, "Person-Level % of Days Missing by Cohort", "% days missing", ncol = 3L), 10, 5)
    fig("longest_gap_by_cohort", plot_metric_violin(long_gap, cfg, "Person-Level Longest Gap by Cohort", "Days", ncol = 3L), 10, 5)
    fig("gap_length_distribution", plot_gap_length_distribution(gaps_all, cfg), 11, 5)
  }

  # ---- 2b. person-level % missing by race / sex / CCI / age group, within each cohort ----
  if (any(c("missingness", "missing_figures") %in% stages)) {
    fig("gap_ecdf", plot_gap_ecdf(pg, cfg), 7, 9)
    ms <- missing_pct_by_stratum(pg, persons)
    res$missing_pct_by_stratum       <- ms$summary
    res$missing_pct_by_stratum_tests <- ms$tests
    res$missing_pct_by_stratum_models <- ms$models
    for (st in c("sex", "age_group", "race_plot", "cci_group")) {
      fig(paste0("missing_pct_by_", st),
          plot_missing_pct_by_stratum(ms$long, st, ms$tests, cfg,
                                      paste("Person-Level % of Days Missing by Cohort and", stratum_label(st))), 13, 5.5)
    }
  }

  # ---- 3. temporal structure ----
  if (has("temporal")) {
    counts <- rbindlist(lapply(streams, function(s) transition_person_counts(pg[stream == s])))
    res$transition <- transition_table(counts, cfg)
    res <- c(res, missing_rate_tables(pg, gaps, persons, streams, cfg))
    dw <- rbindlist(lapply(streams, function(s) cbind(stream = s, dow_counts(pg[stream == s], gaps[[s]]))))
    dow_fits <- lapply(streams, function(s) lapply(c(FALSE, TRUE), function(adj)
      fit_dow_model(dw[stream == s], persons, cfg, adjusted = adj)))
    res$weekday_tests <- rbindlist(lapply(seq_along(streams), function(i)
      rbindlist(lapply(dow_fits[[i]], function(f) cbind(stream = streams[i], f$tests)))))
    res$weekday_odds_ratios <- rbindlist(lapply(seq_along(streams), function(i)
      rbindlist(lapply(dow_fits[[i]], function(f) cbind(stream = streams[i], f$coefs)))))
    fig("transition", plot_transition(res$transition, cfg), 10, 6)
  }

  # ---- 3b. missing share by study week / calendar year / weekday / attained age (one line + shaded band per cohort) ----
  if (any(c("temporal", "missing_figures") %in% stages)) {
    if (!has("temporal")) res <- c(res, missing_rate_tables(pg, gaps, persons, streams, cfg))
    fig("missing_by_study_week", plot_missing_over_time(res$study_week_rates, "week", cfg, "Study week", "Missing Days by Study Week"), 8, 9)
    fig("missing_by_year", plot_missing_over_time(res$calendar_year_rates, "year", cfg, "Calendar year", "Missing Days by Calendar Year"), 8, 9)
    fig("missing_by_weekday", plot_missing_by_weekday(res$weekday_rates, cfg), 8, 9)
    fig("missing_by_age", plot_missing_by_age(res$age_rates, cfg), 8, 9)
  }

  # ---- 4. descriptive statistics of the metric values ----
  if (has("descriptive")) {
    long_all <- rbindlist(long_l)
    res$metric_summary <- suppress_small_cells(metric_summary_table(long_all), cfg)
    for (st in c("sex", "age_group", "race_plot", "cci_group")) {
      res[[paste0("metric_summary_by_", st)]] <- suppress_small_cells(metric_summary_table(long_all, st), cfg)
    }
    for (s in streams) {
      l <- long_all[stream == s]
      fig(paste0("metrics_", s, "_by_cohort"), plot_metric_violin(l, cfg, paste("Patient-Level Mean", nice_label(s), "Metrics by Cohort")), 9, 8)
      for (st in c("sex", "age_group", "race_plot", "cci_group")) {
        fig(paste0("metrics_", s, "_by_", st), plot_metric_by_stratum(l, st, cfg, paste("Mean", nice_label(s), "Metrics by Cohort and", nice_label(st))), 11, 9)
      }
    }
  }
  if (exists("pg")) rm(list = intersect(c("pg", "gaps", "gaps_all", "pg_l", "long_l"), ls()))
  invisible(gc())

  # ---- 5-7. model stages: ONE STREAM AT A TIME, every model its own checkpointed unit ----
  # A stream's tables are built only when some model of that stream really has to run (a stream whose
  # models are all finished is skipped without reading its data) and are dropped before the next stream.
  # Input fingerprints tie every checkpoint to the exact L1 files, so changed data never reuse old fits.
  dkeys <- setNames(lapply(streams, function(s) data_key(store, s)), streams)
  grid  <- make_onset_grid(cfg)[role %in% onset_roles]

  # ---- 5. gap-onset models (primary, A, D) ----
  if (has("onset")) {
    for (s in model_streams) {
      inputs <- lazy(function() stream_inputs(store, s, cfg))
      bundle <- lazy(function() prepare_stream_bundle(s, inputs$get(), persons, cfg))
      units[[paste0("onset:", s)]] <- run_onset_stage(s, bundle, grid, cfg, ckpt, dkeys[[s]])
      bundle$release(); inputs$release()
    }
    res <- c(res, collect_onset_results(ckpt, cfg, dkeys))
    if (!is.null(res$onset_coefficients) && nrow(res$onset_coefficients) > 0L) {
      fig("onset_primary_forest", plot_onset_forest(res$onset_coefficients, cfg), 8, 7)
      d_tab <- res$onset_coefficients[role %in% c("primary", "secondary_threshold") & term == cfg$primary_term,
                                      .(stream, label = paste0("gap >= ", k, " day", ifelse(k > 1, "s", "")), or, or_low, or_high)]
      if (nrow(d_tab) > 0L) fig("forest_gap_definition", plot_forest(d_tab, "Odds of a Gap Starting Tomorrow per 1 SD Below Usual Steps, by Gap Definition", group_col = "stream", style = cfg$plot$interval_style), 10, 5)
    }
    if (!is.null(res$risk_unadjusted) && "risk" %in% names(res$risk_unadjusted)) fig("risk_curve", plot_risk_curve(res$risk_unadjusted, cfg), 7, 9)
    if (nrow(res$cohort_specific) > 0L) fig("cohort_forest", plot_cohort_forest(res$cohort_specific, res$cohort_interaction, cfg), 10, 5)
  }

  # ---- 6. gap-level models (B and C): built from the gaps + one activity predictor, nothing else ----
  if (has("gap_models")) {
    act_col <- cfg$predictors[["act_steps"]]$col
    for (s in model_streams) {
      gap_data <- lazy(function() {
        days <- store$get(s, cols = if (s == "activity") act_col else character())
        act  <- if (s == "activity") days else store$get("activity", cols = act_col)
        build_gap_level_data(days, act, persons, cfg)
      })
      units[[paste0("gap_models:", s)]] <- run_gap_stage(s, gap_data, cfg, ckpt, dkeys[[s]])
      gap_data$release()
    }
    res <- c(res, collect_gap_results(ckpt, cfg, dkeys))
    if (nrow(res$model_B_curve) > 0L) fig("gap_dose_response", plot_gap_dose_response(res$model_B_curve, res$model_B_observed, cfg), 8, 9)
    if (nrow(res$model_C_predicted_probs) > 0L) fig("gap_length_probs", plot_ordinal_probs(res$model_C_predicted_probs, cfg), 10, 4.5)
  }

  # ---- 7. sensitivity analyses of the primary model ----
  if (has("sensitivity")) {
    prim_specs <- make_onset_grid(cfg)[role == "primary" & stream %in% model_streams]
    for (i in seq_len(nrow(prim_specs))) {
      s <- prim_specs$stream[i]
      inputs <- lazy(function() stream_inputs(store, s, cfg, extra_cols = if (s == "sleep") "has_main_sleep"))
      bundle <- lazy(function() prepare_stream_bundle(s, inputs$get(), persons, cfg, trailing = TRUE))   # tw_ columns only needed here
      alt <- if (s == "sleep" && "has_main_sleep" %in% store$columns("sleep"))
        lazy(function() prepare_stream_bundle("sleep", inputs$get(), persons, cfg, observed_filter = "has_main_sleep")) else NULL
      units[[paste0("sensitivity:", s)]] <- run_sensitivity_stage(s, bundle, alt, prim_specs[i], cfg, ckpt, dkeys[[s]])
      bundle$release(); if (!is.null(alt)) alt$release(); inputs$release()
    }
    res <- c(res, collect_sensitivity_results(ckpt, cfg, dkeys))
    if (nrow(res$sensitivity) > 0L) {
      ord <- c("main", "trailing_baseline", "exclude_partial_wear", "first_gap_only", "main_sleep_only", "conditional_logit", "cluster_robust_logit")
      f_tab <- res$sensitivity[!is.na(or), .(stream, label = analysis, or, or_low, or_high)][order(match(label, ord))]
      if (nrow(f_tab) > 0L) fig("forest_sensitivity", plot_forest(f_tab, "Primary Effect Under Different Assumptions", group_col = "stream", style = cfg$plot$interval_style), 10, 5)
    }
  }

  # ---- 8. write everything, and say plainly what ran, what was resumed and what failed ----
  res$unit_run_status  <- unit_status_table(units)
  res$checkpoint_status <- checkpoint_report(ckpt, cfg, dkeys)
  res$unit_failures <- rbindlist(list(res$onset_failures, res$sensitivity_failures, res$gap_failures), fill = TRUE)
  if (!is.null(res$unit_failures) && nrow(res$unit_failures) > 0L) {
    warning(sprintf("%d model unit(s) FAILED - see unit_failures.tsv; run the same cell again to retry them.", nrow(res$unit_failures)), call. = FALSE)
    print(res$unit_failures)
  }
  write_tables(res, outdir)
  if (length(figs)) write_figure_captions(figs, outdir)
  if (any(c("onset", "gap_models", "sensitivity") %in% stages)) write_results_summary(res, cfg, outdir)   # model results only
  if (!identical(stages, "missing_figures")) write_run_metadata(cfg, outdir, store, extra = list(stages = stages, onset_roles = onset_roles))
  invisible(res)
}
