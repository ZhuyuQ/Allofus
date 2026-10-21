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
  # the optional Fitbit device information, if its exports exist (Cell 6c); otherwise run_etl_device() can be run later on its own
  tryCatch(run_etl_device(cfg, out_dir, lister, reader, persons = persons),
           error = function(e) log_info("device information not built (%s)", sub("\\..*$", ".", conditionMessage(e))))
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

# stages: any of "missingness", "temporal", "descriptive", "runin_description", "onset", "gap_models", "sensitivity",
#   "runin", "runin_sensitivity", "missing_figures", "device"
#   "variables"        = person-level % of days missing against every model variable (cohort, sex, watch version, age, between-person
#                        predictor levels): figures outcome_by_variable_*, tables outcome_variable_*. Needs Cell 8b for the watch version.
#   "device"           = missing share and composition of the Fitbit device version by cohort (needs Cells 6c and 8b); no model.
#   "runin_description" = the run-in (first cfg$runin$days days of follow-up): distribution of missing days, the observed gap rate after
#                        the run-in by run-in group (H1), and the steps-gap curve by run-in group (H2); no model.
#   "runin"            = that description first, then the primary onset model and models B and C on the days / gaps AFTER the run-in,
#                        without a run-in term and with the run-in and its interaction (see runin.R).
#   "runin_sensitivity" = everything of "runin", plus the other forms of the run-in (main effect only, categories, linear), with figures of their own.
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
  if (inherits(store, "day_store") && identical(store$kind, "files") && is.null(store$device_person))      # a store opened by an older version of the code: open it again
    store <- open_day_store(cfg, store$dir)
  # the main device family of each person (device.R) joins the person table: it is a covariate of every model (replaces calendar year)
  dp_store <- store$device_person()
  if (!is.null(dp_store))
    persons <- data.table::setkey(merge(data.table::as.data.table(persons), data.table::as.data.table(dp_store)[, .(person_id, device_family)],
                                        by = "person_id", all.x = TRUE, sort = FALSE), person_id)
  if (any(c("onset", "gap_models", "sensitivity", "runin", "runin_sensitivity") %in% stages) && isTRUE(cfg$model$adjust_device) && !"device_family" %in% names(persons))
    stop("the models use the device family as a covariate but no device information is in L1: run Cell 6c (export) and Cell 8b (run_etl_device) first, ",
         "or set cfg$model$adjust_device <- FALSE", call. = FALSE)
  if ("device_family" %in% names(persons)) cfg$device$reference_used <- device_reference(persons, cfg)
  streams <- store$streams
  model_streams <- model_streams %||% streams
  bad <- setdiff(model_streams, streams)
  if (length(bad)) stop("unknown stream(s) in model_streams: ", paste(bad, collapse = ", "), call. = FALSE)
  ckpt <- file.path(outdir, cfg$checkpoint$dirname)
  res <- list(); figs <- character(); units <- list(); early_tables <- character()
  has <- function(s) s %in% stages
  fig <- function(name, p, w = 7, h = 5) {
    if (!save_figures || is.null(p)) return(invisible(NULL))
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
      if (nrow(d_tab) > 0L) fig("forest_gap_definition", plot_forest(d_tab, "Odds of a Gap Starting Tomorrow per +1 SD Higher Steps, by Gap Definition", group_col = "stream", style = cfg$plot$interval_style %||% "se"), 10, 5)
    }
    for (s in model_streams) {                    # one forest plot of ALL predictors, and one observed-rate figure, per stream
      if (!is.null(res$onset_coefficients) && nrow(res$onset_coefficients) > 0L) {
        p <- plot_onset_forest_all(res$onset_coefficients, s, cfg)
        if (!is.null(p)) fig(paste0("onset_forest_other_", s), p)
      }
      if (!is.null(res$onset_deciles) && nrow(res$onset_deciles) > 0L) {
        p <- plot_onset_deciles(res$onset_deciles, s, cfg)
        if (!is.null(p)) fig(paste0("onset_deciles_", s), p, 8, 5)
      }
    }
    for (s in model_streams) {                    # model D (all terms by gap definition) and the exploratory grid, per stream
      oc <- res$onset_coefficients
      if (is.null(oc) || nrow(oc) == 0L) next
      dd <- oc[set_name == "joint" & adjusted & role %in% c("primary", "secondary_threshold")]
      glab <- setNames(sprintf("gap >= %d day%s", sort(unique(dd$k)), ifelse(sort(unique(dd$k)) > 1L, "s", "")), as.character(sort(unique(dd$k))))
      p <- if (nrow(dd) > 0L) plot_effects_by_group(dd, s, cfg, "k", glab, "other", paste0("Gap Definitions: All Terms of the Onset Model (", nice_label(s), " Stream)"), subtitle_or(cfg)) else NULL
      if (!is.null(p)) fig(paste0("gap_definition_all_terms_other_", s), p)
      p <- plot_exploratory_grid(oc, s, cfg); if (!is.null(p)) fig(paste0("exploratory_grid_effects_", s), p)
    }
    if (!is.null(res$risk_unadjusted) && "risk" %in% names(res$risk_unadjusted)) fig("risk_curve", plot_risk_curve(res$risk_unadjusted, cfg), 7, 9)
    if (nrow(res$cohort_specific) > 0L) fig("cohort_forest", plot_cohort_forest(res$cohort_specific, res$cohort_interaction, cfg), 10, 5)
  }

  # ---- 6. gap-level models (B and C): built from the gaps + one activity predictor, nothing else ----
  if (has("gap_models")) {
    act_col <- unique(c(cfg$predictors[["act_steps"]]$col, cfg$predictors[["act_mvpa"]]$col))   # steps (the exposure) and MVPA (the covariate of the *_mvpa units)
    for (s in model_streams) {
      gap_table <- function(standardization) lazy(function() {
        days <- store$get(s, cols = if (s == "activity") act_col else character())
        act  <- if (s == "activity") days else store$get("activity", cols = act_col)
        build_gap_level_data(days, act, persons, cfg, standardization = standardization)
      })
      gap_data <- gap_table("pre_gap"); gap_data_full <- gap_table("full_period")        # built only if a unit needs them
      units[[paste0("gap_models:", s)]] <- run_gap_stage(s, gap_data, cfg, ckpt, dkeys[[s]], gap_data_full)
      gap_data$release(); gap_data_full$release()
    }
    res <- c(res, collect_gap_results(ckpt, cfg, dkeys))
    if (nrow(res$model_B_curve) > 0L) fig("gap_dose_response", plot_gap_dose_response(res$model_B_curve, res$model_B_observed, cfg), 8, 9)
    if (nrow(res$model_C_predicted_probs) > 0L) fig("gap_length_probs", plot_ordinal_probs(res$model_C_predicted_probs, cfg), 10, 4.5)
    for (s in model_streams) {                    # one forest plot of ALL predictors per model and stream
      if (!is.null(res$model_B_all_terms) && nrow(res$model_B_all_terms) > 0L) {
        p <- plot_gapB_forest(res$model_B_all_terms, s, cfg); if (!is.null(p)) fig(paste0("gapB_forest_", s), p)
      }
      if (!is.null(res$model_C_all_terms) && nrow(res$model_C_all_terms) > 0L) {
        p <- plot_gapC_forest(res$model_C_all_terms, s, cfg); if (!is.null(p)) fig(paste0("gapC_forest_", s), p)
      }
      p <- plot_gap_mvpa(res$model_B_all_terms, s, cfg, "B"); if (!is.null(p)) fig(paste0("gapB_mvpa_effect_", s), p)
      p <- plot_gap_mvpa(res$model_C_all_terms, s, cfg, "C"); if (!is.null(p)) fig(paste0("gapC_mvpa_effect_", s), p)
      if (!is.null(res$model_B_observed_trajectory) && nrow(res$model_B_observed_trajectory) > 0L) {
        p <- plot_post_gap_observed(res$model_B_observed_trajectory, res$model_B_curve, s, cfg)
        if (!is.null(p)) fig(paste0("post_gap_activity_observed_", s), p, 10, 8)
      }
    }
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
      alt_full <- lazy(function() prepare_stream_bundle(s, inputs$get(), persons, cfg, standardization = "full_period"))   # person's usual from the whole period
      units[[paste0("sensitivity:", s)]] <- run_sensitivity_stage(s, bundle, alt, prim_specs[i], cfg, ckpt, dkeys[[s]], alt_full = alt_full)
      bundle$release(); alt_full$release(); if (!is.null(alt)) alt$release(); inputs$release()
    }
    res <- c(res, collect_sensitivity_results(ckpt, cfg, dkeys))
    if (nrow(res$sensitivity) > 0L) {
      ord <- c("main", "full_period_standardization", "trailing_baseline", "exclude_partial_wear", "first_gap_only", "main_sleep_only", "conditional_logit", "cluster_robust_logit")
      f_tab <- res$sensitivity[!is.na(or), .(stream, label = analysis, or, or_low, or_high)][order(match(label, ord))]
      if (nrow(f_tab) > 0L) fig("forest_sensitivity", plot_forest(f_tab, "Primary Effect Under Different Assumptions", group_col = "stream", style = cfg$plot$interval_style %||% "se"), 10, 5)
    }
    if (nrow(res$sensitivity_all_terms) > 0L) {
      for (s in model_streams) {                  # every within-person predictor under every analysis, per stream
        p <- plot_sensitivity_forest_stream(res$sensitivity_all_terms, s, cfg)
        if (!is.null(p)) fig(paste0("sensitivity_forest_", s), p)
        sl <- setNames(nice_label(intersect(SENSITIVITY_ORDER, unique(res$sensitivity_all_terms$analysis))), intersect(SENSITIVITY_ORDER, unique(res$sensitivity_all_terms$analysis)))
        p <- plot_effects_by_group(res$sensitivity_all_terms, s, cfg, "analysis", sl, "other",       # ALL terms, the primary analysis (main) against every sensitivity analysis
                                   paste0("Sensitivity Analyses: All Terms (", nice_label(s), " Stream)"), subtitle_or(cfg))
        if (!is.null(p)) fig(paste0("sensitivity_all_terms_other_", s), p)
      }
    }
  }

  # ---- device version (data quality): missing share and composition by cohort; person level, no model of the gap outcome ----
  if (has("device")) {
    dp <- store$device_person()
    if (is.null(dp)) stop("no device information in L1: run Cell 6c (export) and Cell 8b (run_etl_device) first", call. = FALSE)
    ds <- device_summary(dp, persons, cfg)
    res <- c(res, ds)
    fig("device_version_missing_and_composition", plot_device_composition(ds, cfg), 11, 6)
    fig("device_version_within_cohort", plot_device_versions(ds, cfg), 11, 8)
  }

  # ---- person-level % of days missing against every model variable (cohort, sex, watch version, age, between-person predictor levels) ----
  if (has("variables")) {
    dp <- store$device_person()
    if (is.null(dp)) log_info("variables: no device information in L1 (Cell 8b): the watch version is left out")
    ovs <- list()
    for (s in model_streams) {
      d <- store$get(s, cols = character()); pgs <- summarise_person_gaps(d, get_gap_runs(d)); rm(d)
      spec <- outcome_continuous_spec(s, cfg)
      inputs <- stream_inputs(store, s, cfg)
      means <- person_predictor_means(inputs, names(spec), cfg); rm(inputs)
      ovs[[s]] <- outcome_variable_analysis(outcome_variable_data(pgs, persons, dp, means, cfg), s, cfg, spec)
      rm(pgs, means); invisible(gc())
    }
    for (v in unique(unlist(lapply(ovs, function(a) a$tests$variable)))) fig(paste0("outcome_by_variable_", v), plot_outcome_variable(ovs, v, cfg))      # one figure per variable
    res$outcome_variable_tests      <- rbindlist(lapply(ovs, `[[`, "tests"), fill = TRUE)
    res$outcome_variable_contrasts  <- rbindlist(lapply(ovs, `[[`, "contrasts"), fill = TRUE)
    res$outcome_variable_models     <- rbindlist(lapply(ovs, `[[`, "models"), fill = TRUE)
    res$outcome_variable_summary    <- rbindlist(lapply(ovs, `[[`, "summary"), fill = TRUE)
    res$outcome_variable_cells      <- rbindlist(lapply(ovs, `[[`, "cells"), fill = TRUE)
    res$outcome_variable_by_cohort  <- rbindlist(lapply(ovs, `[[`, "by_cohort"), fill = TRUE)
  }

  # ---- 7b. run-in analysis (its own stages, tables and figures, apart from the main analysis) ----
  #   runin_description : no model.   runin : + the primary models (without a run-in term / with the run-in and its interaction), after the run-in.
  #   runin_sensitivity : + the other forms of the run-in (main effect only, categories, linear) with their own figures.
  if (has("runin_description") || has("runin") || has("runin_sensitivity")) {
    with_models <- has("runin") || has("runin_sensitivity"); with_sens <- has("runin_sensitivity")
    for (s in model_streams) {
      inputs <- lazy(function() stream_inputs(store, s, cfg))
      bundle <- lazy(function() prepare_stream_bundle(s, inputs$get(), persons, cfg))
      act_col <- cfg$predictors[["act_steps"]]$col
      gap_data <- if (with_models) lazy(function() {
        days <- store$get(s, cols = if (s == "activity") act_col else character())
        act  <- if (s == "activity") days else store$get("activity", cols = act_col)
        build_gap_level_data(days, act, persons, cfg)
      })
      units[[paste0("runin:", s)]] <- run_runin_stage(s, bundle, gap_data, grid = make_onset_grid(cfg), cfg, ckpt, dkeys[[s]], with_models = with_models, sensitivity = with_sens)
      bundle$release(); inputs$release(); if (!is.null(gap_data)) gap_data$release()
    }
    res <- c(res, collect_runin_results(ckpt, cfg, dkeys))
    res$runin_collinearity <- build_runin_collinearity(res)
    res$runin_term_tests   <- build_runin_term_tests(res)
    fig("runin_distribution", plot_runin_distribution(res$runin_distribution, cfg), 10, 7)
    fig("runin_rate_by_category", plot_runin_rate(res$runin_rate, res$runin_association, cfg), 11, 7)
    for (s in model_streams) fig(paste0("runin_steps_by_group_", s), plot_runin_steps_by_group(res$runin_steps_by_group, s, cfg), 12, 6)
    if (with_models) {
      # primary: the model without a run-in term and the model with the run-in and its interaction
      gl_prim <- stats::setNames(RUNIN_LABELS[RUNIN_PRIMARY_VARIANTS], RUNIN_PRIMARY_VARIANTS)
      fig("runin_primary_effect", plot_runin_primary_effect(res$runin_primary_effect, res$runin_tests, cfg), 11, 6)
      fig("runin_steps_by_missing_days", plot_runin_steps_by_missing(res$runin_effect_by_missing, res$runin_tests, cfg), 11, 7)
      forests <- function(prefix, groups_on, groups_gap, ttl) for (s in model_streams) {
        p <- if (nz(res$runin_onset_coefficients))
          plot_effects_by_group(res$runin_onset_coefficients, s, cfg, "runin_variant", groups_on, "other",
                                paste0("Onset Model After the Run-in", ttl, ": All Terms (", nice_label(s), " Stream)"), subtitle_runin(cfg)) else NULL
        if (!is.null(p)) fig(paste0(prefix, "onset_all_terms_other_", s), p)
        if (nz(res$runin_gapB_all_terms)) {
          p <- plot_effects_by_group(res$runin_gapB_all_terms, s, cfg, "runin_variant", groups_gap, "other",
                                     paste0("Model B After the Run-in", ttl, ": All Terms (", nice_label(s), " Stream)"),
                                     ov_wrap(paste0("Change in activity on the first full day back (SD of the person's usual): per doubling of the gap. ", covariate_note(cfg)), 100), measure = "estimate")
          if (!is.null(p)) fig(paste0(prefix, "gapB_all_terms_", s), p)
        }
        if (nz(res$runin_gapC_all_terms)) {
          p <- plot_effects_by_group(res$runin_gapC_all_terms, s, cfg, "runin_variant", groups_gap, "other",
                                     paste0("Model C After the Run-in", ttl, ": All Terms (", nice_label(s), " Stream)"),
                                     ov_wrap(paste0("Odds of a longer gap (ordinal model): steps before the gap per +1 SD. ", covariate_note(cfg)), 100))
          if (!is.null(p)) fig(paste0(prefix, "gapC_all_terms_", s), p)
        }
      }
      forests("runin_", gl_prim, gl_prim, "")
      # one GVIF figure per run-in model (the model with the run-in and its interaction)
      cfr <- if (nz(res$runin_collinearity)) res$runin_collinearity[figure == TRUE] else NULL
      if (nz(cfr)) for (s in model_streams) for (id in unique(cfr[stream == s, figure_id]))
        fig(paste0("runin_gvif_", id, "_", s), plot_collinearity(cfr[figure_id == id], s, cfg), 9, 2.6 + 0.24 * uniqueN(cfr[stream == s & figure_id == id, term]))
    }
    if (with_sens) {
      # sensitivity: every other form of the run-in, in figures of their own
      fig("runin_sensitivity_primary_effect", plot_runin_primary_effect(res$runin_primary_effect, res$runin_tests, cfg, variants = RUNIN_ONSET_VARIANTS, sensitivity = TRUE), 11, 7)
      fig("runin_sensitivity_shape", plot_runin_shape(res$runin_onset_coefficients, res$runin_onset_models, res$runin_distribution, res$runin_shape, cfg), 11, 6)
      gl_all <- stats::setNames(RUNIN_LABELS[RUNIN_ONSET_VARIANTS], RUNIN_ONSET_VARIANTS); gl_gap <- gl_all[RUNIN_GAP_VARIANTS]
      forests("runin_sensitivity_", gl_all, gl_gap, " (Sensitivity)")
    }
  }

  # ---- 7c. multicollinearity of every fitted model of the main analysis (GVIF, one figure per model) and the joint term tests ----
  if (any(c("onset", "gap_models", "sensitivity") %in% stages)) {
    if (is.null(res$onset_collinearity))       res <- c(res, collect_onset_results(ckpt, cfg, dkeys))
    if (is.null(res$sensitivity_collinearity)) res <- c(res, collect_sensitivity_results(ckpt, cfg, dkeys))
    if (is.null(res$gap_collinearity))         res <- c(res, collect_gap_results(ckpt, cfg, dkeys))
    res$model_collinearity <- build_collinearity_table(res)
    res$model_term_tests   <- build_term_test_table(res)
    cf <- collinearity_for_figure(res$model_collinearity)
    if (nz(cf)) for (s in model_streams) for (id in unique(cf[stream == s, figure_id])) {
      d <- cf[stream == s & figure_id == id]
      fig(paste0("collinearity_gvif_", id, "_", s), plot_collinearity(d, s, cfg), 9, 2.6 + 0.24 * uniqueN(d$term))
    }
  }

  # ---- 8. write everything, and say plainly what ran, what was resumed and what failed ----
  res$unit_run_status  <- unit_status_table(units)
  res$checkpoint_status <- checkpoint_report(ckpt, cfg, dkeys)
  res$unit_failures <- rbindlist(list(res$onset_failures, res$sensitivity_failures, res$gap_failures, res$runin_failures), fill = TRUE)
  if (nrow(res$unit_failures) > 0L && "stream" %in% names(res$unit_failures))      # only the streams of THIS run: an earlier failure of another stream is retried when that stream is run
    res$unit_failures <- res$unit_failures[stream %in% model_streams]
  if (!is.null(res$unit_failures) && nrow(res$unit_failures) > 0L) {
    warning(sprintf("%d model unit(s) FAILED - see unit_failures.tsv; run the same cell again to retry them.", nrow(res$unit_failures)), call. = FALSE)
    print(res$unit_failures)
  }
  write_tables(res[setdiff(names(res), early_tables)], outdir)
  write_table_index(outdir)                       # tables/table_index.tsv: each table with its analysis step
  if (length(figs)) write_figure_captions(figs, outdir, cfg)
  if (any(c("onset", "gap_models", "sensitivity", "runin", "runin_sensitivity") %in% stages)) write_results_summary(res, cfg, outdir)   # model results only
  if (!identical(stages, "missing_figures")) write_run_metadata(cfg, outdir, store, extra = list(stages = stages, onset_roles = onset_roles))
  invisible(res)
}
