# ==============================================================================
# report.R
# ------------------------------------------------------------------------------
# Writing results and the information needed to reproduce a run.
#   write_tables()          every data.table in a named list -> <outdir>/tables/<name>.tsv
#   write_run_metadata()    config, package versions, input fingerprints (md5 of L1), git commit
#   write_figure_captions() what each figure shows and how to read it
#   write_results_summary() the main numbers in plain language, from the saved tables only
#   unit_status_table()     which models ran, were resumed from a checkpoint, or failed in this call
# ==============================================================================

write_tables <- function(tables, outdir) {
  dir.create(file.path(outdir, "tables"), recursive = TRUE, showWarnings = FALSE)
  for (nm in names(tables)) {
    x <- tables[[nm]]
    if (is.data.frame(x) && nrow(x) > 0L) {
      fwrite(x, file.path(outdir, "tables", paste0(nm, ".tsv")), sep = "\t")
    }
  }
  invisible(names(tables))
}

# Input manifest: which L1 files (md5) or in-memory tables the run used, so a result can be tied to
# the exact input it came from. Works from the day store; nothing big is read for file-backed stores.
input_manifest <- function(store) {
  rbindlist(lapply(c(store$streams, "persons"), function(nm) {
    fp <- store$fingerprint(nm)
    data.table(input = nm, md5 = fp$md5 %||% NA_character_, bytes = fp$bytes %||% NA_real_,
               n_rows = fp$n_rows %||% NA_integer_, n_persons = fp$n_persons %||% NA_integer_)
  }))
}

git_commit <- function() {
  tryCatch(trimws(system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE)[1L]),
           error = function(e) NA_character_, warning = function(w) NA_character_)
}

write_run_metadata <- function(cfg, outdir, store = NULL, extra = list()) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    jsonlite::write_json(c(list(config = cfg), extra), file.path(outdir, "run_config.json"),
                         auto_unbox = TRUE, pretty = TRUE, null = "null")
  } else {
    saveRDS(c(list(config = cfg), extra), file.path(outdir, "run_config.rds"))
  }
  writeLines(c(paste("run time:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
               paste("git commit:", git_commit()), "", capture.output(sessionInfo())),
             file.path(outdir, "session_info.txt"))
  if (!is.null(store)) fwrite(input_manifest(store), file.path(outdir, "input_manifest.tsv"), sep = "\t")
  invisible(outdir)
}

# ---- execution summary: what is done, what failed, what was skipped ---------------------------------------

# The analysis step and the model a figure belongs to, in the order of the workflow:
#   1 descriptive statistics -> 2 temporal patterns -> 3 data quality / wear-consistency characterization -> 4 primary models
#   -> 5 stratified analyses -> 6 sensitivity analyses -> 7 model diagnostics
# (the exploratory specifications are listed with the primary models, step 4). Returns data.table(step, analysis, model).
figure_analysis <- function(name) {
  rules <- list(
    c("^device_",                              "3", "Data quality: device version",            "device version (descriptive, person level)"),
    c("^outcome_by_variable",                  "1", "Descriptive statistics: outcome by model variable", "descriptive (person level)"),
    c("^wear_consistency_",                    "3", "Data quality: wear consistency",          "wear coverage (descriptive)"),
    c("^wear_strata_",                         "5", "Stratified analyses: wear consistency",   "gap onset / model B / model C within strata"),
    c("^collinearity_",                        "7", "Model diagnostics",                       "all models (GVIF)"),
    c("^sensitivity_|^forest_sensitivity",     "6", "Sensitivity analyses",                    "gap onset (primary model, alternative assumptions)"),
    c("^onset_|^risk_curve|^cohort_forest|^forest_gap_definition|^gap_definition_|^exploratory_grid_", "4", "Primary models: gap onset", "gap onset (primary, model D, exploratory)"),
    c("^gapB_|^gap_dose_response|^post_gap_activity_observed_", "4", "Primary models: gap-level activity after a gap", "model B"),
    c("^gapC_|^gap_length_probs",              "4", "Primary models: gap-level gap length",    "model C"),
    c("^transition$|^missing_by_",             "2", "Temporal patterns of missingness",        "descriptive / temporal"),
    c("^heatmap_|^missing_pct_|^longest_gap|^gap_ecdf|^gap_length_distribution|^metrics_", "1", "Descriptive statistics", "descriptive"))
  for (r in rules) if (grepl(r[1L], name)) return(data.table(step = as.integer(r[2L]), analysis = r[3L], model = r[4L]))
  data.table(step = NA_integer_, analysis = NA_character_, model = NA_character_)
}

# The same for the result tables: which analysis step each table in <outdir>/tables belongs to.
table_analysis <- function(name) {
  rules <- list(
    c("collinearity|term_tests|^unit_|^checkpoint_|failures",  "7", "Model diagnostics and run status"),
    c("^device_",                                              "3", "Data quality: device version"),
    c("^outcome_variable",                                     "1", "Descriptive statistics: outcome by model variable"),
    c("^wear_consistency_",                                    "3", "Data quality: wear consistency"),
    c("^wear_strata_",                                         "5", "Stratified analyses: wear consistency"),
    c("^sensitivity",                                          "6", "Sensitivity analyses"),
    c("^onset_|^primary_test|^vif$|^cohort_interaction|^cohort_specific|^risk_unadjusted", "4", "Primary models: gap onset"),
    c("^model_B",                                              "4", "Primary models: gap-level activity after a gap (model B)"),
    c("^model_C",                                              "4", "Primary models: gap-level gap length (model C)"),
    c("^transition|^study_week|^calendar_year|^weekday|^missing_trend|^missing_weekday|^age_rates", "2", "Temporal patterns of missingness"),
    c("^cohort_|^missing_|^metric_|^selection_check|^gap_",     "1", "Descriptive statistics"))
  for (r in rules) if (grepl(r[1L], name)) return(data.table(step = as.integer(r[2L]), analysis = r[3L]))
  data.table(step = NA_integer_, analysis = NA_character_)
}

write_table_index <- function(outdir) {
  tdir <- file.path(outdir, "tables")
  nm <- sort(tools::file_path_sans_ext(list.files(tdir, pattern = "\\.tsv$")))
  nm <- setdiff(nm, "table_index")
  if (!length(nm)) return(invisible(NULL))
  idx <- rbindlist(lapply(nm, function(n) cbind(data.table(table = paste0(n, ".tsv")), table_analysis(n))))
  fwrite(idx[order(step, table)], file.path(tdir, "table_index.tsv"), sep = "\t")
  invisible(idx)
}

# One line per figure PRESENT in <outdir>/figures (not only the ones made in this call): file name + what it
# shows and how to read it. figure_index.tsv lists the same figures with the analysis step and model they belong to.
write_figure_captions <- function(names = NULL, outdir, cfg = NULL) {
  fdir <- file.path(outdir, "figures")
  dir.create(fdir, recursive = TRUE, showWarnings = FALSE)
  names <- sort(unique(c(names, tools::file_path_sans_ext(list.files(fdir, pattern = "\\.png$")))))
  cap <- vapply(names, caption_for, character(1), cfg = cfg)
  writeLines(c("# Figure captions", "", paste0("* **", names, ".png** - ", ifelse(is.na(cap), "(no caption)", cap))),
             file.path(fdir, "figure_captions.md"))
  idx <- rbindlist(lapply(names, function(n) cbind(data.table(figure = paste0(n, ".png")), figure_analysis(n))))
  idx[, caption := unname(cap)]
  fwrite(idx[order(step, figure)], file.path(fdir, "figure_index.tsv"), sep = "\t")
  invisible(cap)
}

# Which units ran, were resumed from a checkpoint, or failed in this call (a data.table).
unit_status_table <- function(status_list) {
  if (length(status_list) == 0L) return(data.table(stream = character(), unit = character(), status = character()))
  rbindlist(lapply(names(status_list), function(s) data.table(stream = s, unit = names(status_list[[s]]), status = unname(status_list[[s]]))))
}

# Short plain-language summary of the main findings, written from the saved result tables only.
write_results_summary <- function(res, cfg, outdir) {
  L <- c("# Results summary", "", "Written from the saved result tables (no model is refitted).", "")
  fmt_or <- function(or, lo, hi) sprintf("%.2f (%.2f-%.2f)", or, lo, hi)
  if (!is.null(res$primary_test) && nrow(res$primary_test) > 0L) {
    L <- c(L, "## Gap onset (primary model, gap >= 1 day)",
           "Odds that a gap starts tomorrow per +1 SD higher than the person's usual steps yesterday (OR < 1 = higher activity before fewer gaps; OR > 1 = higher activity before more gaps).",
           res$primary_test[, sprintf("* %s: OR %s, p=%s, Holm p=%s [%s, %s events]", stream, fmt_or(or, or_low, or_high),
                                      signif(p, 3), signif(p_holm, 3), status, format(n_events, big.mark = ","))], "")
  }
  if (!is.null(res$sensitivity) && nrow(res$sensitivity) > 0L) {
    r <- res$sensitivity[!is.na(or), .(lo = min(or), hi = max(or), n = .N), by = stream]
    L <- c(L, "## Sensitivity analyses", "Range of the primary-term odds ratio over the estimable sensitivity analyses:",
           r[, sprintf("* %s: %.2f to %.2f over %d analyses", stream, lo, hi, n)], "")
  }
  if (!is.null(res$sensitivity) && nrow(res$sensitivity) > 0L && "full_period_standardization" %in% res$sensitivity$analysis) {
    f <- merge(res$sensitivity[analysis == "main", .(stream, or_main = or)], res$sensitivity[analysis == "full_period_standardization", .(stream, or_full = or)], by = "stream")
    L <- c(L, "## Standardization: pre-gap (primary) vs full observation period (sensitivity)",
           "The person's usual is estimated only from days observed before each day (primary) or from the whole observation period (look-ahead). Odds ratio per +1 SD higher steps:",
           f[, sprintf("* %s: primary %.2f, full period %.2f", stream, or_main, or_full)], "")
  }
  if (!is.null(res$model_collinearity) && nrow(res$model_collinearity) > 0L) {
    fl <- res$model_collinearity[figure == TRUE & flag != "ok"]
    L <- c(L, "## Multicollinearity (VIF-equivalent, models drawn in collinearity_gvif_<stream>)",
           if (nrow(fl) == 0L) "No term reached a VIF-equivalent of 5." else c("Terms with a VIF-equivalent >= 5 (high) or >= 10 (severe); nothing was removed because of it. They widen that term's interval and do not bias the others:",
             fl[, .(max_vif = max(vif_equiv)), by = .(stream, term, flag)][order(stream, -max_vif)][, sprintf("* %s: %s (%s, max %.1f)", stream, term, flag, max_vif)]), "")
  }
  if (!is.null(res$wear_strata_comparison) && nrow(res$wear_strata_comparison) > 0L) {
    L <- c(L, "## Wear consistency (consistent vs inconsistent wear; reference: consistent)",
           "Primary term by stratum and the test that the strata differ (not the significance of either stratum):",
           res$wear_strata_comparison[, sprintf("* %s, %s: %.2f (%.2f-%.2f) vs %.2f (%.2f-%.2f); p(difference) = %s", model, stream, effect_consistent, low_consistent, high_consistent,
                                                effect_inconsistent, low_inconsistent, high_inconsistent, signif(p_difference, 2))], "")
    if (!is.null(res$wear_strata_quality) && nrow(res$wear_strata_quality[usable == FALSE]) > 0L)
      L <- c(L, "Strata that are sparse (< events per coefficient threshold) or did not converge cleanly are listed in wear_strata_quality.tsv.", "")
  }
  B <- res$model_B_post_gap_activity
  if (!is.null(B) && nrow(B) > 0L) {
    L <- c(L, "## Model B (exploratory): activity on the first full day back, per doubling of gap length (SD units)",
           B[model_id == "B_back2_adj", sprintf("* %s: %.3f (%.3f to %.3f), p=%s [%s]", stream, estimate, conf_low, conf_high,
                                                signif(p, 3), status)], "")
  }
  C <- res$model_C_gap_length
  if (!is.null(C) && nrow(C) > 0L) {
    L <- c(L, "## Model C (exploratory): odds of a longer gap per +1 SD higher than usual steps (OR > 1 = longer gaps after more active days)",
           C[model_id == "C_adj", sprintf("* %s: %s [%s]", stream,
                                          ifelse(is.na(or), paste0("not estimated (", note, "); proportional-odds p=", signif(po_p, 2)), sprintf("OR %s", fmt_or(or, or_low, or_high))), status)], "")
  }
  M <- res$gap_mvpa_effects
  if (!is.null(M) && nrow(M) > 0L && any(M$term == "w_pre_day_mvpa")) {
    m <- M[term == "w_pre_day_mvpa"]
    L <- c(L, "## MVPA before the gap as a covariate of models B and C (exploratory)",
           "MVPA minutes (log) on the last observed day before the gap, against the person's usual before that day, per +1 SD higher. Model B: change in activity on the first full day back (SD of the person's usual steps); model C: odds ratio of a longer gap.",
           m[, sprintf("* model %s %s, %s: %s (95%% CI %s to %s), p=%s [%s]", model, model_id, stream,
                       ifelse(model == "C", sprintf("OR %.2f", or), sprintf("%.3f", estimate)),
                       ifelse(model == "C", sprintf("%.2f", exp(conf_low)), sprintf("%.3f", conf_low)),
                       ifelse(model == "C", sprintf("%.2f", exp(conf_high)), sprintf("%.3f", conf_high)), signif(p, 3), status)], "")
  }
  fl <- rbindlist(list(res$onset_failures, res$sensitivity_failures, res$gap_failures), fill = TRUE)
  if (nrow(fl) > 0L)
    L <- c(L, "## FAILED units (re-run the same cell to retry them)", fl[, sprintf("* %s / %s / %s: %s", kind, stream, id, error)], "")
  writeLines(L, file.path(outdir, "results_summary.md"))
  invisible(L)
}
