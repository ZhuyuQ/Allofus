# ==============================================================================
# notebook_cells.R - the notebook cells, in order (each "# ==== CELL n" block = one notebook cell)
# These were rehearsed top to bottom on a simulated workbench (local exports in the bucket layout,
# a stand-in gsutil, fake final_pid / CCI files); only Cell 6 (BigQuery) could not be run for real.
# Edit BASE (Cells 1 and 2) if your workspace path differs.
#
# AFTER A KERNEL CRASH OR RESTART:  run Cell 2 (setup), Cell 10 (open the saved day tables), then the
# cell you were in. Every finished model is a checkpoint in OUT_DIR/checkpoints and is skipped; only
# missing, failed or outdated models run again. You never need to repeat Cells 5-9 (ETL) unless L1 is gone.
#
# RESOURCES (cfg$resources, defaults for 4 CPU cores / 32 GB RAM): one heavy model at a time in its own
# worker process, at most 24 GB used by a model and 4 GB always left free. Change them in Cell 2.
# ==============================================================================
# ==== CELL 1
BASE <- "/home/jupyter/workspace/rw-migration-aou-rw-6f1dba72/workspaces"
unzip(file.path(BASE, "temporaldigitalhealthphenotyping/missingness.zip"),
      exdir = file.path(BASE, "temporaldigitalhealthphenotyping"))
print(list.files(file.path(BASE, "temporaldigitalhealthphenotyping/missingness")))

# ==== CELL 2
BASE <- "/home/jupyter/workspace/rw-migration-aou-rw-6f1dba72/workspaces"
MISS <- file.path(BASE, "temporaldigitalhealthphenotyping/missingness")
source(file.path(MISS, "load_missingness.R")); load_missingness(MISS)

EXPORT_DATE <- "20261006"
FINAL_PID   <- file.path(BASE, "result/final/final_pid.txt")
CCI_FILE    <- file.path(BASE, "result/pre_final/pre_final_cci.txt")
L1_DIR      <- file.path(BASE, "result/final/missingness/L1")
OUT_DIR     <- file.path(BASE, "result/final/missingness/results")
OVERLAP_POLICY <- "exclude"

cfg <- make_config(list(export  = list(date = EXPORT_DATE),
                        overlap = list(policy = OVERLAP_POLICY)))
# Only needed for Cell 6 (BigQuery export). If the CDR below is blank or not "project.dataset"
# (e.g. after the workspace migration), uncomment and fill in these two lines:
# cfg$export$cdr     <- "<project>.<dataset>"
# cfg$export$billing <- "<your-google-project>"
# Machine budget (the defaults suit 4 CPU cores / 32 GB RAM; change only if your runtime differs):
# cfg$resources$memory_budget_gb  <- 24   # most memory one model fit may use
# cfg$resources$memory_reserve_gb <- 4    # always keep this much free
# cfg$resources$light_workers      <- 3    # models B / C of a stream run side by side (each < 1 GB; capped by free memory)
cat("bucket :", cfg$export$bucket, "\nowner  :", cfg$export$owner, "\n")
cat("CDR    :", cfg$export$cdr, "\nbilling:", cfg$export$billing, "\n")
stopifnot("final_pid.txt not found - fix FINAL_PID"   = file.exists(FINAL_PID),
          "CCI file not found - fix CCI_FILE"         = file.exists(CCI_FILE),
          "WORKSPACE_BUCKET / OWNER_EMAIL not set"    = nzchar(cfg$export$bucket) && nzchar(cfg$export$owner))

# ==== CELL 3
need <- c("data.table", "lme4", "ordinal", "survival", "sandwich", "ggplot2", "rstatix", "ggpubr")
ok <- sapply(need, requireNamespace, quietly = TRUE); print(ok)
stopifnot("install the packages that are FALSE (install.packages(\"name\")), then run this cell again" = all(ok))
if (!requireNamespace("bigrquery", quietly = TRUE)) message("note: bigrquery is not installed; it is only needed for Cell 6 (BigQuery export)")

# ==== CELL 4
cfg_demo <- make_config(list(boot = list(B = 50L)))
inp_demo <- simulate_inputs(cfg_demo)
res_demo <- run_missingness_pipeline(inp_demo$day_tables, inp_demo$persons, cfg_demo,
                                     file.path(tempdir(), "demo"),
                                     stages = c("missingness", "temporal", "descriptive", "onset"),
                                     onset_roles = "primary")
print(res_demo$primary_test[, .(stream, or = round(or, 2), p_holm)])

# ==== CELL 5
final_pid <- read_final_pid(FINAL_PID)
stopifnot("cohort names in final_pid must be Insomnia / Diabetes / Healthy / Prostate" =
            all(final_pid$cohort %in% cfg$cohorts))
print(final_pid[, .N, by = cohort])
ov <- resolve_cohort_overlap(final_pid, OVERLAP_POLICY, cfg$overlap$priority)
cat("unique people:", ov$n_unique_people, "| in more than one cohort:", length(ov$overlap_ids),
    "| kept with policy", OVERLAP_POLICY, ":", ov$n_kept, "\n")
print(ov$overlap)

# ==== CELL 6a
exports <- run_bq_exports(cfg, final_pid, cohorts = "Prostate")

# ==== CELL 6b
exports <- run_bq_exports(cfg, final_pid)

# ==== CELL 7
for (co in cfg$cohorts) for (s in c("activity", "sleep", "hr", "meta")) {
  n <- tryCatch(length(gsutil_ls(export_path(cfg, co, s))), error = function(e) 0)
  cat(sprintf("%-9s %-9s %d file(s)\n", co, s, n))
}

# ==== CELL 8
etl <- run_etl(cfg, L1_DIR, final_pid_path = FINAL_PID, cci_path = CCI_FILE)   # one stream in memory at a time
print(etl$summary)                                   # rows / people / date range of each day table
print(etl$dropped_columns)                           # columns removed because they are mostly NA or constant
print(etl$persons[, .N, by = cohort])

# ==== CELL 9
store <- open_day_store(cfg, L1_DIR)
hr <- store$get("hr")
hr_cols <- intersect(c("total_hr_min", "elevated_min", "mins_out_of_range", "mins_fat_burn",
                       "mins_cardio", "mins_peak", "total_cal"), names(hr))
print(summary(hr[, ..hr_cols]))
rm(etl, store, hr); invisible(gc())                  # the day tables are saved on disk (L1); Cell 10 opens them

# ==== CELL 10
inp <- load_inputs(cfg, L1_DIR)                      # nothing big is read: each stream is read when a stage needs it
print(data.frame(file = basename(inp$day_tables$paths), MB = round(file.size(inp$day_tables$paths) / 1e6, 1)))
cat("people:", nrow(inp$persons), "\n")

# ==== CELL 10b   (optional: memory right now)
cat(sprintf("this R session: %.1f GB | free: %.1f GB | a model may use up to: %.1f GB\n",
            rss_gb(), available_memory_gb(), memory_headroom_gb(cfg)))

# ==== CELL 11
res <- run_missingness_pipeline(inp$day_tables, inp$persons, cfg, OUT_DIR,
                                stages = c("missingness", "temporal", "descriptive"))
print(res$cohort_missingness)
print(res$selection_check)

# ==== CELL 11b   (optional, ~1-2 min: redraws ONLY the missingness-summary figures; no models, nothing else is rerun)
# missing % by sex / age group / race / CCI group, the gap ECDF, and missing share by study week, year, weekday and age;
# the statistics behind them are saved in results/tables (missing_pct_by_stratum_tests / _models, missing_trend_tests,
# missing_trend_heterogeneity, missing_weekday_tests, missing_weekday_odds_ratios)
res <- run_missingness_pipeline(inp$day_tables, inp$persons, cfg, OUT_DIR, stages = "missing_figures")
print(res$missing_pct_by_stratum_tests[family == "groups_within_cohort" & test != "none",
      .(stream, stratum_var, cohort = group, test, n_total, statistic = round(statistic, 2), df, p = signif(p, 3),
        p_holm = signif(p_holm, 3), effect = round(effect, 3), effect_name)])
print(res$missing_trend_tests[, .(stream, variable, cohort, unit, or = round(odds_ratio, 3), or_low = round(or_low, 3), or_high = round(or_high, 3),
      p = signif(p, 3), p_holm = signif(p_holm, 3), status)])

# ==== CELL 12   (Cells 12-14: onset models; new figures per stream: onset_forest_year_<stream> + onset_forest_other_<stream> = all predictors in two plots, onset_deciles_<stream> = observed gap rate by decile.
#                 Finished models are skipped; only the small "deciles" unit is computed on the first run after the update)
res <- run_missingness_pipeline(inp$day_tables, inp$persons, cfg, OUT_DIR, stages = "onset",
                                onset_roles = c("primary", "secondary_threshold"),
                                model_streams = "activity")
print(res$primary_test)

# ==== CELL 13
res <- run_missingness_pipeline(inp$day_tables, inp$persons, cfg, OUT_DIR, stages = "onset",
                                onset_roles = c("primary", "secondary_threshold"),
                                model_streams = "sleep")
print(res$primary_test)

# ==== CELL 14
res <- run_missingness_pipeline(inp$day_tables, inp$persons, cfg, OUT_DIR, stages = "onset",
                                onset_roles = c("primary", "secondary_threshold"),
                                model_streams = "hr")
print(res$primary_test)

# ==== CELL 15
# New figures: gapB_forest_<stream> and gapC_forest_<stream> (all predictors of each model). They need every coefficient, so models
# B and C saved by an older version of the code are refitted ONCE (checkpoint schema 2); model C is the slow one.
# Models B and C: 6 independent units per stream (4 x B, 2 x C), each saved on its own. They are built from
# the gaps and one activity predictor only (no person-day table), so they need little memory.
res <- run_missingness_pipeline(inp$day_tables, inp$persons, cfg, OUT_DIR, stages = "gap_models",
                                model_streams = c("activity", "sleep", "hr"))
print(res$model_B_post_gap_activity[, .(stream, model_id, estimate = round(estimate, 3), p = signif(p, 3), status)])
print(res$model_C_gap_length[, .(stream, model_id, or = round(or, 2), po_p = signif(po_p, 2), status)])
cat("failed units:", nrow(res$unit_failures), "\n")      # 0 = nothing failed; otherwise run this cell again to retry them
if (nrow(res$unit_failures) > 0L) print(res$unit_failures)

# ==== CELL 16   (new figure per stream: sensitivity_forest_<stream> = every within-person predictor under every analysis; sensitivity units
#                saved by an older version of the code are refitted ONCE (checkpoint schema 2) - this is the long cell)
res <- run_missingness_pipeline(inp$day_tables, inp$persons, cfg, OUT_DIR, stages = "sensitivity",
                                model_streams = c("activity", "sleep", "hr"))
print(res$sensitivity[, .(stream, analysis, term, or = round(or, 2), p = signif(p, 3), status)])

# ==== CELL 17
cfg_fast <- make_config(list(export  = list(date = EXPORT_DATE),
                             overlap = list(policy = OVERLAP_POLICY),
                             model   = list(exploratory_sample_persons = 3000L)))
res <- run_missingness_pipeline(inp$day_tables, inp$persons, cfg_fast, OUT_DIR, stages = "onset",
                                onset_roles = "exploratory", model_streams = c("activity", "sleep", "hr"))
print(table(res$onset_models$role, res$onset_models$stream))

# ==== CELL 18
CKPT <- file.path(OUT_DIR, cfg$checkpoint$dirname)
print(checkpoint_report(CKPT, cfg)[, .N, by = .(kind, state)])     # valid / stale / failed / corrupt checkpoints
cat(readLines(file.path(OUT_DIR, "results_summary.md")), sep = "\n")
print(list.files(file.path(OUT_DIR, "tables")))
print(list.files(file.path(OUT_DIR, "figures")))                   # figure_captions.md explains every figure
system2("gsutil", c("-m", "cp", "-r", OUT_DIR, file.path(Sys.getenv("WORKSPACE_BUCKET"), "missingness_results/")))
