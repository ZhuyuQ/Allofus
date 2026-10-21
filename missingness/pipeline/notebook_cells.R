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
# ORDER OF THE ANALYSES (the reading order of README.md / STATISTICS.md):
#   descriptive statistics (Cell 11, 11b) -> data quality (device, Cell 11d) and the run-in description (Cell 11c) -> primary models (Cells 12-15) ->
#   sensitivity analyses (Cell 16) -> exploratory specifications (Cell 17) -> the run-in analysis (Cell 18; sensitivity Cell 18b) ->
#   export (Cell 19). Multicollinearity diagnostics (GVIF) are written by every model cell; Cell 18 prints their summary.
# STANDARDIZATION: every within- and between-person z-score uses only data observed BEFORE each day / gap (cfg$standardization);
#   Cell 16 adds the full-observation-period standardization as a sensitivity analysis.
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

# ==== CELL 6c   (OPTIONAL, separate from the exports above: the Fitbit device version of each person; nothing else is changed)
# 1) is there a device table with a device_version column?  2) does every person have a device version?  3) export it (4 files, one per cohort).
probe <- bq_probe_device_table(cfg)
print(probe$columns_check); print(probe$sample)
stopifnot("the device table has no device_version column (see probe$table_columns)" = probe$all_found)
dev_cov <- bq_device_coverage(cfg, final_pid)
print(bq_device_versions(cfg, final_pid), nrows = 100)   # which devices exist (device_type SCALE = a Fitbit scale, not a wearable)
print(dev_cov)      # pct_activity_with_version = share of people with activity data who have a device version; pct_days_with_version = share of their activity days
dev_exports <- run_bq_exports(cfg, final_pid, streams = "device")

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

# ==== CELL 8b   (OPTIONAL, ~1 min: device version -> L1. Reads the 4 device exports of Cell 6c; the day tables and persons.rds are NOT touched,
#                no model result is invalidated. Writes device_person.rds, device_timeline.rds, device_status.tsv into L1_DIR.)
dev <- run_etl_device(cfg, L1_DIR)
print(dev$summary)                                   # people per cohort and status (wearable version known / scale only / no device record ...)

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

# ==== CELL 11c   (run-in: description of the wear in the FIRST 28 days, no model; builds each stream's tables once, a few minutes per stream)
# Run-in = the first cfg$runin$days (28) days of each person's follow-up in the stream. Missing days in the run-in is the wear measure of the run-in
# analysis; everything modelled later (Cell 18) starts AFTER the run-in, so the measure always comes first (whole-period coverage is never used:
# it contains the gaps that are the outcome). Hypotheses: H1 = more missing days in the run-in -> more gaps later; H2 = the steps effect depends on
# the run-in. Figures: runin_distribution, runin_rate_by_category (H1), runin_steps_by_group_<stream> (H2); tables: runin_flow, runin_distribution,
# runin_rate, runin_steps_by_group, runin_association.
res <- run_missingness_pipeline(inp$day_tables, inp$persons, cfg, OUT_DIR, stages = "runin_description",
                                model_streams = c("activity", "sleep", "hr"))
print(res$runin_flow[cohort == "All cohorts", .(stream, n_people, n_complete_runin, n_analysed, n_events_after_runin)])      # who is in the run-in analysis
print(res$runin_distribution[cohort == "All cohorts", .(stream, category, n_people, pct = round(pct, 1))])
print(res$runin_association[, .(stream, cohort, n_people, rho = round(rho, 2), p_holm = signif(p_holm, 2))])                 # person level: run-in missing days vs later gap rate

# ==== CELL 11d   (data quality: missing Fitbit device version and device composition by cohort; person level, no model; needs Cells 6c and 8b; seconds)
# WATCH data only: scales, clip-on trackers (One, Zip, Ultra, Classic) and phone tracking (MobileTrack) are NOT a watch version (cfg$device$forms / wrist_forms; the
# default counts smartwatches and wristbands; use cfg$device$wrist_forms <- "smartwatch" for smartwatches only, then rerun Cells 8b and 11d).
# Missing = no device record, scale only (weight only) or no version; "Clip-on or phone" = activity data but no watch (its own category, not missing).
# Figures device_version_missing_and_composition (one bar per cohort: grey = missing, brown = clip-on or phone, colours = watch family) and device_version_within_cohort
# (the device versions within each cohort, with Wilson intervals and a test per version), and the tables device_*:
# missing share with Wilson intervals, chi-square / Cramer's V / Fisher tests, the model adjusted for age and sex, residuals, full version list.
res <- run_missingness_pipeline(inp$day_tables, inp$persons, cfg, OUT_DIR, stages = "device")
print(res$device_missing_by_cohort)
print(res$device_tests[, .(test, chisq = round(chisq, 1), df, p_asymptotic = signif(p_asymptotic, 3), p_monte_carlo = signif(p_monte_carlo, 3), cramers_v = round(cramers_v, 3))])
print(res$device_missing_vs_reference)
print(res$device_version_tests[, .(device_version, group, n_people, pct_min = round(pct_min, 1), pct_max = round(pct_max, 1), cramers_v = round(cramers_v, 3), p_holm = signif(p_holm, 3))])   # per version: does its share differ between cohorts?
show_figure(OUT_DIR, "device_version_missing_and_composition")   # the figures are saved in <OUT_DIR>/figures; these lines also show them here
show_figure(OUT_DIR, "device_version_within_cohort")

# ==== CELL 11e   (person-level % of days missing against EVERY model variable, ONE figure per variable; ~3-5 min, no model of the gap outcome)
# Variables: cohort, sex, CCI group, watch version, age, CCI score, and the between-person level of each model predictor. Each figure has the streams side by side: a continuous variable is one line per cohort
# with a band over 10 bins, a categorical one is one violin per cohort within each level; the statistics are written below the figure. Figures outcome_by_variable_<variable>;
# tables outcome_variable_tests (overall test, effect size with bootstrap interval, Holm over the variables of the stream, cohort-difference test), _by_cohort (rank correlation within each cohort),
# _contrasts (each level against the reference), _models (odds ratio of a missing day, unadjusted and adjusted), _summary, _cells (the numbers behind the figures).
# The watch version needs Cells 6c and 8b; without them it is left out.
res <- run_missingness_pipeline(inp$day_tables, inp$persons, cfg, OUT_DIR, stages = "variables")
print(res$outcome_variable_tests[, .(stream, variable, test, n, p_holm = signif(p_holm, 3), effect_name, effect = round(effect, 3), low = round(effect_low, 3), high = round(effect_high, 3),
                                     p_cohorts_differ_holm = signif(p_cohort_difference_holm, 3))])
for (v_ in unique(res$outcome_variable_tests$variable)) show_figure(OUT_DIR, paste0("outcome_by_variable_", v_))     # shows each saved figure here

# ==== CELL 12   (Cells 12-14: onset models. Figures per stream: onset_forest_other_<stream> = all predictors (incl. device family and CCI score; no calendar year),
#                 onset_deciles_<stream> = observed gap rate by decile, gap_definition_all_terms_*_<stream> (model D), exploratory_grid_effects_<stream>,
#                 collinearity_gvif_<model>_<stream> (one figure per model: onset_primary, onset_full_period, model_B, model_C). Predictors are standardized with data observed BEFORE each day (pre-gap); models saved by an older
#                 version of the code are refitted once.)
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
# Figures: gap_dose_response, gapB_forest_<stream>, gapC_forest_<stream> (all terms), post_gap_activity_observed_<stream> (observed steps after gaps
# with the Model B estimate drawn apart). Days after a gap are standardized against the person's usual BEFORE the gap. Models saved by an older
# version of the code are refitted ONCE; model C is the slow one. B_back2_adj_fullperiod / C_adj_fullperiod are the full-period sensitivity fits.
# MVPA before the gap as an extra covariate (units B_back2_adj_mvpa, B_back2_adj_pre_mvpa, C_adj_mvpa; the existing units are NOT refitted): table
# gap_mvpa_effects, figures gapB_mvpa_effect_<stream>, gapC_mvpa_effect_<stream>.
# Models B and C: independent units per stream (the B and C specifications and their variants), each saved on its own. They are built from
# the gaps and one activity predictor only (no person-day table), so they need little memory.
res <- run_missingness_pipeline(inp$day_tables, inp$persons, cfg, OUT_DIR, stages = "gap_models",
                                model_streams = c("activity", "sleep", "hr"))
print(res$model_B_post_gap_activity[, .(stream, model_id, estimate = round(estimate, 3), p = signif(p, 3), status)])
print(res$model_C_gap_length[, .(stream, model_id, or = round(or, 2), po_p = signif(po_p, 2), status)])
# MVPA before the gap (per +1 SD higher than the person's usual): model B = SD change in activity on the first full day back, model C = odds ratio of a longer gap
print(res$gap_mvpa_effects[term == "w_pre_day_mvpa", .(model, stream, model_id, estimate = round(estimate, 3), low = round(conf_low, 3), high = round(conf_high, 3),
                                                       or = round(or, 2), or_low = round(or_low, 2), or_high = round(or_high, 2), n_gaps, vif_equiv = round(vif_equiv, 2))])
cat("failed units:", nrow(res$unit_failures), "\n")      # 0 = nothing failed; otherwise run this cell again to retry them
if (nrow(res$unit_failures) > 0L) print(res$unit_failures)

# ==== CELL 16   (sensitivity analyses, incl. the NEW full-period standardization; figures per stream: sensitivity_forest_<stream> = within-person
#                predictors under every analysis, sensitivity_all_terms_other_<stream> = ALL terms (incl. the CCI refits cci_removed / cci_as_groups / cci_within_cohort); units saved by an older version are
#                refitted ONCE - this is the long cell)
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

# ==== CELL 18   (run-in analysis, PRIMARY: the onset model, model B and model C on the days / gaps AFTER the first 28 days, without a run-in term and with the run-in AND its interaction; long)
# The description of Cell 11c runs first inside this call. Same rows in both models. Run-in = log2(1 + missing days in the first 28 days), centred, per doubling;
# its interaction is with the steps term (onset, C) / the gap length (B). H1 = run-in effect, H2 = interaction (the confirmatory test, Holm over the streams),
# H3 = does the primary steps effect change. Figures (all prefixed runin_, apart from the main analysis): runin_primary_effect, runin_steps_by_missing_days,
# runin_onset_all_terms_other_<stream>, runin_gapB_all_terms_<stream>, runin_gapC_all_terms_<stream>, and ONE GVIF figure per model: runin_gvif_onset_<stream>,
# runin_gvif_model_B_<stream>, runin_gvif_model_C_<stream>. Other forms of the run-in are in Cell 18b (sensitivity).
res <- run_missingness_pipeline(inp$day_tables, inp$persons, cfg, OUT_DIR, stages = "runin",
                                model_streams = c("activity", "sleep", "hr"))
print(res$runin_tests[analysis == "primary", .(model, stream, hypothesis, term, effect = round(effect, 3), low = round(low, 3), high = round(high, 3), p = signif(p, 3), p_holm = signif(p_holm, 3))])
print(res$runin_primary_effect[variant %in% RUNIN_PRIMARY_VARIANTS, .(model, stream, variant, effect = round(effect, 3))])   # H3: the steps effect without / with the run-in terms
print(res$runin_quality[variant %in% RUNIN_PRIMARY_VARIANTS, .(model, stream, variant, n_people, n_events, events_per_coefficient = round(events_per_coefficient, 1), status, sparse, usable)])
print(res$runin_collinearity[variant == "interaction" & flag != "ok", .(stream, analysis, term, vif_equiv = round(vif_equiv, 1), flag)])   # run-in terms with a high / severe VIF-equivalent
cat("failed units:", nrow(res$unit_failures), "\n")
for (f_ in list.files(file.path(OUT_DIR, "figures"), pattern = "^runin_(primary|steps_by_missing|gvif)")) show_figure(OUT_DIR, sub("\\.png$", "", f_))

# ==== CELL 18b  (run-in analysis, SENSITIVITY: the other forms of the run-in - main effect only, run-in as categories (shape check), run-in linear; long)
# Runs the primary units first (finished ones are skipped). Figures (prefix runin_sensitivity_): runin_sensitivity_primary_effect, runin_sensitivity_shape,
# runin_sensitivity_onset_all_terms_other_<stream>, runin_sensitivity_gapB_all_terms_<stream>, runin_sensitivity_gapC_all_terms_<stream>.
# Tests (runin_tests, analysis = "sensitivity"): the likelihood-ratio test of H2 against the main-effect model, H1 / H2 with the linear coding, the shape of the coding.
res <- run_missingness_pipeline(inp$day_tables, inp$persons, cfg, OUT_DIR, stages = "runin_sensitivity",
                                model_streams = c("activity", "sleep", "hr"))
print(res$runin_tests[analysis == "sensitivity", .(model, stream, hypothesis, form, effect = round(effect, 3), p = signif(p, 3), p_holm = signif(p_holm, 3))])
print(res$runin_shape[, .(stream, variant, n_rows, aic = round(aic, 1), delta_aic_vs_main = round(delta_aic_vs_main, 1))])      # coding of the run-in: AIC of every form
cat("failed units:", nrow(res$unit_failures), "\n")

# ==== CELL 19
CKPT <- file.path(OUT_DIR, cfg$checkpoint$dirname)
print(checkpoint_report(CKPT, cfg)[, .N, by = .(kind, state)])     # valid / stale / failed / corrupt checkpoints
cat(readLines(file.path(OUT_DIR, "results_summary.md")), sep = "\n")
print(list.files(file.path(OUT_DIR, "tables")))
print(list.files(file.path(OUT_DIR, "figures")))                   # figure_captions.md explains every figure
system2("gsutil", c("-m", "cp", "-r", OUT_DIR, file.path(Sys.getenv("WORKSPACE_BUCKET"), "missingness_results/")))
