# ==============================================================================
# run_exports.R   -  STEP 0: BigQuery -> CSV files in the workspace bucket (16 exports)
# ------------------------------------------------------------------------------
# Needed only when the exports do not exist (e.g. after the workspace migration).
# Runs the same SQL as the earlier notebooks for every cohort x stream
# (activity, sleep, hr, meta). Finished exports are skipped, so after a failure
# just run it again. Use the SAME EXPORT_DATE in run_etl.R.
# ==============================================================================
BASE            <- "/home/jupyter/workspace/rw-migration-aou-rw-6f1dba72/workspaces"
MISSINGNESS_DIR <- file.path(BASE, "temporaldigitalhealthphenotyping/missingness")       # <- where you put the folder
FINAL_PID       <- file.path(BASE, "result/final/final_pid.txt")
EXPORT_DATE     <- "20261006"      # folder name for ALL 16 exports - use the same value in run_etl.R

source(file.path(MISSINGNESS_DIR, "load_missingness.R")); load_missingness(MISSINGNESS_DIR)
cfg <- make_config(list(export = list(date = EXPORT_DATE)))     # bucket / owner from WORKSPACE_BUCKET / OWNER_EMAIL
cat("CDR:", cfg$export$cdr, "| billing project:", cfg$export$billing, "\n")   # blank / not project.dataset? set cfg$export$cdr and cfg$export$billing
cat("exporting to:", file.path(cfg$export$bucket, "bq_exports", cfg$export$owner, cfg$export$date), "\n")

final_pid <- read_final_pid(FINAL_PID)
print(final_pid[, .N, by = cohort])                    # should be 1383 / 1298 / 250 / 6381 (9,312 in total)
exports <- run_bq_exports(cfg, final_pid)

# If an export fails with "exceeds 1 GB / use a wildcard" for sleep or hr, make that stream
# sharded for BOTH the export and the ETL:  make_config(list(..., streams = list(hr = list(sharded = TRUE))))
