# ==============================================================================
# run_etl.R   -  RUN ONCE per export date (All of Us workbench, R kernel)
# ------------------------------------------------------------------------------
# Raw BigQuery CSV exports (L0)  ->  compact day tables + persons table (L1).
#   * reads each cohort x stream with ONE reader (no more 16 copy-pasted readers)
#   * activity / sleep / hr -> one row per person-day; meta + cohort + CCI -> persons
#   * writes day_activity / day_sleep / day_hr (+ persons.rds) to L1_DIR
# Re-run only when the raw export changes. All later analyses read L1_DIR.
# ==============================================================================
BASE            <- "/home/jupyter/workspace/rw-migration-aou-rw-6f1dba72/workspaces"
MISSINGNESS_DIR <- file.path(BASE, "temporaldigitalhealthphenotyping/missingness")       # <- where you put the folder
L1_DIR          <- file.path(BASE, "result/final/missingness/L1")
FINAL_PID       <- file.path(BASE, "result/final/final_pid.txt")
CCI_FILE        <- file.path(BASE, "result/pre_final/pre_final_cci.txt")
EXPORT_DATE     <- "20261006"      # folder date used by run_exports.R (the date the 16 exports were written)

source(file.path(MISSINGNESS_DIR, "load_missingness.R")); load_missingness(MISSINGNESS_DIR)
# People who belong to more than one cohort: "exclude" (default) leaves them out, "priority" keeps them in the
# first of their cohorts listed in cfg$overlap$priority (Prostate > Diabetes > Insomnia > Healthy).
cfg <- make_config(list(export  = list(date = EXPORT_DATE),      # bucket / owner come from WORKSPACE_BUCKET / OWNER_EMAIL
                        overlap = list(policy = "exclude")))

etl <- run_etl(cfg, L1_DIR, final_pid_path = FINAL_PID, cci_path = CCI_FILE)

# What did the "mostly NA / constant" rule remove? (also written to L1/dropped_columns.tsv)
print(etl$dropped_columns); print(etl$summary)
