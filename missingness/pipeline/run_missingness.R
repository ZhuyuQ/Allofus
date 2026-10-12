# ==============================================================================
# run_missingness.R   -  the analysis (reads L1, writes L2)
# ------------------------------------------------------------------------------
# Stages (any subset, in any order of re-running):
#   missingness  cohort tables + comparisons + figures
#   temporal     onset vs continuation, study week, calendar year, day of week
#   descriptive  summary statistics of the metric values by cohort and strata
#   onset        gap-onset mixed models (primary, model A, model D, exploratory)
#   gap_models   models B and C (use gap length)
#   sensitivity  true sensitivity analyses of the primary model
# ==============================================================================
BASE            <- "/home/jupyter/workspace/rw-migration-aou-rw-6f1dba72/workspaces"
MISSINGNESS_DIR <- file.path(BASE, "temporaldigitalhealthphenotyping/missingness")       # <- where you put the folder
L1_DIR          <- file.path(BASE, "result/final/missingness/L1")
OUT_DIR         <- file.path(BASE, "result/final/missingness/results")

STAGES      <- c("missingness", "temporal", "descriptive", "onset", "gap_models", "sensitivity")
ONSET_ROLES <- c("primary", "secondary_threshold", "exploratory")
# Speed knob: fit only the exploratory (one-predictor / one-set) models on a seeded sample of
# people. The primary, model-D and sensitivity fits always use everyone. NULL = use everyone.
EXPLORATORY_SAMPLE_PERSONS <- 3000L

source(file.path(MISSINGNESS_DIR, "load_missingness.R")); load_missingness(MISSINGNESS_DIR)
cfg <- make_config(list(model = list(exploratory_sample_persons = EXPLORATORY_SAMPLE_PERSONS)))
# Resources (defaults suit 4 CPU cores / 32 GB RAM): cfg$resources$memory_budget_gb, memory_reserve_gb, light_workers
inp <- load_inputs(cfg, L1_DIR)
# The model stages (onset, gap_models, sensitivity) handle one stream at a time and save every finished
# model in OUT_DIR/checkpoints, so a crashed run can simply be started again: finished models are skipped,
# missing / failed / outdated ones are run. To run one stream per call:
#   run_missingness_pipeline(..., stages = "onset", model_streams = "activity")
res <- run_missingness_pipeline(inp$day_tables, inp$persons, cfg, OUT_DIR,
                                stages = STAGES, onset_roles = ONSET_ROLES)

print(res$cohort_missingness)
cat("failed units:", nrow(res$unit_failures), "\n")   # 0 = no model failed (otherwise start this script again to retry them)
print(res$primary_test)    # the pre-specified test: within-person prior-day steps, one row per stream, Holm-adjusted
