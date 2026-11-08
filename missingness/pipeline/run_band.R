# ==============================================================================
# run_band.R - the missingness pipeline on ANOTHER compliance group, and its comparison with the >= 70% run
# (each "# ==== CELL n" block = one notebook cell; run them in order)
#
# Compliance 50-60% in ALL THREE streams, every other filter as in the baseline selection (>= 21 record days and longest gap <= 14 days in
# each stream, CCI <= 10, age and CCI at the first record). The selection itself is the code of ../selection (nothing there or in the model code
# is changed): it only reads the wear metrics, persons and CCI of the full-pool selection (result/selection/20261031) and applies the new band.
# The people of this group are NOT in the >= 70% list, so their 16 exports, day tables and results are new, in their own folders:
#   list      result/selection/<NEW_DATE>/final_pid_base50_60.txt (+ attrition, cci.txt)
#   exports   <bucket>/bq_exports/<owner>/<NEW_DATE>/
#   L1        result/final_missing/L1/5060            results   result/final_missing/results/5060
#   compare   result/final_missing/results/compare_70_vs_5060   (tables/ and figures/)
# Another band: change BAND_TAG, BAND_FILE, COMP_MIN, COMP_MAX, NEW_DATE (a new export date for a new list) and the label of the comparison.
#
# AFTER A KERNEL RESTART: run Cell 1, then the cell you were in. Exports skip folders that exist; every finished model is a checkpoint.
# ==============================================================================
# ==== CELL 1   (setup)
library(data.table)
BASE <- "/home/jupyter/workspace/rw-migration-aou-rw-6f1dba72/workspaces"
ROOT <- file.path(BASE, "temporaldigitalhealthphenotyping")
MISS <- file.path(ROOT, "final_missing")
source(file.path(MISS, "load_missingness.R")); load_missingness(MISS)
source(file.path(ROOT, "functions/cci_function.R"))
source(file.path(ROOT, "selection/selection_sql.R")); source(file.path(ROOT, "selection/selection_core.R"))
source(file.path(MISS, "pipeline/compare_bands.R"))

BAND_TAG  <- "5060"                  # sub-folder of L1 and results
BAND_FILE <- "base50_60"             # final_pid_<BAND_FILE>.txt
COMP_MIN  <- 50; COMP_MAX <- 60      # compliance band [COMP_MIN, COMP_MAX) percent, in every stream
BAND_LABELS <- c(">=70%", "50-60%")  # reference band first
NEW_DATE  <- "20261104"              # export date of this list (the 16 exports); a new list needs a new date
OLD_DATE  <- "20261031"              # the full-pool selection and its exports

SEL_OLD <- file.path(BASE, "result/selection", OLD_DATE)       # wear_metrics.tsv, persons.tsv, cci.txt, final_pid_base.txt
SEL_NEW <- file.path(BASE, "result/selection", NEW_DATE)
FINAL_PID_BAND <- file.path(SEL_NEW, sprintf("final_pid_%s.txt", BAND_FILE))
CCI_BAND       <- file.path(SEL_NEW, "cci.txt")
L1_REF  <- file.path(BASE, "result/final_missing/L1");      OUT_REF  <- file.path(BASE, "result/final_missing/results")
L1_BAND <- file.path(L1_REF, BAND_TAG);                     OUT_BAND <- file.path(OUT_REF, BAND_TAG)
CMP_DIR <- file.path(OUT_REF, "compare_70_vs_5060")

# the two runs use the same settings; only the export date (and the compliance check of the group) differ
cfg_ref  <- make_config(list(export = list(date = OLD_DATE), overlap = list(policy = "exclude")))
cfg_band <- make_config(list(export = list(date = NEW_DATE), overlap = list(policy = "exclude"),
                             selection = list(min_pct_observed = COMP_MIN)))       # `selection_check` compares with the band's lower limit
# If the CDR is blank or not "project.dataset" (Cell 3 below needs it), set cfg_band$export$cdr and $billing as in notebook_cells.R Cell 2.
stopifnot("full-pool selection files not found" = all(file.exists(file.path(SEL_OLD, c("wear_metrics.tsv", "persons.tsv", "cci.txt", "final_pid_base.txt")))),
          "the >= 70% results are not there (run the pipeline on the baseline list first)" = dir.exists(file.path(OUT_REF, "tables")),
          "WORKSPACE_BUCKET / OWNER_EMAIL not set" = nzchar(cfg_band$export$bucket) && nzchar(cfg_band$export$owner))

# ==== CELL 2   (the people of the band: the baseline rules with only the compliance band changed; no BigQuery)
wear <- fread(file.path(SEL_OLD, "wear_metrics.tsv")); cci <- fread(file.path(SEL_OLD, "cci.txt"))
pool <- fread(file.path(SEL_OLD, "persons.tsv"))[, .(cohort, person_id)]
sel  <- sel_select(wear, cci, pool, comp_min = COMP_MIN, comp_max = COMP_MAX)
dir.create(SEL_NEW, recursive = TRUE, showWarnings = FALSE)
sel_write(sel, SEL_NEW, BAND_FILE)
file.copy(file.path(SEL_OLD, "cci.txt"), CCI_BAND, overwrite = TRUE)         # the CCI of the whole pool (age at the first record) covers these people
print(sel$attrition)
ref_pid <- fread(file.path(SEL_OLD, "final_pid_base.txt"))
print(sel_compare(sel$pids, ref_pid)[, .(cohort, n_band = n_new, n_ref = n_old, n_in_both = n_both)])
stopifnot("the two lists share a person: the compliance bands must not overlap" = !any(sel$pids$person_id %in% ref_pid$person_id))
final_pid <- read_final_pid(FINAL_PID_BAND)
ov <- resolve_cohort_overlap(final_pid, "exclude", cfg_band$overlap$priority)
cat("memberships:", nrow(final_pid), "| unique people:", ov$n_unique_people, "| in more than one cohort:", length(ov$overlap_ids), "| kept (policy exclude):", ov$n_kept, "\n")
print(ov$overlap)
# A cohort with few people may give "not_estimable" models (fewer than cfg$model$min_events events or min_persons_with_events people with events).

# ==== CELL 3   (BigQuery: the 16 exports of these people + the device version; finished exports are skipped)
exports     <- run_bq_exports(cfg_band, final_pid)
dev_exports <- run_bq_exports(cfg_band, final_pid, streams = "device")

# ==== CELL 4   (ETL: raw exports -> day tables of this group in L1_BAND)
dir.create(L1_BAND, recursive = TRUE, showWarnings = FALSE)
etl <- run_etl(cfg_band, L1_BAND, final_pid_path = FINAL_PID_BAND, cci_path = CCI_BAND)
print(etl$summary); print(etl$persons[, .N, by = cohort])
dev <- run_etl_device(cfg_band, L1_BAND); print(dev$summary)
rm(etl, dev); invisible(gc())

# ==== CELL 5   (missingness description of the group: the tables the comparison reads; ~ minutes)
inp_b <- load_inputs(cfg_band, L1_BAND)
res <- run_missingness_pipeline(inp_b$day_tables, inp_b$persons, cfg_band, OUT_BAND, stages = c("missingness", "temporal"))
print(res$cohort_missingness); print(res$selection_check)       # selection_check: n below the band's limit should be 0

# ==== CELL 6   (the PRIMARY onset model, three streams; one stream at a time, every model a checkpoint)
res <- run_missingness_pipeline(inp_b$day_tables, inp_b$persons, cfg_band, OUT_BAND, stages = "onset",
                                onset_roles = "primary", model_streams = c("activity", "sleep", "hr"))
print(res$primary_test)
print(res$onset_models[role == "primary", .(stream, n_rows, n_persons, n_events, n_persons_event, status)])

# ==== CELL 7   (the sensitivity analyses of the primary model, incl. the 16 / 20 / 22 h wear analyses; the long cell)
res <- run_missingness_pipeline(inp_b$day_tables, inp_b$persons, cfg_band, OUT_BAND, stages = "sensitivity",
                                model_streams = c("activity", "sleep", "hr"))
print(res$sensitivity[, .(stream, analysis, term, or = round(or, 2), p = signif(p, 3), n_events, status)])
cat("failed units:", nrow(res$unit_failures), "\n")                 # 0 = nothing failed; otherwise run this cell again to retry them

# ==== CELL 8   (the comparison: tables and figures with both groups; reads the saved tables of the two runs, no model is refitted)
# figures: compare_primary_effect (the primary term per stream), compare_sensitivity (every sensitivity analysis of both groups),
#          compare_missing_pct / compare_longest_gap / compare_gaps_per_person (person level), compare_gap_length_distribution, compare_transition
# tables : compare_primary (estimate, SE, OR and size of each model), compare_primary_ratio (ratio of the two odds ratios with a 95% interval;
#          independent groups, no p-value), compare_within_sd (what "+1 SD" is in each run), compare_sensitivity, compare_missingness
out_dirs <- setNames(list(OUT_REF, OUT_BAND), BAND_LABELS)
stores   <- setNames(list(open_day_store(cfg_ref, L1_REF), open_day_store(cfg_band, L1_BAND)), BAND_LABELS)
cmp <- cmp_run(out_dirs, cfg_ref, CMP_DIR, stores)
print(cmp$ratio[, .(stream, n_events_ref, n_events_cmp, or_ref = round(or_ref, 3), or_cmp = round(or_cmp, 3),
                    ratio = round(ratio, 3), low = round(ratio_low, 3), high = round(ratio_high, 3),
                    sd_ref = round(within_sd_ref, 3), sd_cmp = round(within_sd_cmp, 3))])
print(cmp$missingness)
show_figure(CMP_DIR, "compare_primary_effect")
