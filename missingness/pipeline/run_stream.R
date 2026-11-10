# ==============================================================================
# run_stream.R - which streams are present on each day (the five states), Parts A-G; see STREAM.md
# (each "# ==== CELL n" block = one notebook cell; run them in order)
#
# Uses the people of the existing selection and the day tables already in L1 (nothing is re-selected, nothing of the other analyses is changed); the only
# new input is the RAW device export (clip-on and phone trackers are not in the L1 watch timeline). Descriptive: estimates with person-clustered intervals, no tests.
# Output folder: <results>/stream/  (tables/, figures/, methods.md; the check without unknown-device days in stream/no_unknown_device/).
# Parts A-C run first and end with a gate: if the date check (Part C) or the minutes check (Part G) fails, Cell 4 stops and parts D-G are not run.
# ==============================================================================
# ==== CELL 1   (setup)
library(data.table)
BASE <- "/home/jupyter/workspace/rw-migration-aou-rw-6f1dba72/workspaces"
MISS <- file.path(BASE, "temporaldigitalhealthphenotyping/final_missing")
source(file.path(MISS, "load_missingness.R")); load_missingness(MISS)
source(file.path(MISS, "R/stream_states.R")); source(file.path(MISS, "R/plots_stream_states.R"))
EXPORT_DATE <- "20261031"                                                   # the export of the people of the selection (as in notebook_cells.R)
L1_DIR      <- file.path(BASE, "result/final_missing/L1")
OUT_DIR     <- file.path(BASE, "result/final_missing/results")
OUT_STREAM  <- file.path(OUT_DIR, "stream"); dir.create(OUT_STREAM, recursive = TRUE, showWarnings = FALSE)
cfg <- make_config(list(export = list(date = EXPORT_DATE), overlap = list(policy = "exclude")))
st  <- st_config(list(seed = cfg$seed))             # settings of this analysis (cut-offs, gate limits, the list of versions without heart rate)
stopifnot("WORKSPACE_BUCKET / OWNER_EMAIL not set" = nzchar(cfg$export$bucket) && nzchar(cfg$export$owner))

# ==== CELL 2   (people, day tables and the devices of every day; the raw device export is read once)
inp <- load_inputs(cfg, L1_DIR); store <- inp$day_tables
persons <- copy(inp$persons); dp <- store$device_person()
stopifnot("run Cell 8b of notebook_cells.R first (device_person.rds is missing)" = !is.null(dp))
persons[dp, on = "person_id", device_family := i.device_family]
stopifnot("the person table needs cohort and age_group" = all(c("cohort", "age_group") %in% names(persons)))
raw_dev <- rbindlist(lapply(cfg$cohorts, function(co) read_export_files(export_path(cfg, co, "device"), c(device_version = "character"))), use.names = TRUE, fill = TRUE)
trk <- st_trackers(raw_dev, persons, cfg); rm(raw_dev)
lookup <- st_device_lookup(trk, st)                 # CHECK THIS TABLE: which versions are treated as having no heart-rate sensor (product knowledge, not data)
dir.create(file.path(OUT_STREAM, "tables"), showWarnings = FALSE); fwrite(lookup, file.path(OUT_STREAM, "tables", "stream_device_hr_lookup.tsv"), sep = "\t")
print(lookup, nrows = 60)

# ==== CELL 3   (the person-day table: one row per person-day of every window; ~ a minute)
D <- st_build_days(store, persons, trk, st, cfg)
writeLines(st_methods_text(st), file.path(OUT_STREAM, "methods.md"))
print(D[, .N, by = .(cohort, reason)][order(cohort, reason)])

# ==== CELL 4   (PARTS A-C and the gate: data flow, the activity filter, patterns and nesting, date check, minutes check)
abc <- st_run_abc(D, st)
st_write_tables(abc, OUT_STREAM); st_save_figures(st_figures(abc, NULL, cfg, st), OUT_STREAM)
print(abc$flow[, .(cohort, step, days = n_days_at_step, share_pct = round(share_pct, 2), low = round(low, 2), high = round(high, 2))], nrows = 60)
print(abc$activity_filter[, .(cohort, step, n_days, n_events, share_pct = round(share_pct, 3))])
print(abc$combos[stage == "after the filter", .(cohort, label, category, share_pct = round(share_pct, 3))], nrows = 40)
print(abc$gate); print(abc$date_check)
st_require_gate(abc$gate)                           # stops here if the date check or the minutes check fails: fix the data preparation first

# ==== CELL 5   (PARTS D-G: what state 3 is, composition, what 'night missing' contains, partial wear)
dg <- st_run_dg(D, st)
st_write_tables(dg, OUT_STREAM); st_save_figures(st_figures(abc, dg, cfg, st), OUT_STREAM)
print(dg$state3_vs_state1); print(dg$night_habit); print(dg$probably_worn); print(dg$hr_wear_below)

# ==== CELL 6   (the check: the same analysis WITHOUT every day of unknown device; results in stream/no_unknown_device/)
st_chk <- modifyList(st, list(drop_unknown_days = TRUE))
D2 <- st_assign_reason(copy(D), st_chk)
chk <- st_run_variant(D2, st_chk, cfg, file.path(OUT_STREAM, "no_unknown_device"))
print(chk$abc$flow[step == "analysed days", .(cohort, days = n_days_at_step, share_pct = round(share_pct, 2))])
