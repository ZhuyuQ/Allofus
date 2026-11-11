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

# ==== CELL 4b   (ONLY IF THE GATE FAILED ON THE MINUTES: what the records above the length of the day look like; read-only, nothing is changed)
# ratio of the minutes to the length of the day (about 2 = rows counted twice), how concentrated in people, in which years
diag <- st_minutes_diagnosis(D); print(diag)
# the raw heart-rate export, one cohort at a time: the same person-days under three ways of combining the rows of a zone
tot <- rbindlist(lapply(cfg$cohorts, function(co) {
  raw <- read_export_files(export_path(cfg, co, "hr"), c(zone_name = "character"))
  r <- st_hr_raw_totals(raw[person_id %in% persons$person_id], cfg$etl$hr_zones); rm(raw); invisible(gc()); r
}))
print(st_hr_raw_check(tot))                       # the rule with ~0% above the length of the day is the right way to combine the rows

# ==== CELL 4c   (ONLY IF THE GATE FAILED ON THE MINUTES; choose ONE option from what Cell 4b showed, then run Cell 4 again; L1 and the other analyses are untouched)
# OPTION 1 - rows counted twice (Cell 4b: ratios near 2, one combining rule removes the problem): use the corrected minutes of that rule ("t_distinct" or "t_max")
# D <- st_use_hr_minutes(D, tot, "t_distinct", st)
# OPTION 2 - the minutes are slightly above the length of the day in the raw export itself (Cell 4b: ratios 1.0-1.1, no duplicated zone rows, spread over many people): the day was longer
#   than 24 hours (travel across time zones, daylight-saving dates outside the US rule). Records up to 10% above the length of the day are accepted by the check, and the heart-rate wear of Parts F and G is
#   capped at the length of the day. Records further above still count against the 0.1% limit.
st <- modifyList(st, list(minutes_tolerance = 0.1, cap_hr_minutes = TRUE))
D <- st_assign_reason(st_flags(D, st), st)
fwrite(st_hr_raw_check(tot), file.path(OUT_STREAM, "tables", "stream_hr_minutes_raw_check.tsv"), sep = "\t")
# then run Cell 4 again: the gate must pass before Cell 5

# ==== CELL 5   (PARTS D-G: what state 3 is, composition, what 'night missing' contains, partial wear)
dg <- st_run_dg(D, st)
st_write_tables(dg, OUT_STREAM); st_save_figures(st_figures(abc, dg, cfg, st), OUT_STREAM)
print(dg$state3_vs_state1); print(dg$night_habit); print(dg$probably_worn); print(dg$hr_wear_below)

# ==== CELL 6   (the check: the same analysis WITHOUT every day of unknown device; results in stream/no_unknown_device/)
st_chk <- modifyList(st, list(drop_unknown_days = TRUE))
D2 <- st_assign_reason(copy(D), st_chk)
chk <- st_run_variant(D2, st_chk, cfg, file.path(OUT_STREAM, "no_unknown_device"))
print(chk$abc$flow[step == "analysed days", .(cohort, days = n_days_at_step, share_pct = round(share_pct, 2))])

# ==============================================================================
# STREAM MODEL (STREAM.md, "Stream model"): does a high-step day d predict the status of day d + 1?  Run AFTER Cell 4 passed (D = the person-day table, with the minutes option of Cell 4c if it was needed).
# Model A: one mixed logistic model per category (1-4 against 0); model B: one multinomial logit with person-clustered errors. Every fit is a checkpoint; all rows are used (nothing subsampled).
# ==============================================================================
# ==== CELL 7   (setup of the stream model, risk set, flow and event counts - look at the event counts BEFORE the models)
source(file.path(MISS, "R/stream_model.R")); source(file.path(MISS, "R/plots_stream_model.R"))
sc <- stm_config(list(seed = cfg$seed)); OUT_MODEL <- file.path(OUT_STREAM, "model"); dir.create(OUT_MODEL, showWarnings = FALSE)
CKPT <- file.path(OUT_DIR, cfg$checkpoint$dirname)
rs <- stm_risk_set(stm_days(D, sc), persons, cfg, sc); R <- rs$R
dkey <- list(export = EXPORT_DATE, n_rows = nrow(R), n_people = uniqueN(R$person_id))
print(rs$flow); print(rs$scales)
ev <- list(main = stm_event_counts(R, sc)); print(ev$main[group == "all"])
print(stm_cat4_composition(R))                      # which original patterns make up category 4

# ==== CELL 8   (main run: model A for categories 1-4 and the secondary split of category 3, model B; the longest cell - a checkpoint per model, so it can be re-run after a restart)
out_main <- stm_run_models(R, "main", cfg, sc, CKPT, dkey, secondary = TRUE, multinomial = TRUE)
res_main <- stm_collect(out_main, sc, "main")
print(res_main$tests[, .(model, or, or_low, or_high, p_wald, p_tost, p_tost_holm, equivalent, negative_control, n_events)])
print(res_main$contrasts); print(res_main$risk_diff)

# ==== CELL 9   (sensitivity: more than 16 h / 22 h of heart-rate wear on d, first 14 days excluded, first transition only; model A only)
runs <- list(wear_ge16h = stm_risk_set(stm_days(D, stm_config(list(wear_hours = 16, seed = cfg$seed))), persons, cfg, sc)$R,
             wear_ge22h = stm_risk_set(stm_days(D, stm_config(list(wear_hours = 22, seed = cfg$seed))), persons, cfg, sc)$R,
             first14 = stm_variant_rows(R, "first14", sc), first_transition = stm_variant_rows(R, "first_transition", sc))
res_var <- lapply(names(runs), function(v) { ev[[v]] <<- stm_event_counts(runs[[v]], sc)
  stm_collect(stm_run_models(runs[[v]], v, cfg, sc, CKPT, c(dkey, list(n = nrow(runs[[v]]))), secondary = FALSE, multinomial = FALSE), sc, v) })
names(res_var) <- names(runs)

# ==== CELL 10   (cohort-specific fits, only for cohorts with enough category-1 events; checks of the engine)
big <- ev$main[category == 1L & group != "all" & events >= sc$min_events_effect, group]
res_coh <- lapply(big, function(co) stm_collect(stm_run_models(R[cohort == co], paste0("cohort_", co), cfg, sc, CKPT, c(dkey, list(co = co)), secondary = FALSE, multinomial = FALSE), sc, paste0("cohort_", co)))
checks <- list(clogit = stm_clogit(R, sc), engine = stm_verify_engine(R, cfg, sc)); print(checks)

# ==== CELL 11   (tables and figures)
all_res <- c(list(main = res_main), res_var, setNames(res_coh, paste0("cohort_", big)))
bindr <- function(nm) rbindlist(lapply(all_res, `[[`, nm), fill = TRUE)
res <- list(info = bindr("info"), coefficients = bindr("coefficients"), tests = bindr("tests"), multinomial_coefficients = res_main$multinomial_coefficients,
            contrasts = res_main$contrasts, risk_curve = res_main$risk_curve, risk_diff = res_main$risk_diff)
ec <- rbindlist(lapply(names(ev), function(v) ev[[v]][, variant := v]))
figs <- stm_figures(res, stm_deciles(R, sc), sc, ec)
stm_write(c(res, list(event_counts = ec, flow = rs$flow, scales = rs$scales, category4_composition = stm_cat4_composition(R), checks_clogit = checks$clogit, checks_engine = checks$engine)), figs, OUT_MODEL)
