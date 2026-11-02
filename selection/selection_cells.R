# ==============================================================================
# selection_cells.R - the notebook cells, in order (each "# ==== CELL n" block = one notebook cell)
# Full-pool selection on CDR C2024Q3R9: >= 21 days, longest gap <= 14 days, compliance >= 70% in activity,
# sleep AND hr, CCI <= 10 with age at the person's first record. Nothing in missingness/ is changed.
# Upload the repo folder Allofus/ (missingness/, function/, selection/) to ROOT below.
# After a crash or kernel restart: run Cells 1-2 again (the CDR is an environment variable of the R session), then the cell you were in; finished queries are read from OUT_DIR.
# Another compliance band, other rules: Cell 8. Export of the selected people: Cell 9.
# ==============================================================================
# ==== CELL 1   (CDR and bucket of this workspace)
cdr    <- trimws(tail(system("wb resource resolve --id=C2024Q3R9", intern = TRUE), 1))
bucket <- trimws(tail(system("wb resource resolve --id=rw-migration-aou-rw-6f1dba72", intern = TRUE), 1))
Sys.setenv(WORKSPACE_CDR = cdr, WORKSPACE_BUCKET = bucket)
cat("CDR   :", cdr, "\nbucket:", bucket, "\n")        # wb-silky-artichoke-2408.C2024Q3R9 | gs://rw-migration-aou-rw-6f1dba72

# ==== CELL 2   (also run after every kernel restart: Sys.setenv of Cell 1 lasts only for the R session)
library(data.table)
stopifnot("run Cell 1 first: WORKSPACE_CDR is '' or malformed" = grepl("^[^.[:space:]]+\\.[^.[:space:]]+$", Sys.getenv("WORKSPACE_CDR")))
BASE        <- "/home/jupyter/workspace/rw-migration-aou-rw-6f1dba72/workspaces"
ROOT        <- file.path(BASE, "temporaldigitalhealthphenotyping")           # = Allofus/ of the repo
EXPORT_DATE <- "20261031"
OUT_DIR     <- file.path(BASE, "result/selection", EXPORT_DATE)
OLD_PID     <- file.path(BASE, "result/final/final_pid.txt")                 # the earlier selection (R8), for comparison
dir.create(file.path(OUT_DIR, "units"), recursive = TRUE, showWarnings = FALSE)
source(file.path(ROOT, "final_missing/R/utils.R"))                            # age_on(), log_info()
source(file.path(ROOT, "functions/cci_function.R"))                          # cci_calc_func()
source(file.path(ROOT, "selection/selection_sql.R"))
source(file.path(ROOT, "selection/selection_core.R"))
stopifnot("GOOGLE_PROJECT is empty: set it (billing project)" = nzchar(Sys.getenv("GOOGLE_PROJECT")))

# ==== CELL 3   (the pool of every cohort and date of birth)
persons <- rbindlist(lapply(names(SEL_COHORTS), function(co) {
  log_info("pool %s", co); x <- sel_fetch(sel_sql_person(co)); x[, cohort := co]
}))
persons[, `:=`(dob = as.Date(date_of_birth), date_of_birth = NULL)]
print(persons[, .N, by = cohort])

# ==== CELL 4   (wear metrics, SQL; one file per cohort x stream)
wear <- rbindlist(lapply(names(SEL_COHORTS), function(co) rbindlist(lapply(names(SEL_STREAMS), function(s) {
  f <- file.path(OUT_DIR, "units", sprintf("wear_%s_%s.tsv", co, s))
  if (!file.exists(f)) { log_info("wear %s / %s", co, s); fwrite(sel_fetch(sel_sql_wear(co, s)), f, sep = "\t") }
  x <- fread(f); x[, `:=`(cohort = co, stream = s)]
}))))
wear <- sel_prepare_wear(wear)
fwrite(wear, file.path(OUT_DIR, "wear_metrics.tsv"), sep = "\t")
print(wear[, .(n_people = .N, median_days = as.numeric(median(n_days)), median_pct_wear = round(median(pct_wear), 1)),
           by = .(cohort, stream)])

# ==== CELL 5   (check: the SQL result equals the metrics recomputed from the raw dates of a random sample)
sel_verify(wear, function(s, ids) sel_fetch(sel_sql_dates(s, ids)), n = 30L)

# ==== CELL 6   (age at the first record over the three streams, conditions, CCI)
persons <- sel_persons(persons, wear)
print(persons[, .(n = .N, no_record = sum(is.na(first_date)), median_age = as.numeric(median(age_first, na.rm = TRUE))), by = cohort])
codes <- sel_cci_codes(); cat(length(codes), "SNOMED codes in the CCI lookup\n")
conds <- rbindlist(lapply(names(SEL_COHORTS), function(co) {
  log_info("conditions %s", co); x <- sel_fetch(sel_sql_conditions(co, codes)); x[, cohort := co]
}))
cci <- sel_cci(persons, conds)
fwrite(persons, file.path(OUT_DIR, "persons.tsv"), sep = "\t")
fwrite(cci, file.path(OUT_DIR, "cci.txt"), sep = "\t", quote = FALSE)       # same columns as pre_final_cci.txt
print(cci[, .(n = .N, median_cci = median(cci_score), n_above_10 = sum(cci_score > 10)), by = cohort])

# ==== CELL 7   (the selection with the default rules, and its comparison with the earlier final_pid.txt)
pool <- persons[, .(cohort, person_id)]
sel  <- sel_write(sel_select(wear, cci, pool), OUT_DIR, "base")
print(sel$attrition)
old <- fread(OLD_PID)
print(sel_compare(sel$pids, old))
# Do the earlier people pass the wear rules (CCI left out)? only_old = people the wear rules do not reproduce.
print(sel_compare(sel_select(wear, cci, pool, cci_max = Inf)$pids, old)[, .(cohort, n_old, only_old)])

# ==== CELL 8   (another rule set: only the argument you name changes; example = compliance 50-60%, i.e. 50 <= % < 60)
# sel_write(sel_select(wear, cci, pool, comp_min = 50, comp_max = 60), OUT_DIR, "comp50-60")
# other arguments: min_days = 21, gap_max = 14, comp_min = 70, comp_max = Inf, cci_max = 10

# ==== CELL 9   (export the selected people; the missingness pipeline is only read)
source(file.path(ROOT, "final_missing/load_missingness.R")); load_missingness(file.path(ROOT, "final_missing"))
cfg <- make_config(list(export = list(date = EXPORT_DATE)))
final_pid <- read_final_pid(file.path(OUT_DIR, "final_pid_base.txt"))
print(final_pid[, .N, by = cohort])
exports <- run_bq_exports(cfg, final_pid)       # to <bucket>/bq_exports/<owner>/20261031/
# The ETL then runs with EXPORT_DATE "20261031", FINAL_PID = final_pid_base.txt and CCI_FILE = cci.txt of OUT_DIR.

# ==== CELL 10   (summary tables and Venn diagram of the selection in OUT_DIR; needs Cells 1-2 and the files of Cells 4-7)
# tables/table1_demographics_<tag>_{all,excl_overlap}.tsv, tables/table_fitbit_<tag>_{all,excl_overlap}.tsv, tables/venn_counts_<tag>.tsv,
# figures/venn_<tag>.png / .pdf.   all = every cohort membership (N = rows of final_pid); excl_overlap = people in several cohorts removed.
source(file.path(ROOT, "final_missing/R/persons.R"))                        # race_group_from_name()
source(file.path(ROOT, "selection/summary_tables.R"))
TAG  <- "base"                                                              # or "comp50-60" etc.
meta <- sel_fetch(st_sql_meta(unique(fread(file.path(OUT_DIR, sprintf("final_pid_%s.txt", TAG)))$person_id)))
res  <- st_run(OUT_DIR, TAG, meta)
print(res$all$table1); print(res$all$fitbit); print(res$venn)
