# ==============================================================================
# bq_export.R
# ------------------------------------------------------------------------------
# Step 0 of the pipeline (needed only when the exports do not exist, e.g. after a
# workspace migration): run the BigQuery queries and save the results as CSV files in
# the workspace bucket - one export per cohort x stream (4 x 4 = 16).
#
# The SQL is the SAME as in the earlier notebooks (copied column for column), only
# written once instead of 16 times. The destination path comes from export_path(), the
# same function the ETL uses to read, so reading and writing can never disagree.
#
# bigrquery is used only inside bq_executor(); everything else is plain R so the
# SQL and the paths can be tested without BigQuery access.
# ==============================================================================

# ---- SQL -------------------------------------------------------------------------------

# Columns of the simple (one-table) streams, as in the earlier export notebooks.
BQ_STREAM_TABLES <- list(
  activity = list(table = "activity_summary",
                  cols = c("person_id", "date", "activity_calories", "calories_bmr", "calories_out",
                           "elevation", "fairly_active_minutes", "floors", "lightly_active_minutes",
                           "marginal_calories", "sedentary_minutes", "steps", "very_active_minutes")),
  sleep    = list(table = "sleep_daily_summary",
                  cols = c("person_id", "sleep_date", "is_main_sleep", "minute_in_bed", "minute_asleep",
                           "minute_after_wakeup", "minute_awake", "minute_restless", "minute_deep",
                           "minute_light", "minute_rem", "minute_wake")),
  hr       = list(table = "heart_rate_summary",
                  cols = c("person_id", "date", "zone_name", "min_heart_rate", "max_heart_rate",
                           "minute_in_zone", "calorie_count"))
)

# "123, 456, 789" - the person ids of one cohort, as used in the WHERE ... IN (...) clause.
pid_in_clause <- function(person_ids) paste(person_ids, collapse = ", ")

# SQL text for one stream and a set of person ids.
bq_stream_sql <- function(stream, person_ids) {
  ids <- pid_in_clause(person_ids)
  if (stream == "meta") {
    return(paste0("
    SELECT
        person.person_id,
        person.birth_datetime as date_of_birth,
        person.race_concept_id,
        p_race_concept.concept_name as race,
        person.ethnicity_concept_id,
        p_ethnicity_concept.concept_name as ethnicity,
        person.sex_at_birth_concept_id,
        p_sex_at_birth_concept.concept_name as sex_at_birth,
        person.self_reported_category_concept_id,
        p_self_reported_category_concept.concept_name as self_reported_category
    FROM
        `person` person
    LEFT JOIN
        `concept` p_race_concept
            ON person.race_concept_id = p_race_concept.concept_id
    LEFT JOIN
        `concept` p_ethnicity_concept
            ON person.ethnicity_concept_id = p_ethnicity_concept.concept_id
    LEFT JOIN
        `concept` p_sex_at_birth_concept
            ON person.sex_at_birth_concept_id = p_sex_at_birth_concept.concept_id
    LEFT JOIN
        `concept` p_self_reported_category_concept
            ON person.self_reported_category_concept_id = p_self_reported_category_concept.concept_id
    WHERE
        person.person_id IN (", ids, ")"))
  }
  spec <- BQ_STREAM_TABLES[[stream]]
  if (is.null(spec)) stop("unknown stream for export: ", stream, call. = FALSE)
  paste0("
    SELECT
        ", paste0(spec$table, ".", spec$cols, collapse = ",\n        "), "
    FROM
        `", spec$table, "` ", spec$table, "
    WHERE
        ", spec$table, ".PERSON_ID IN (", ids, ")")
}

# ---- running a query and saving it to the bucket ---------------------------------------------

# Default executor: run the SQL in the CDR dataset and write CSV file(s) to `path`
# (the same three calls as in the earlier notebooks, with the dataset and billing project taken
# from cfg$export$cdr / cfg$export$billing, which default to WORKSPACE_CDR / GOOGLE_PROJECT).
bq_executor <- function(sql, path, cfg) {
  cdr <- cfg$export$cdr; billing <- cfg$export$billing
  if (!nzchar(billing)) {
    stop("no billing project: GOOGLE_PROJECT is empty. Set cfg$export$billing <- \"<your google project>\"", call. = FALSE)
  }
  if (!grepl("^[^.[:space:]]+\\.[^.[:space:]]+$", cdr)) {
    stop("the CDR dataset must look like \"project.dataset\" (two parts, one dot) but is '", cdr,
         "'. WORKSPACE_CDR is blank or malformed in this workspace: set cfg$export$cdr <- \"<project>.<dataset>\"",
         call. = FALSE)
  }
  bigrquery::bq_table_save(bigrquery::bq_dataset_query(cdr, sql, billing = billing),
                           path, destination_format = "CSV")
  invisible(path)
}

# The executor used when none is given: bq_executor with the settings of `cfg`.
default_executor <- function(cfg) function(sql, path) bq_executor(sql, path, cfg)

# Export one cohort x stream. Returns a one-row status table.
#   skipped  : files already exist at the destination (set overwrite = TRUE to redo)
#   exported : the query ran and files now exist
export_cohort_stream <- function(cfg, cohort, stream, person_ids, executor = default_executor(cfg),
                                 lister = gsutil_ls, overwrite = FALSE) {
  path <- export_path(cfg, cohort, stream)
  status_row <- function(status, note = "") data.table(cohort = cohort, stream = stream, n_people = length(person_ids),
                                                 path = path, status = status, note = note)
  if (length(person_ids) == 0L) return(status_row("failed", "no person ids for this cohort in final_pid"))

  existing <- tryCatch(lister(path), error = function(e) character())
  if (length(existing) > 0L && !overwrite) {
    log_info("export %s/%s: %d file(s) already there - skipped", cohort, stream, length(existing))
    return(status_row("skipped", paste(length(existing), "existing file(s)")))
  }

  log_info("export %s/%s: running query for %d people", cohort, stream, length(person_ids))
  res <- tryCatch({
    executor(bq_stream_sql(stream, person_ids), path)
    n <- length(lister(path))                     # stops if the export produced no file
    status_row("exported", paste(n, "file(s)"))
  }, error = function(e) status_row("failed", conditionMessage(e)))
  res
}

# Export every cohort x stream. `final_pid` has person_id and cohort. A failure in one
# export does not stop the others, but the function stops at the end if any failed, so a
# partial set of exports is never carried into the ETL.
run_bq_exports <- function(cfg, final_pid, streams = c("activity", "sleep", "hr", "meta"),
                           cohorts = cfg$cohorts, executor = default_executor(cfg), lister = gsutil_ls,
                           overwrite = FALSE) {
  out <- rbindlist(lapply(cohorts, function(co) {
    ids <- final_pid[cohort == co, person_id]
    rbindlist(lapply(streams, function(s) export_cohort_stream(cfg, co, s, ids, executor, lister, overwrite)))
  }))
  print(out[, .(cohort, stream, n_people, status, note)])
  if (any(out$status == "failed")) {
    why <- unique(out[status == "failed", note])
    stop(sum(out$status == "failed"), " export(s) failed. Reason(s): ",
         paste(utils::head(substr(why, 1, 400), 3L), collapse = " || "),
         "\nFix and re-run (finished exports are skipped).", call. = FALSE)
  }
  invisible(out)
}
