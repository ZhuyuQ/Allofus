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

# Optional extra export (not part of the default run): the Fitbit device table, one row per person x device x day.
# Columns asked for: the device version, plus the keys needed to attach it to a person and a day, and device_type
# (the table also contains Fitbit scales, which are not wearables).
# Whether the CDR has this table and these columns is checked first with bq_probe_device_table(); whether every person has a
# device version is checked with bq_device_coverage().
BQ_DEVICE_TABLE <- "device"
BQ_DEVICE_COLS  <- c("person_id", "device_date", "device_version", "device_type")   # device_type: the table also holds scales (SCALE)

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
  if (stream == "device") {
    return(paste0("
    SELECT
        ", paste0(BQ_DEVICE_TABLE, ".", BQ_DEVICE_COLS, collapse = ",\n        "), "
    FROM
        `", BQ_DEVICE_TABLE, "` ", BQ_DEVICE_TABLE, "
    WHERE
        ", BQ_DEVICE_TABLE, ".PERSON_ID IN (", ids, ")"))
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

# Is the device export feasible in this CDR? Reads 5 rows of the device table and reports which of the wanted columns exist.
# Nothing is written. `sampler(sql)` returns a data.frame; the default runs the query in the CDR (billing project and dataset
# from cfg$export, as bq_executor does).
bq_probe_device_table <- function(cfg, columns = BQ_DEVICE_COLS, sampler = NULL) {
  if (is.null(sampler)) sampler <- function(sql) {
    cdr <- cfg$export$cdr; billing <- cfg$export$billing
    if (!nzchar(billing)) stop("no billing project: set cfg$export$billing", call. = FALSE)
    bigrquery::bq_table_download(bigrquery::bq_dataset_query(cdr, sql, billing = billing), n_max = 5L)
  }
  s <- sampler(paste0("SELECT * FROM `", BQ_DEVICE_TABLE, "` LIMIT 5"))
  have <- tolower(names(s))
  out <- data.table(column = columns, in_table = tolower(columns) %in% have)
  list(columns_check = out, all_found = all(out$in_table), table_columns = names(s), sample = utils::head(as.data.frame(s), 5L))
}

# Does every person have a device version? One query per cohort, nothing is written.
# The CDR device table holds Fitbit SCALES (device_type SCALE, e.g. "Aria Air") as well as wearable trackers, so a person who only owns
# a scale has a device_version but no wearable. Everything below that says "version" counts non-scale devices only.
#   n_people                  people of the cohort in final_pid
#   n_with_activity           ... with at least one activity_summary day (the people the analysis can use)
#   n_with_device_row         ... with at least one row in the device table (any device type)
#   n_with_version            ... with a non-empty device_version on a non-scale device
#   n_activity_with_version   people with activity AND such a device version (share: pct_activity_with_version, the number that matters)
#   n_activity_days / n_activity_days_with_version   activity days, and those with a non-scale device_version row on the same date
#                             (rows without a device_date cannot be matched to a day, so this share is a lower bound)
bq_device_coverage_sql <- function(person_ids) {
  ids <- pid_in_clause(person_ids)
  paste0("
    WITH a AS (
        SELECT person_id, CAST(date AS DATE) AS day FROM `activity_summary` WHERE person_id IN (", ids, ")
    ), dp AS (
        SELECT person_id, CAST(device_date AS DATE) AS day FROM `", BQ_DEVICE_TABLE, "`
        WHERE person_id IN (", ids, ") AND device_version IS NOT NULL AND TRIM(CAST(device_version AS STRING)) != ''
          AND (device_type IS NULL OR UPPER(CAST(device_type AS STRING)) != 'SCALE')
    ), d AS (
        SELECT DISTINCT person_id, day FROM dp WHERE day IS NOT NULL
    )
    SELECT
        (SELECT COUNT(DISTINCT person_id) FROM a) AS n_with_activity,
        (SELECT COUNT(DISTINCT person_id) FROM `", BQ_DEVICE_TABLE, "` WHERE person_id IN (", ids, ")) AS n_with_device_row,
        (SELECT COUNT(DISTINCT person_id) FROM dp) AS n_with_version,
        (SELECT COUNT(DISTINCT a.person_id) FROM a JOIN dp ON a.person_id = dp.person_id) AS n_activity_with_version,
        (SELECT COUNT(*) FROM a) AS n_activity_days,
        (SELECT COUNT(*) FROM a JOIN d ON a.person_id = d.person_id AND a.day = d.day) AS n_activity_days_with_version")
}

# Which devices are there? People and rows per device_type and device_version (all cohorts together), most common first.
bq_device_versions_sql <- function(person_ids) {
  paste0("
    SELECT device_type, device_version, COUNT(DISTINCT person_id) AS n_people, COUNT(*) AS n_rows
    FROM `", BQ_DEVICE_TABLE, "`
    WHERE person_id IN (", pid_in_clause(person_ids), ")
    GROUP BY device_type, device_version
    ORDER BY n_people DESC")
}

bq_device_versions <- function(cfg, final_pid, sampler = NULL) {
  if (is.null(sampler)) sampler <- default_bq_sampler(cfg)
  as.data.table(sampler(bq_device_versions_sql(unique(final_pid$person_id))))
}

default_bq_sampler <- function(cfg) function(sql) {
  cdr <- cfg$export$cdr; billing <- cfg$export$billing
  if (!nzchar(billing)) stop("no billing project: set cfg$export$billing", call. = FALSE)
  bigrquery::bq_table_download(bigrquery::bq_dataset_query(cdr, sql, billing = billing))
}

bq_device_coverage <- function(cfg, final_pid, cohorts = cfg$cohorts, sampler = NULL) {
  if (is.null(sampler)) sampler <- default_bq_sampler(cfg)
  out <- rbindlist(lapply(cohorts, function(co) {
    ids <- final_pid[cohort == co, person_id]
    r <- as.data.table(sampler(bq_device_coverage_sql(ids)))
    cbind(data.table(cohort = co, n_people = length(ids)), r[, lapply(.SD, as.numeric)])
  }))
  out[, `:=`(pct_activity_with_version = round(100 * n_activity_with_version / n_with_activity, 1),
             pct_days_with_version     = round(100 * n_activity_days_with_version / n_activity_days, 1))]
  out[]
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
