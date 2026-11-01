# ==============================================================================
# selection_sql.R
# ------------------------------------------------------------------------------
# SQL (CDR C2024Q3R9) for the full-pool selection. Pure string builders, plus sel_fetch(),
# the only function that talks to BigQuery.
#   pool        : the cohort definitions of the dataset-builder exports (export/*.ipynb)
#   wear        : per person and stream - distinct days, first / last date, longest interior gap
#   conditions  : (person, SNOMED code) for the codes of the CCI lookup
# ==============================================================================

SEL_STREAMS <- list(
  activity = list(table = "activity_summary",    date_col = "date"),
  sleep    = list(table = "sleep_daily_summary", date_col = "sleep_date"),
  hr       = list(table = "heart_rate_summary",  date_col = "date")
)

# Fitbit data flags every cohort needs (Prostate also needs has_fitbit_device_information).
SEL_FLAGS <- c("has_fitbit", "has_fitbit_activity_summary", "has_fitbit_plus_device_information",
               "has_fitbit_heart_rate_summary", "has_fitbit_heart_rate_level", "has_fitbit_steps_intraday",
               "has_fitbit_sleep_daily_summary", "has_fitbit_sleep_level")

# include / exclude : parent concepts (the concept and its descendants) a person must / must not have.
SEL_COHORTS <- list(
  Diabetes = list(flags = SEL_FLAGS, include = 201820),
  Insomnia = list(flags = SEL_FLAGS, include = c(436962, 42538607, 37110488, 4102985, 439013)),
  Prostate = list(flags = c(SEL_FLAGS, "has_fitbit_device_information"), include = c(4163261, 200970),
                  exclude_sex = 45878463),                       # Female
  Healthy  = list(flags = SEL_FLAGS,
                  exclude = c(135930, 4101796, 4117930, 255919, 444112, 432453, 4024567, 4274025, 376106,
                              4011630, 4213101, 4038502, 4093991, 40480457, 4041284, 4042141, 43021226,
                              4041283, 4227253, 4023995, 4102111, 138239, 4293175, 4212577, 4302537, 4021667))
)

sel_num <- function(x) paste(sprintf("%.0f", x), collapse = ", ")

# person_ids with an event in any descendant of the parent concepts (dataset-builder "rank1" expansion).
sel_sql_events <- function(parent_ids) paste0("
    SELECT DISTINCT person_id FROM `cb_search_all_events`
    WHERE is_standard = 1 AND concept_id IN (
        SELECT DISTINCT c.concept_id FROM `cb_criteria` c
        JOIN (SELECT CAST(cr.id AS STRING) AS id FROM `cb_criteria` cr
              WHERE cr.concept_id IN (", sel_num(parent_ids), ") AND cr.full_text LIKE '%_rank1]%') a
          ON (c.path LIKE CONCAT('%.', a.id, '.%') OR c.path LIKE CONCAT('%.', a.id)
              OR c.path LIKE CONCAT(a.id, '.%') OR c.path = a.id)
        WHERE c.is_standard = 1 AND c.is_selectable = 1)")

# One column of person_ids: everyone who belongs to the cohort.
sel_sql_pool <- function(cohort) {
  s <- SEL_COHORTS[[cohort]]
  if (is.null(s)) stop("unknown cohort: ", cohort, call. = FALSE)
  w <- paste0(s$flags, " = 1")
  if (!is.null(s[["include"]])) w <- c(w, paste0("person_id IN (", sel_sql_events(s[["include"]]), ")"))   # [[ ]]: exact names, `$` would match exclude_sex
  if (!is.null(s[["exclude"]])) w <- c(w, paste0("person_id NOT IN (", sel_sql_events(s[["exclude"]]), ")"))
  if (!is.null(s[["exclude_sex"]]))
    w <- c(w, paste0("person_id NOT IN (SELECT person_id FROM `person` WHERE sex_at_birth_concept_id IN (",
                     sel_num(s[["exclude_sex"]]), "))"))
  paste0("SELECT DISTINCT person_id FROM `cb_search_person` WHERE ", paste(w, collapse = "\n      AND "))
}

sel_sql_person <- function(cohort) paste0(
  "SELECT person_id, birth_datetime AS date_of_birth FROM `person` WHERE person_id IN (", sel_sql_pool(cohort), ")")

# Per person: n_days (distinct dates), first_date, last_date, max_gap (longest run of missing dates
# between first and last date; 0 if none) - the definition of filter_pid_func.
sel_sql_wear <- function(cohort, stream) {
  st <- SEL_STREAMS[[stream]]
  paste0("
    WITH d AS (
      SELECT DISTINCT person_id, CAST(", st$date_col, " AS DATE) AS dt FROM `", st$table, "`
      WHERE ", st$date_col, " IS NOT NULL AND person_id IN (", sel_sql_pool(cohort), ")),
    g AS (
      SELECT person_id, dt, DATE_DIFF(dt, LAG(dt) OVER (PARTITION BY person_id ORDER BY dt), DAY) - 1 AS gap FROM d)
    SELECT person_id, COUNT(*) AS n_days, MIN(dt) AS first_date, MAX(dt) AS last_date, IFNULL(MAX(gap), 0) AS max_gap
    FROM g GROUP BY person_id")
}

# Distinct dates of given persons - only to verify sel_sql_wear on a sample.
sel_sql_dates <- function(stream, person_ids) {
  st <- SEL_STREAMS[[stream]]
  paste0("SELECT DISTINCT person_id, CAST(", st$date_col, " AS DATE) AS date FROM `", st$table, "`
    WHERE ", st$date_col, " IS NOT NULL AND person_id IN (", sel_num(person_ids), ")")
}

sel_sql_conditions <- function(cohort, codes) paste0("
    SELECT DISTINCT co.person_id, c.concept_code AS standard_concept_code
    FROM `condition_occurrence` co JOIN `concept` c ON co.condition_concept_id = c.concept_id
    WHERE c.concept_code IN (", paste0("'", codes, "'", collapse = ", "), ")
      AND co.person_id IN (", sel_sql_pool(cohort), ")")

# Run a query in the CDR and return a data.table (BIGINT as integer).
sel_fetch <- function(sql, cdr = Sys.getenv("WORKSPACE_CDR"), billing = Sys.getenv("GOOGLE_PROJECT")) {
  if (!grepl("^[^.[:space:]]+\\.[^.[:space:]]+$", cdr))
    stop("WORKSPACE_CDR must look like project.dataset but is '", cdr, "': run Cell 1 in this R session", call. = FALSE)
  if (!nzchar(billing)) stop("GOOGLE_PROJECT is empty", call. = FALSE)
  data.table::as.data.table(
    bigrquery::bq_table_download(bigrquery::bq_dataset_query(cdr, sql, billing = billing), bigint = "integer"))
}
