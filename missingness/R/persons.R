# ==============================================================================
# persons.R
# ------------------------------------------------------------------------------
# One row per person: cohort, sex, age, race group and CCI group. Used both as
# adjustment covariates (age, sex, cohort) and as strata for the descriptive plots.
#
# Differences from the earlier code (deliberate, documented):
#   * age is measured at the person's FIRST Fitbit record, not at Sys.Date(),
#     so results do not change with the day the code is run;
#   * race is grouped from the race NAME (e.g. "Native Hawaiian or Other Pacific
#     Islander") instead of a concept-id map that used the key "8657" twice.
# ==============================================================================

# Group the race name into the categories of Table 1.
race_group_from_name <- function(x) {
  x0 <- tolower(trimws(as.character(x)))
  fcase(
    grepl("american indian|alaska native", x0), "American Indian or Alaska Native",
    x0 == "asian",                              "Asian",
    grepl("black",          x0),                "Black or African American",
    grepl("middle eastern", x0),                "Middle Eastern or North African",
    grepl("native hawaiian|pacific islander", x0), "Native Hawaiian or Other Pacific Islander",
    x0 == "white",                              "White",
    grepl("more than one",  x0),                "More than one population",
    default = "Others")
}

# For plots, small groups are merged into "Other" (each is < 2% of any cohort).
race_group_for_plots <- function(race_group) {
  keep <- c("White", "Black or African American", "Asian", "More than one population")
  fifelse(race_group %in% keep, race_group, "Other")
}

# Some people belong to more than one cohort (e.g. insomnia and diabetes). The cohort
# comparisons assume independent groups, so each person must count once. `policy` decides how:
#   "exclude"  (default) people in more than one cohort are left out of the analysis
#   "priority" they are kept in the first of their cohorts listed in `priority`
#   "error"    stop, so the overlap is looked at before anything is run
# Returns the one-row-per-person assignment, the number of people per overlap combination and
# the ids that were affected. The exports are NOT affected: they were made per cohort.
resolve_cohort_overlap <- function(final_pid, policy = "exclude",
                                   priority = c("Prostate", "Diabetes", "Insomnia", "Healthy")) {
  m <- unique(data.table::as.data.table(final_pid)[, .(person_id, cohort = as.character(cohort))])
  per <- m[, .(n_cohorts = .N, combination = paste(sort(cohort), collapse = " + ")), by = person_id]
  overlap <- per[n_cohorts > 1L, .(n_people = .N), by = combination][order(-n_people)]
  multi_ids <- per[n_cohorts > 1L, person_id]
  assignment <- m[!person_id %in% multi_ids]

  if (length(multi_ids) > 0L) {
    if (policy == "error") {
      stop(length(multi_ids), " people belong to more than one cohort (", paste0(overlap$combination, ": ",
           overlap$n_people, collapse = "; "), "). Set cfg$overlap$policy to \"exclude\" or \"priority\".",
           call. = FALSE)
    } else if (policy == "priority") {
      missing_c <- setdiff(unique(m$cohort), priority)
      if (length(missing_c) > 0L) stop("priority list lacks cohort(s): ", paste(missing_c, collapse = ", "), call. = FALSE)
      kept <- m[person_id %in% multi_ids][order(match(cohort, priority))][, .SD[1L], by = person_id]
      assignment <- rbind(assignment, kept)
    } else if (policy != "exclude") stop("unknown overlap policy: ", policy, call. = FALSE)
  }
  list(assignment = assignment[order(person_id)], overlap = overlap,
       overlap_ids = multi_ids, policy = policy,
       n_unique_people = nrow(per), n_kept = nrow(assignment))
}

# meta        : person_id, date_of_birth, race, sex_at_birth   (all cohorts stacked)
# final_pid   : person_id, cohort   ONE row per person (use resolve_cohort_overlap()$assignment)
# cci         : person_id, cci_score
# first_dates : person_id, first_date  (earliest record over the three streams)
build_person_table <- function(meta, final_pid, cci, first_dates, cfg) {
  meta <- data.table::as.data.table(meta)
  assert_columns(meta, c("person_id", "date_of_birth", "race", "sex_at_birth"), "meta")
  assert_unique_keys(meta, "person_id", "meta")
  assert_unique_keys(final_pid, "person_id", "final_pid (resolve cohort overlap first)")
  assert_columns(first_dates, c("person_id", "first_date"), "first_dates")

  p <- merge(final_pid[, .(person_id, cohort)], meta, by = "person_id", all.x = TRUE)
  p <- merge(p, first_dates, by = "person_id", all.x = TRUE)
  p <- merge(p, cci[, .(person_id, cci_score)], by = "person_id", all.x = TRUE)

  n_no_meta <- p[is.na(date_of_birth), .N]
  n_no_cci  <- p[is.na(cci_score), .N]
  if (n_no_meta > 0L) log_info("persons: %d person(s) without meta data", n_no_meta)
  if (n_no_cci  > 0L) log_info("persons: %d person(s) without a CCI score", n_no_cci)

  p[, dob := as.Date(substr(as.character(date_of_birth), 1, 10))]
  p[, age_first := age_on(dob, first_date)]
  p[, age_group := cut(age_first, cfg$strata$age_breaks, cfg$strata$age_labels)]
  p[, sex := factor(fcase(sex_at_birth == "Female", "Female",
                          sex_at_birth == "Male",   "Male",
                          default = "Other"),
                    levels = c("Female", "Male", "Other"))]
  p[, race_group := race_group_from_name(race)]
  p[, race_plot  := factor(race_group_for_plots(race_group),
                           levels = names(cfg$palette$race))]
  p[, cci_group := cut(cci_score, cfg$strata$cci_breaks, cfg$strata$cci_labels)]
  p[, cohort := factor(cohort, levels = cfg$cohorts)]
  p[, c("date_of_birth", "dob", "race", "sex_at_birth") := NULL]
  setkey(p, person_id)
  p[]
}

# Earliest record date per person over a list of day tables.
first_record_dates <- function(day_tables) {
  rbindlist(lapply(day_tables, function(d) d[, .(first_date = min(date)), by = person_id]))[
    , .(first_date = min(first_date)), by = person_id]
}
