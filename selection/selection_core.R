# ==============================================================================
# selection_core.R
# ------------------------------------------------------------------------------
# Selection rules in plain R (no BigQuery). Needs age_on() from missingness/R/utils.R and
# cci_calc_func() from function/cci_function.R, both sourced unchanged.
#
#   rules   : per stream (activity, sleep, hr) days >= min_days, longest gap <= gap_max, compliance
#             in [comp_min, comp_max) percent of the first-to-last span; a person must pass in ALL streams;
#             then CCI <= cci_max. Defaults = the earlier rule (21 days, 14 days, 70%, CCI 10).
#   age     : whole years at the person's earliest record over the three streams.
#   compliance is compared exactly (n_days * 100 >= comp_min * followup_days), not through a rounded ratio.
# ==============================================================================

suppressPackageStartupMessages(library(data.table))

# R twin of sel_sql_wear(): used to verify the SQL on a sample and in the tests.
# dt: person_id, date.
wear_from_dates <- function(dt) {
  d <- unique(as.data.table(dt)[, .(person_id, date = as.Date(date))])
  setorder(d, person_id, date)
  d[, gap := c(NA_integer_, as.integer(diff(date)) - 1L), by = person_id]
  d[, .(n_days = .N, first_date = date[1L], last_date = date[.N], max_gap = max(c(0L, gap), na.rm = TRUE)),
    by = person_id]
}

# Adds followup_days (first to last date, inclusive) and pct_wear to the SQL result.
sel_prepare_wear <- function(w) {
  w <- as.data.table(w)
  w[, `:=`(first_date = as.Date(first_date), last_date = as.Date(last_date))]
  w[, followup_days := as.integer(last_date - first_date) + 1L]
  w[, pct_wear := 100 * n_days / followup_days]
  stopifnot("n_days exceeds the follow-up span" = all(w$n_days <= w$followup_days))
  w[]
}

# Person table with the age at the first record over the three streams.
# persons: cohort, person_id, dob ; wear: cohort, person_id, first_date (any stream).
sel_persons <- function(persons, wear) {
  first <- wear[, .(first_date = min(first_date)), by = .(cohort, person_id)]
  p <- merge(as.data.table(persons), first, by = c("cohort", "person_id"), all.x = TRUE)
  p[, age_first := age_on(dob, first_date)]
  p[]
}

# The SNOMED codes of the CCI lookup, read from cci_calc_func() itself so the two never differ.
sel_cci_codes <- function() {
  hit <- Filter(function(e) is.call(e) && identical(e[[1L]], as.name("<-")) &&
                  identical(e[[2L]], as.name("cci_lookup")), as.list(body(cci_calc_func)))
  stopifnot("cci_lookup not found in cci_calc_func" = length(hit) == 1L)
  unique(as.character(eval(hit[[1L]][[3L]])$standard_concept_code))
}

# CCI with cci_calc_func() and the age at the first record. persons: cohort, person_id, age_first;
# conds: cohort, person_id, standard_concept_code. People without condition rows get age points only.
sel_cci <- function(persons, conds) {
  p <- as.data.table(persons)[!is.na(age_first), .(cohort, person_id, age = age_first)]
  has <- merge(unique(as.data.table(conds)[, .(cohort, person_id, standard_concept_code)]), p,
               by = c("cohort", "person_id"))
  df <- rbind(has, cbind(p, standard_concept_code = NA_character_), use.names = TRUE)
  as.data.table(cci_calc_func(df = as.data.frame(df), path.out = NULL))
}

# Apply the rules. wear: sel_prepare_wear() output with cohort, stream; cci: cohort, person_id, cci_score;
# pool: cohort, person_id of everyone in the cohort. Returns the selected people and the attrition.
sel_select <- function(wear, cci, pool, min_days = 21L, gap_max = 14L, comp_min = 70, comp_max = Inf,
                       cci_max = 10, streams = c("activity", "sleep", "hr")) {
  per <- wear[stream %in% streams,
              .(n_streams = .N,
                ok_days = all(n_days >= min_days),
                ok_gap  = all(max_gap <= gap_max),
                ok_comp = all(n_days * 100 >= comp_min * followup_days & n_days * 100 < comp_max * followup_days)),
              by = .(cohort, person_id)]
  a <- merge(unique(pool[, .(cohort, person_id)]), per, by = c("cohort", "person_id"), all.x = TRUE)
  a <- merge(a, unique(cci[, .(cohort, person_id, cci_score)]), by = c("cohort", "person_id"), all.x = TRUE)
  s1 <- a$n_streams %in% length(streams)
  s2 <- s1 & a$ok_days %in% TRUE
  s3 <- s2 & a$ok_gap %in% TRUE
  s4 <- s3 & a$ok_comp %in% TRUE
  s5 <- s4 & !is.na(a$cci_score) & a$cci_score <= cci_max
  lab <- c(s0 = "pool", s1 = paste0("data in ", paste(streams, collapse = "+")),
           s2 = sprintf("days >= %d", min_days), s3 = sprintf("longest gap <= %d", gap_max),
           s4 = sprintf("compliance in [%g, %g)", comp_min, comp_max), s5 = sprintf("CCI <= %g", cci_max))
  att <- melt(data.table(cohort = a$cohort, s0 = TRUE, s1, s2, s3, s4, s5)[, lapply(.SD, sum), by = cohort],
              id.vars = "cohort", variable.name = "step", value.name = "n")
  att[, step := factor(lab[as.character(step)], levels = lab)]
  list(pids = a[s5, .(cohort, person_id)][order(cohort, person_id)], attrition = att[order(cohort, step)])
}

# Selected people of two lists, per cohort.
sel_compare <- function(new, old) {
  rbindlist(lapply(union(new$cohort, old$cohort), function(co) {
    n <- new[cohort == co]$person_id; o <- old[cohort == co]$person_id
    data.table(cohort = co, n_new = length(n), n_old = length(o), n_both = length(intersect(n, o)),
               only_new = length(setdiff(n, o)), only_old = length(setdiff(o, n)))
  }))
}

# final_pid_<tag>.txt (cohort, person_id: the format of final_pid.txt) and attrition_<tag>.tsv.
sel_write <- function(sel, out_dir, tag) {
  fwrite(sel$pids, file.path(out_dir, sprintf("final_pid_%s.txt", tag)), sep = "\t", quote = FALSE)
  fwrite(sel$attrition, file.path(out_dir, sprintf("attrition_%s.tsv", tag)), sep = "\t")
  invisible(sel)
}

# Recompute the wear metrics of a random sample of people from their raw dates and compare them with
# the SQL result. fetch(stream, person_ids) returns person_id, date. Stops on any difference.
sel_verify <- function(wear, fetch, n = 30L, seed = 1L) {
  set.seed(seed)
  res <- rbindlist(lapply(unique(wear$stream), function(s) {
    ws <- unique(wear[stream == s, .(person_id, n_days, first_date, last_date, max_gap)])
    ws <- ws[sample.int(.N, min(n, .N))]
    raw <- wear_from_dates(fetch(s, ws$person_id))
    m <- merge(ws, raw, by = "person_id", all.x = TRUE, suffixes = c("", "_raw"))
    m[, .(stream = s, n_checked = .N,
          n_mismatch = sum(is.na(n_days_raw) | n_days != n_days_raw | first_date != first_date_raw |
                             last_date != last_date_raw | max_gap != max_gap_raw))]
  }))
  print(res)
  stopifnot("SQL wear metrics differ from the raw dates" = all(res$n_mismatch == 0L))
  invisible(res)
}
