# ==============================================================================
# summary_tables.R
# ------------------------------------------------------------------------------
# Table 1 (demographics), the Fitbit characteristics table and the cohort-overlap Venn diagram
# of a selection, from its own outputs: final_pid_<tag>.txt, persons.tsv (age at first record),
# cci.txt, wear_metrics.tsv, plus one query of the person table (race, ethnicity, sex).
# Needs format_med_iqr() (missingness/R/utils.R), race_group_from_name() (missingness/R/persons.R),
# sel_fetch() / sel_num() (selection_sql.R).
#   variant "all"          : every cohort membership, N = rows of final_pid (as in the earlier tables)
#   variant "excl_overlap" : people in more than one cohort removed (the pipeline's "exclude" policy)
# ==============================================================================

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

ST_COHORTS <- c("Insomnia", "Diabetes", "Prostate", "Healthy")
ST_LABELS  <- c(Insomnia = "Insomnia", Diabetes = "Diabetes", Prostate = "Prostate cancer", Healthy = "Healthy")
ST_COLORS  <- c(Diabetes = "#E74C3C", Prostate = "#3498DB", Insomnia = "#F39C12", Healthy = "#27AE60")

st_sql_meta <- function(person_ids) paste0("
    SELECT p.person_id, p.ethnicity_concept_id, r.concept_name AS race, s.concept_name AS sex_at_birth
    FROM `person` p
    LEFT JOIN `concept` r ON p.race_concept_id = r.concept_id
    LEFT JOIN `concept` s ON p.sex_at_birth_concept_id = s.concept_id
    WHERE p.person_id IN (", sel_num(person_ids), ")")

# "n (x%)": "-" for 0, "<1%" below one percent, thousands separator.
st_n_pct <- function(n, N) {
  p <- round(100 * n / N)
  ifelse(n == 0, "-", paste0(format(n, big.mark = ","), ifelse(p < 1, " (<1%)", paste0(" (", p, "%)"))))
}

# People in exactly one cohort.
st_exclude_overlap <- function(final_pid) {
  final_pid[person_id %in% final_pid[, .N, by = person_id][N == 1L, person_id]]
}

# One row per membership with the variables of Table 1.
# final_pid: cohort, person_id | persons: cohort, person_id, age_first | cci: person_id, cci_score
# meta: person_id, ethnicity_concept_id, race, sex_at_birth
st_members <- function(final_pid, persons, cci, meta) {
  m <- merge(final_pid[, .(cohort, person_id)], persons[, .(cohort, person_id, age_first)],
             by = c("cohort", "person_id"), all.x = TRUE)
  m <- merge(m, unique(cci[, .(person_id, cci_score)], by = "person_id"), by = "person_id", all.x = TRUE)
  m <- merge(m, unique(meta[, .(person_id, ethnicity_concept_id, race, sex_at_birth)], by = "person_id"),
             by = "person_id", all.x = TRUE)
  m[, age_group := cut(age_first, c(18, 45, 65, Inf), c("18-44", "45-64", ">=65"), right = FALSE)]
  m[, ethnicity := fcase(ethnicity_concept_id == 38003563, "Hispanic or Latino",
                         ethnicity_concept_id == 38003564, "Not Hispanic or Latino", default = "Others")]
  m[, race_group := race_group_from_name(race)]
  m[, sex := fcase(sex_at_birth == "Female", "Female", sex_at_birth == "Male", "Male", default = "Others")]
  m[, cci_group := cut(cci_score, c(-Inf, 0, 2, 4, Inf), c("0", "1-2", "3-4", ">=5"))]
  m[]
}

# Table 1: one row per category, one column per cohort.
st_table1 <- function(m) {
  blocks <- list(
    list("Age",          "age_group",  c("18-44", "45-64", ">=65")),
    list("Ethnicity",    "ethnicity",  c("Hispanic or Latino", "Not Hispanic or Latino", "Others")),
    list("Race",         "race_group", c("American Indian or Alaska Native", "Asian", "Black or African American",
                                         "Middle Eastern or North African", "Native Hawaiian or Other Pacific Islander",
                                         "White", "More than one population", "Others")),
    list("Sex at birth", "sex",        c("Female", "Male", "Others")),
    list("Charlson comorbidity index (CCI)", "cci_group", c("0", "1-2", "3-4", ">=5")))
  N <- m[, .N, by = cohort]
  col <- function(f) vapply(ST_COHORTS, function(co) f(m[cohort == co], N[cohort == co, N]), character(1))
  rows <- list(c(Demographic = "N", col(function(d, n) format(n, big.mark = ","))))
  for (b in blocks) {
    rows[[length(rows) + 1L]] <- c(Demographic = b[[1]], setNames(rep("", length(ST_COHORTS)), ST_COHORTS))
    for (lv in b[[3]])
      rows[[length(rows) + 1L]] <- c(Demographic = paste0("  ", lv),
                                     col(function(d, n) st_n_pct(sum(d[[b[[2]]]] == lv, na.rm = TRUE), n)))
  }
  out <- as.data.table(do.call(rbind, rows))
  setnames(out, ST_COHORTS, unname(ST_LABELS[ST_COHORTS]))
  out[]
}

# Fitbit characteristics: median (IQR) per stream and cohort.
# wear: cohort, person_id, stream, n_days, followup_days, pct_wear, max_gap
st_fitbit_table <- function(wear, final_pid) {
  w <- merge(wear, final_pid[, .(cohort, person_id)], by = c("cohort", "person_id"))
  w[, missing_pct := 100 - pct_wear]
  metrics <- list(c("Median follow-up days (IQR)", "followup_days"), c("Median record days (IQR)", "n_days"),
                  c("Median missing rate (%, IQR)", "missing_pct"),
                  c("Median max continuous non-wear gap days (IQR)", "max_gap"))
  streams <- c(activity = "Activity", sleep = "Sleep", hr = "Heart rate")
  rows <- list()
  for (mt in metrics) for (s in names(streams))
    rows[[length(rows) + 1L]] <- c(`Fitbit characteristics` = mt[1], Metrics = streams[[s]],
      vapply(ST_COHORTS, function(co) format_med_iqr(w[cohort == co & stream == s][[mt[2]]]), character(1)))
  out <- as.data.table(do.call(rbind, rows))
  setnames(out, ST_COHORTS, unname(ST_LABELS[ST_COHORTS]))
  out[]
}

# Counts of the 7 regions of the three disease cohorts (Healthy never overlaps them).
st_venn_counts <- function(final_pid) {
  sets <- lapply(c("Diabetes", "Insomnia", "Prostate"), function(co) final_pid[cohort == co, person_id])
  names(sets) <- c("Diabetes", "Insomnia", "Prostate")
  ids <- unique(unlist(sets))
  key <- vapply(ids, function(i) paste(names(sets)[vapply(sets, function(s) i %in% s, logical(1))], collapse = " + "), "")
  data.table(region = names(table(key)), n_people = as.integer(table(key)))[order(-n_people)]
}

# Three-circle Venn (circles of equal size; regions are NOT proportional to the counts). Healthy is
# shown as a separate box because no one is in Healthy and another cohort.
st_venn_plot <- function(final_pid) {
  cnt <- st_venn_counts(final_pid)
  n <- function(r) { v <- cnt[region == r, n_people]; if (length(v)) v else 0L }
  ct <- list(Diabetes = c(0, 0), Insomnia = c(1, 0), Prostate = c(0.5, 0.866))
  circ <- rbindlist(lapply(names(ct), function(k) {
    a <- seq(0, 2 * pi, length.out = 200)
    data.table(cohort = k, x = ct[[k]][1] + cos(a), y = ct[[k]][2] + sin(a))
  }))
  tot <- final_pid[, .N, by = cohort]
  lab_names <- data.table(cohort = c("Diabetes", "Insomnia", "Prostate"),
                          x = c(-0.6, 1.6, 0.5), y = c(-1.3, -1.3, 2.0),
                          label = vapply(c("Diabetes", "Insomnia", "Prostate"), function(k)
                            sprintf("%s\nN = %s", ST_LABELS[[k]], format(tot[cohort == k, N], big.mark = ",")), ""))
  reg <- data.table(
    x = c(-0.3, 1.3, 0.5, 0.5, 0.0, 1.0, 0.5),
    y = c(-0.2, -0.2, 1.45, -0.4, 0.75, 0.75, 0.3),
    n = c(n("Diabetes"), n("Insomnia"), n("Prostate"), n("Diabetes + Insomnia"),
          n("Diabetes + Prostate"), n("Insomnia + Prostate"), n("Diabetes + Insomnia + Prostate")))
  nh <- tot[cohort == "Healthy", N]
  ggplot() +
    geom_polygon(data = circ, aes(x, y, fill = cohort, group = cohort), alpha = 0.35, colour = NA) +
    geom_path(data = circ, aes(x, y, colour = cohort, group = cohort), linewidth = 0.9) +
    geom_text(data = reg, aes(x, y, label = format(n, big.mark = ",")), size = 4.2, fontface = "bold") +
    geom_text(data = lab_names, aes(x, y, label = label, colour = cohort), size = 3.8, fontface = "bold", lineheight = 0.95) +
    annotate("label", x = 2.95, y = 0.45, size = 3.6, colour = ST_COLORS[["Healthy"]], fontface = "bold",
             label = sprintf("%s\nN = %s\nno overlap with the\nother cohorts", ST_LABELS[["Healthy"]], format(nh, big.mark = ","))) +
    scale_fill_manual(values = ST_COLORS, guide = "none") + scale_colour_manual(values = ST_COLORS, guide = "none") +
    coord_fixed(xlim = c(-1.3, 3.7), ylim = c(-1.7, 2.3)) + theme_void() +
    labs(title = "Overlap of the cohorts (people with more than one cohort)")
}

# Everything for one selection (tag), written to <out_dir>/tables and <out_dir>/figures.
# meta: the result of sel_fetch(st_sql_meta(...)); passed in so the function itself needs no BigQuery.
st_run <- function(out_dir, tag, meta) {
  dir.create(file.path(out_dir, "tables"), showWarnings = FALSE); dir.create(file.path(out_dir, "figures"), showWarnings = FALSE)
  final_pid <- fread(file.path(out_dir, sprintf("final_pid_%s.txt", tag)))
  persons   <- fread(file.path(out_dir, "persons.tsv"))
  cci       <- fread(file.path(out_dir, "cci.txt"))
  wear      <- fread(file.path(out_dir, "wear_metrics.tsv"))
  res <- list()
  for (v in c("all", "excl_overlap")) {
    fp <- if (v == "all") final_pid else st_exclude_overlap(final_pid)
    t1 <- st_table1(st_members(fp, persons, cci, meta)); t2 <- st_fitbit_table(wear, fp)
    fwrite(t1, file.path(out_dir, "tables", sprintf("table1_demographics_%s_%s.tsv", tag, v)), sep = "\t")
    fwrite(t2, file.path(out_dir, "tables", sprintf("table_fitbit_%s_%s.tsv", tag, v)), sep = "\t")
    res[[v]] <- list(table1 = t1, fitbit = t2)
  }
  vc <- st_venn_counts(final_pid)
  fwrite(vc, file.path(out_dir, "tables", sprintf("venn_counts_%s.tsv", tag)), sep = "\t")
  p <- st_venn_plot(final_pid)
  ggsave(file.path(out_dir, "figures", sprintf("venn_%s.png", tag)), p, width = 8, height = 5.5, dpi = 200, bg = "white")
  ggsave(file.path(out_dir, "figures", sprintf("venn_%s.pdf", tag)), p, width = 8, height = 5.5)
  res$venn <- vc
  invisible(res)
}
