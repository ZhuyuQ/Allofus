# Run from selection/:   Rscript tests/test_summary.R      (writes the example figure to tests/venn_example.png)
suppressMessages({ library(testthat); library(data.table) })
root <- normalizePath("..")
# folder names: this repo has missingness/, the Workbench copy final_missing/
miss <- if (dir.exists(file.path(root, "final_missing"))) "final_missing" else "missingness"
source(file.path(root, miss, "R/utils.R")); source(file.path(root, miss, "R/persons.R"))
source(file.path(root, "selection/selection_sql.R")); source(file.path(root, "selection/summary_tables.R"))

# the sizes of the 2026-10-31 selection: D 1306, I 1388, P 253, H 6381; overlaps DI 296, DP 44, IP 43, DIP 3
mk_pid <- function() {
  id <- 1000000L
  take <- function(n) { x <- id + seq_len(n); id <<- id + n; x }
  di <- take(296); dp <- take(44); ip <- take(43); dip <- take(3)
  d <- c(take(1306 - 296 - 44 - 3), di, dp, dip); i <- c(take(1388 - 296 - 43 - 3), di, ip, dip)
  p <- c(take(253 - 44 - 43 - 3), dp, ip, dip); h <- take(6381)
  data.table(cohort = rep(c("Diabetes", "Insomnia", "Prostate", "Healthy"), lengths(list(d, i, p, h))), person_id = c(d, i, p, h))
}
set.seed(3)
pid <- mk_pid()
persons <- copy(pid)[, age_first := sample(18:90, .N, TRUE)]
cci <- unique(copy(pid)[, .(person_id)])[, cci_score := sample(0:10, .N, TRUE)]
wear <- rbindlist(lapply(c("activity", "sleep", "hr"), function(s) copy(pid)[, `:=`(
  stream = s, followup_days = sample(40:600, .N, TRUE))]))
wear[, n_days := pmax(21L, as.integer(followup_days * runif(.N, 0.7, 1)))][, `:=`(pct_wear = 100 * n_days / followup_days, max_gap = sample(0:14, .N, TRUE))]
meta <- unique(pid[, .(person_id)])[, `:=`(
  ethnicity_concept_id = sample(c(38003563, 38003564, 0), .N, TRUE, c(.1, .85, .05)),
  race = sample(c("White", "Asian", "Black or African American", "More than one population", "Native Hawaiian or Other Pacific Islander", "PMI: Skip"), .N, TRUE),
  sex_at_birth = sample(c("Female", "Male", "PMI: Skip"), .N, TRUE, c(.5, .49, .01)))]

test_that("overlap handling and Venn counts", {
  vc <- st_venn_counts(pid)
  expect_equal(vc[region == "Diabetes + Insomnia", n_people], 296L)
  expect_equal(vc[region == "Diabetes + Insomnia + Prostate", n_people], 3L)
  expect_equal(sum(vc$n_people), uniqueN(pid[cohort != "Healthy", person_id]))
  ex <- st_exclude_overlap(pid)
  expect_equal(nrow(ex), nrow(pid) - 296 * 2 - 44 * 2 - 43 * 2 - 3 * 3)
  expect_equal(ex[, .N, by = cohort][order(cohort), N], c(963L, 6381L, 1046L, 163L))
})

test_that("Table 1 and the Fitbit table have the shape of the earlier tables", {
  m <- st_members(pid, persons, cci, meta)
  t1 <- st_table1(m)
  expect_equal(names(t1), c("Demographic", "Insomnia", "Diabetes", "Prostate cancer", "Healthy"))
  expect_equal(t1[Demographic == "N", `Prostate cancer`], "253")
  expect_equal(t1[Demographic == "N", Insomnia], "1,388")
  ages <- t1[Demographic %in% c("  18-44", "  45-64", "  >=65"), Healthy]
  expect_equal(sum(as.integer(gsub(",", "", sub(" .*", "", ages)))), 6381L)       # every person in one age group
  expect_equal(st_n_pct(0L, 10L), "-"); expect_equal(st_n_pct(1L, 1000L), "1 (<1%)"); expect_equal(st_n_pct(2500L, 5000L), "2,500 (50%)")
  t2 <- st_fitbit_table(wear, pid)
  expect_equal(nrow(t2), 12L); expect_equal(t2[1, Metrics], "Activity")
  expect_match(t2[1, Healthy], "^\\d+ \\(\\d+-\\d+\\)$")
})

test_that("st_run writes the tables and the figure", {
  d <- file.path(tempdir(), "sel"); dir.create(d)
  fwrite(pid, file.path(d, "final_pid_base.txt"), sep = "\t"); fwrite(persons, file.path(d, "persons.tsv"), sep = "\t")
  fwrite(cci, file.path(d, "cci.txt"), sep = "\t"); fwrite(wear, file.path(d, "wear_metrics.tsv"), sep = "\t")
  r <- st_run(d, "base", meta)
  expect_true(all(file.exists(file.path(d, c("tables/table1_demographics_base_all.tsv", "tables/table1_demographics_base_excl_overlap.tsv",
                                            "tables/table_fitbit_base_all.tsv", "tables/venn_counts_base.tsv", "figures/venn_base.png")))))
  expect_equal(r$excl_overlap$table1[Demographic == "N", Diabetes], "963")
  file.copy(file.path(d, "figures/venn_base.png"), "tests/venn_example.png", overwrite = TRUE)
  print(r$all$table1); print(r$all$fitbit)
})
cat("summary tests done\n")
