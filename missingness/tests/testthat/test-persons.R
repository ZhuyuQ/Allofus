test_that("age_on handles birthdays", {
  expect_equal(age_on(as.Date("1960-06-15"), as.Date("2024-06-14")), 63L)
  expect_equal(age_on(as.Date("1960-06-15"), as.Date("2024-06-15")), 64L)
})

test_that("race is grouped by name (NHPI is its own group) and plots merge small groups", {
  x <- c("White", "Native Hawaiian or Other Pacific Islander", "None of these", "Black or African American",
         "American Indian or Alaska Native", "More than one population")
  g <- race_group_from_name(x)
  expect_equal(g[2], "Native Hawaiian or Other Pacific Islander"); expect_equal(g[3], "Others")
  expect_equal(race_group_for_plots(g), c("White", "Other", "Other", "Black or African American", "Other", "More than one population"))
})

test_that("person table builds age at first record, groups and CCI bins", {
  cfg <- cfg_test()
  meta <- data.table(person_id = 1:3, date_of_birth = c("1990-01-01 00:00:00 UTC", "1970-01-01 00:00:00 UTC", "1940-01-01 00:00:00 UTC"),
                     race = c("White", "Asian", "White"), sex_at_birth = c("Female", "Male", "PMI: Skip"))
  fp <- data.table(person_id = 1:3, cohort = c("Healthy", "Diabetes", "Prostate"))
  cci <- data.table(person_id = 1:3, cci_score = c(0, 2, 7))
  fd <- data.table(person_id = 1:3, first_date = as.Date("2022-06-01"))
  p <- build_person_table(meta, fp, cci, fd, cfg)
  expect_equal(p$age_first, c(32L, 52L, 82L))
  expect_equal(as.character(p$age_group), c("18-44", "45-64", ">=65"))
  expect_equal(as.character(p$cci_group), c("0", "1-2", ">=5")); expect_equal(as.character(p$sex), c("Female", "Male", "Other"))
})

test_that("overlap policies: exclude, priority, error; no overlap is a no-op", {
  fp <- data.table(person_id = c(1, 2, 3, 3, 4, 4, 4),
                   cohort = c("Healthy", "Insomnia", "Insomnia", "Diabetes", "Prostate", "Diabetes", "Insomnia"))
  ex <- resolve_cohort_overlap(fp, "exclude")
  expect_equal(ex$assignment$person_id, c(1, 2)); expect_equal(ex$n_unique_people, 4L)
  expect_setequal(ex$overlap_ids, c(3, 4)); expect_equal(ex$overlap[combination == "Diabetes + Insomnia", n_people], 1L)
  pr <- resolve_cohort_overlap(fp, "priority")
  expect_equal(pr$assignment[person_id == 3, cohort], "Diabetes"); expect_equal(pr$assignment[person_id == 4, cohort], "Prostate")
  expect_equal(nrow(pr$assignment), 4L)
  expect_error(resolve_cohort_overlap(fp, "error"), "belong to more than one cohort")
  expect_error(resolve_cohort_overlap(fp, "priority", "Healthy"), "lacks cohort")
  expect_equal(nrow(resolve_cohort_overlap(fp[1:2])$overlap), 0L)
  expect_error(build_person_table(data.table(person_id = 1, date_of_birth = "1990-01-01", race = "White", sex_at_birth = "Male"),
                                  fp, data.table(person_id = 1, cci_score = 0), data.table(person_id = 1, first_date = Sys.Date()), make_config()),
               "resolve cohort overlap first")
})
