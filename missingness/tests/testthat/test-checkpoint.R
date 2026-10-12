cfg_ck <- function() cfg_test()
unit_a <- list(kind = "onset", stream = "activity", id = "spec A")
key_a  <- function(cfg = cfg_ck(), data = list(days = list(md5 = "x")), spec = list(k = 1L), sample = NULL)
  ckpt_key(cfg, "onset", data, spec, sample)

test_that("run_unit computes once, then resumes, and never reuses a result made with other settings, data or spec", {
  dir <- tempfile(); n <- 0
  compute <- function() { n <<- n + 1; data.table(v = n) }
  quiet <- function(expr) suppressMessages(expr)

  r1 <- quiet(run_unit(dir, unit_a, key_a(), compute)); expect_equal(r1$status, "computed")
  r2 <- quiet(run_unit(dir, unit_a, key_a(), compute)); expect_equal(r2$status, "resumed")
  expect_equal(n, 1); expect_equal(r2$value$v, 1)

  cfg2 <- cfg_ck(); cfg2$onset$days_since_cap <- 30L                       # result-changing setting
  expect_equal(quiet(run_unit(dir, unit_a, key_a(cfg2), compute))$status, "computed")
  expect_equal(quiet(run_unit(dir, unit_a, key_a(data = list(days = list(md5 = "y"))), compute))$status, "computed")   # new data
  expect_equal(quiet(run_unit(dir, unit_a, key_a(spec = list(k = 2L)), compute))$status, "computed")                     # new spec
  expect_equal(quiet(run_unit(dir, unit_a, key_a(sample = 100L), compute))$status, "computed")                         # new sample
  expect_equal(n, 5)
})

test_that("resource settings never invalidate a checkpoint, result settings always do", {
  cfg <- cfg_ck(); cfg_r <- cfg
  cfg_r$resources$memory_budget_gb <- 3; cfg_r$resources$isolate_fits <- FALSE; cfg_r$resources$light_workers <- 1L
  expect_identical(ckpt_cfg_part(cfg, "onset"), ckpt_cfg_part(cfg_r, "onset"))
  cfg_m <- cfg; cfg_m$model$nAGQ <- 1L
  expect_false(identical(ckpt_cfg_part(cfg, "onset"), ckpt_cfg_part(cfg_m, "onset")))
  cfg_b <- cfg; cfg_b$gap_models$min_bin_n <- 5L                            # gap-model settings do not touch onset results
  expect_identical(ckpt_cfg_part(cfg, "onset"), ckpt_cfg_part(cfg_b, "onset"))
  expect_false(identical(ckpt_cfg_part(cfg, "gapC"), ckpt_cfg_part(cfg_b, "gapC")))
})

test_that("writes are atomic: no temporary file is left, and an unreadable file is recomputed", {
  dir <- tempfile(); n <- 0; compute <- function() { n <<- n + 1; n }
  suppressMessages(run_unit(dir, unit_a, key_a(), compute))
  expect_length(list.files(dir, pattern = "tmp"), 0L)
  f <- list.files(dir, full.names = TRUE)[1]
  writeBin(as.raw(1:20), f)                                                 # a half-written / damaged file
  r <- suppressMessages(run_unit(dir, unit_a, key_a(), compute))
  expect_equal(r$status, "computed"); expect_equal(n, 2)
  expect_equal(suppressMessages(run_unit(dir, unit_a, key_a(), compute))$status, "resumed")
})

test_that("a failed unit is recorded with its error, is retried, and never touches other units", {
  dir <- tempfile(); other <- list(kind = "onset", stream = "activity", id = "spec B")
  suppressMessages(run_unit(dir, other, key_a(), function() data.table(v = 99)))
  boom <- suppressMessages(run_unit(dir, unit_a, key_a(), function() stop("model exploded")))
  expect_equal(boom$status, "failed"); expect_match(boom$error, "model exploded"); expect_null(boom$value)
  expect_equal(suppressMessages(run_unit(dir, other, key_a(), function() stop("must not run")))$value$v, 99)   # still valid

  cu <- collect_units(dir, "onset", cfg_ck())
  expect_equal(nrow(cu$units), 1L); expect_equal(cu$failures$id, "spec A"); expect_match(cu$failures$error, "exploded")
  rep <- checkpoint_report(dir, cfg_ck())
  expect_setequal(rep$state, c("valid", "failed"))

  ok <- suppressMessages(run_unit(dir, unit_a, key_a(), function() data.table(v = 1)))           # retried -> success
  expect_equal(ok$status, "computed")
  expect_equal(nrow(collect_units(dir, "onset", cfg_ck())$failures), 0L)                         # failure record removed
  expect_error(suppressMessages(run_unit(dir, unit_a, key_a(cfg = within(cfg_ck(), resources$stop_on_error <- TRUE), sample = 1L),
                                         function() stop("hard stop"), cfg = within(cfg_ck(), resources$stop_on_error <- TRUE))), "hard stop")
})

test_that("a model that reports status 'failed' is not stored as a valid result", {
  dir <- tempfile(); n <- 0
  bad <- function() { n <<- n + 1; list(info = data.table(status = "failed", note = "did not converge at all")) }
  r <- suppressMessages(run_unit(dir, unit_a, key_a(), bad))
  expect_equal(r$status, "failed"); expect_match(r$error, "did not converge")
  suppressMessages(run_unit(dir, unit_a, key_a(), bad)); expect_equal(n, 2)                      # re-run, not resumed
  expect_equal(length(collect_units(dir, "onset", cfg_ck())$values), 0L)
  # legitimate outcomes (not estimable / not converged) ARE results and are reused
  ne <- function() { n <<- n + 1; list(info = data.table(status = "not_estimable", note = "too few events")) }
  suppressMessages(run_unit(dir, unit_a, key_a(), ne)); suppressMessages(run_unit(dir, unit_a, key_a(), ne))
  expect_equal(n, 3)
})

test_that("collect_units ignores results made with other settings or other data", {
  dir <- tempfile(); dk <- list(activity = list(days = list(md5 = "x")))
  suppressMessages(run_unit(dir, unit_a, key_a(), function() data.table(v = 1)))
  expect_equal(length(collect_units(dir, "onset", cfg_ck(), dk)$values), 1L)
  cfg2 <- cfg_ck(); cfg2$seed <- 1L
  expect_equal(length(collect_units(dir, "onset", cfg2, dk)$values), 0L)
  expect_equal(length(collect_units(dir, "onset", cfg_ck(), list(activity = list(days = list(md5 = "CHANGED"))))$values), 0L)
  rep <- checkpoint_report(dir, cfg_ck(), list(activity = list(days = list(md5 = "CHANGED"))))
  expect_equal(rep$state, "stale"); expect_match(rep$detail, "data changed")
  pr <- prune_checkpoints(dir, cfg_ck(), list(activity = list(days = list(md5 = "CHANGED"))), dry_run = TRUE)
  expect_equal(nrow(pr), 1L); expect_length(list.files(dir), 1L)                                  # dry run deletes nothing
  prune_checkpoints(dir, cfg_ck(), list(activity = list(days = list(md5 = "CHANGED"))), dry_run = FALSE)
  expect_length(list.files(dir), 0L)
})

test_that("checkpoint files record what is needed to judge them later", {
  dir <- tempfile(); suppressMessages(run_unit(dir, unit_a, key_a(), function() { warning("careful"); 1 }))
  o <- readRDS(list.files(dir, full.names = TRUE)[1])
  expect_true(all(c("key", "unit", "value", "created", "elapsed_sec", "warnings", "r_version", "packages") %in% names(o)))
  expect_equal(o$warnings, "careful"); expect_true(all(c("cfg", "data", "spec", "sample") %in% names(o$key)))
  expect_match(basename(list.files(dir)[1]), "^onset__activity__spec_A__[0-9a-f]{8}\\.rds$")      # stable, descriptive name
})

test_that("the day store reads from L1 files, selects columns, and fingerprints change with the data", {
  cfg <- cfg_test(storage_format = "rds"); dir <- tempfile()
  inp <- suppressMessages(simulate_inputs(cfg, n_per_cohort = c(Insomnia = 10, Diabetes = 10, Healthy = 20, Prostate = 8), seed = 2))
  for (s in names(inp$day_tables)) write_day_table(inp$day_tables[[s]], dir, s, "rds")
  saveRDS(inp$persons, file.path(dir, "persons.rds"))
  st <- open_day_store(cfg, dir)
  expect_setequal(names(st$get("activity", "steps")), c("person_id", "date", "steps"))
  expect_equal(nrow(st$get("sleep")), nrow(inp$day_tables$sleep))
  expect_error(st$get("activity", "no_such_column"), "no column")
  expect_true("steps" %in% st$columns("activity"))
  k1 <- data_key(st, "sleep"); expect_true(!is.null(k1$act) && is.null(data_key(st, "activity")$act))
  d <- copy(inp$day_tables$activity); d[1, steps := steps + 1]
  write_day_table(d, dir, "activity", "rds")
  expect_false(identical(k1, data_key(open_day_store(cfg, dir), "sleep")))                         # new store sees the change
  expect_error(open_day_store(cfg, tempfile()), "not found")
  mem <- as_day_store(inp$day_tables, inp$persons)
  expect_false(identical(mem$fingerprint("activity"), as_day_store(modifyList(inp$day_tables, list(activity = d)), inp$persons)$fingerprint("activity")))
})

test_that("lazy builds once, on demand, and can be released", {
  n <- 0; z <- lazy(function() { n <<- n + 1; 42 })
  expect_equal(n, 0); expect_false(z$built()); expect_equal(z$get(), 42); expect_equal(z$get(), 42); expect_equal(n, 1)
  z$release(); expect_false(z$built()); z$get(); expect_equal(n, 2)
})

test_that("worker count is capped by free memory and at least one", {
  cfg <- cfg_test(); cfg$resources$memory_budget_gb <- 10; cfg$resources$memory_reserve_gb <- 0
  expect_equal(n_workers(cfg, per_worker_gb = 100, requested = 4L), 1L)
  expect_lte(n_workers(cfg, per_worker_gb = 0.001, requested = 3L), 3L)
})
