# ==============================================================================
# utils.R
# ------------------------------------------------------------------------------
# Helpers used by every other file: validation, numeric helpers, computing resources (memory
# diagnostics, isolated workers, lazy objects) and the checkpoint store (run_unit and friends).
# Rule of thumb: validation (`assert_*`) is used only where bad input would
# silently give a wrong answer; everything else trusts its inputs.
# ==============================================================================

suppressPackageStartupMessages(library(data.table))

# `a %||% b` : b when a is NULL
`%||%` <- function(a, b) if (is.null(a)) b else a

# Timestamped progress message (goes to stderr, so it never pollutes outputs).
log_info <- function(...) {
  message(format(Sys.time(), "%H:%M:%S"), " | ", sprintf(...))
}

# ---- validation (used sparingly) ---------------------------------------------

# Stop if required columns are missing. `what` names the table in the message.
assert_columns <- function(dt, cols, what = "table") {
  missing_cols <- setdiff(cols, names(dt))
  if (length(missing_cols) > 0L) {
    stop(sprintf("%s is missing required column(s): %s",
                 what, paste(missing_cols, collapse = ", ")), call. = FALSE)
  }
  invisible(TRUE)
}

# Stop if (keys) do not identify rows uniquely.
assert_unique_keys <- function(dt, keys, what = "table") {
  n_dup <- anyDuplicated(dt, by = keys)
  if (n_dup > 0L) {
    stop(sprintf("%s has duplicated (%s) rows; first duplicate at row %d",
                 what, paste(keys, collapse = ", "), n_dup), call. = FALSE)
  }
  invisible(TRUE)
}

# ---- numeric helpers ---------------------------------------------------------

# sum() that returns NA (not 0) when every value is NA.
safe_sum <- function(x) if (all(is.na(x))) NA_real_ else sum(x, na.rm = TRUE)
safe_min <- function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
safe_max <- function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)

# Row-wise sum of several equal-length vectors; NA only if ALL are NA.
row_sum_na <- function(...) {
  m <- cbind(...)
  out <- rowSums(m, na.rm = TRUE)
  out[rowSums(!is.na(m)) == 0L] <- NA_real_
  out
}

# Named transforms used for predictors (see cfg$predictors).
apply_transform <- function(x, transform) {
  switch(transform,
         identity = x,
         log1p    = log1p(pmax(x, 0)),
         stop("unknown transform: ", transform, call. = FALSE))
}

# "median (Q1-Q3)" as in the paper's Table 2, e.g. "112 (97-469)".
format_med_iqr <- function(x, digits = 0) {
  x <- x[!is.na(x)]
  if (length(x) == 0L) return(NA_character_)
  q <- stats::quantile(x, c(0.5, 0.25, 0.75), names = FALSE)
  sprintf("%s (%s-%s)", formatC(q[1], format = "f", digits = digits),
          formatC(q[2], format = "f", digits = digits),
          formatC(q[3], format = "f", digits = digits))
}

# Whole years between date of birth and a reference date (vectorised).
# Subtracts one year when the birthday has not yet occurred in the reference year.
age_on <- function(dob, ref) {
  dob <- as.Date(dob); ref <- as.Date(ref)
  y <- as.integer(format(ref, "%Y")) - as.integer(format(dob, "%Y"))
  had_birthday <- format(ref, "%m%d") >= format(dob, "%m%d")
  y - as.integer(!had_birthday)
}

# Evaluate `expr` with a fixed RNG seed and restore the previous RNG state after.
with_seed <- function(seed, expr) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) assign(".Random.seed", old, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(seed)
  force(expr)
}

# Seeded random sample of person ids (used when cfg$model$sample_persons is set).
sample_person_ids <- function(ids, n, seed) {
  ids <- unique(ids)
  if (is.null(n) || n >= length(ids)) return(ids)
  with_seed(seed, sample(ids, n))
}

# Cluster (person-level) bootstrap of a ratio sum(num)/sum(den).
# `dt` has one row per (group, person) with counts `num`, `den`. Resamples
# persons within each group. Returns a data.table of estimate and 95% interval.
# Cheap because it works on per-person counts, never on day-level rows.
cluster_boot_ratio <- function(dt, group_cols, num = "num", den = "den",
                               id = "person_id", B = 200L, seed = 1L) {
  stopifnot(all(c(group_cols, num, den, id) %in% names(dt)))
  one_group <- function(sub) {
    n <- nrow(sub)
    est <- sum(sub[[num]]) / sum(sub[[den]])
    if (n < 2L || B < 2L) return(list(estimate = est, lo = NA_real_, hi = NA_real_, n_persons = n))
    nv <- sub[[num]]; dv <- sub[[den]]
    draws <- vapply(seq_len(B), function(b) {
      i <- sample.int(n, n, replace = TRUE)
      sum(nv[i]) / sum(dv[i])
    }, numeric(1))
    q <- stats::quantile(draws, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
    list(estimate = est, lo = q[1], hi = q[2], n_persons = n)
  }
  with_seed(seed, dt[, one_group(.SD), by = group_cols])
}

# Attach a package for the duration of `expr` and detach it afterwards (only if it was
# not already attached). survival::clogit() looks up coxph() on the search path.
with_attached <- function(pkg, expr) {
  name <- paste0("package:", pkg)
  if (!name %in% search()) {
    suppressPackageStartupMessages(attachNamespace(pkg))
    on.exit(detach(name, character.only = TRUE), add = TRUE)
  }
  force(expr)
}

# Grouped sum that keeps "all NA stays NA" but stays fast on millions of groups.
# (A per-group R function such as safe_sum() is far slower than data.table's optimised
# sum, so non-missing counts are summed alongside and used to restore the NAs.)
grouped_sum_na <- function(dt, by, cols) {
  nn <- paste0(".nn_", cols)
  tmp <- dt[, c(by, cols), with = FALSE]
  tmp[, (nn) := lapply(.SD, function(x) as.integer(!is.na(x))), .SDcols = cols]
  sums <- tmp[, lapply(.SD, sum, na.rm = TRUE), by = by, .SDcols = cols]
  cnts <- tmp[, lapply(.SD, sum), by = by, .SDcols = nn]
  for (i in seq_along(cols)) set(sums, which(cnts[[nn[i]]] == 0L), cols[i], NA_real_)
  sums[]
}

# ==============================================================================
# Computing resources
# ------------------------------------------------------------------------------
# R does not return memory to the operating system after a large model, so the pipeline
#   (1) keeps one stream's big tables alive at a time (see the lazy bundle in pipeline.R),
#   (2) runs every heavy fit in a short-lived worker process (run_isolated), and
#   (3) checks BEFORE a fit whether it can fit in the free memory (check_fit_memory).
# All limits live in cfg$resources (config.R).
# ==============================================================================

# Memory (GB) held by THIS R process right now (resident set). NA when /proc is unavailable.
rss_gb <- function() {
  st <- if (file.exists("/proc/self/status")) readLines("/proc/self/status") else character()
  l <- st[grepl("^VmRSS:", st)]
  if (length(l)) as.numeric(gsub("[^0-9]", "", l)) / 1024^2 else NA_real_
}

# Memory (GB) that is still free for this session: the smaller of the system's available memory
# and what the container limit (cgroup) still allows. NA when it cannot be read.
available_memory_gb <- function() {
  num <- function(path) if (file.exists(path)) suppressWarnings(as.numeric(readLines(path, n = 1L))) else NA_real_
  out <- NA_real_
  mi <- if (file.exists("/proc/meminfo")) readLines("/proc/meminfo") else character()
  l <- mi[grepl("^MemAvailable:", mi)]
  if (length(l)) out <- as.numeric(gsub("[^0-9]", "", l)) / 1024^2                    # kB -> GB
  lim <- num("/sys/fs/cgroup/memory.max"); use <- num("/sys/fs/cgroup/memory.current")       # cgroup v2
  if (is.na(lim)) { lim <- num("/sys/fs/cgroup/memory/memory.limit_in_bytes")                # cgroup v1
                    use <- num("/sys/fs/cgroup/memory/memory.usage_in_bytes") }
  if (!is.na(lim) && !is.na(use) && lim < 1e15) out <- suppressWarnings(min(out, (lim - use) / 1024^3, na.rm = TRUE))
  if (is.infinite(out)) NA_real_ else out
}

# Memory (GB) one more heavy computation may use: free memory minus the reserve, capped by the budget.
memory_headroom_gb <- function(cfg) {
  avail <- available_memory_gb()
  if (is.na(avail)) return(cfg$resources$memory_budget_gb)
  max(0, min(avail - cfg$resources$memory_reserve_gb, cfg$resources$memory_budget_gb))
}

# After a big worker exits the kernel needs a moment to take its memory back, so free memory looks too
# small for a few seconds. Wait (up to cfg$resources$memory_wait_sec) for `need_gb` of headroom.
wait_for_headroom <- function(cfg, need_gb) {
  h <- memory_headroom_gb(cfg); t0 <- Sys.time()
  while (h < need_gb && as.numeric(difftime(Sys.time(), t0, units = "secs")) < cfg$resources$memory_wait_sec) {
    Sys.sleep(1); h <- memory_headroom_gb(cfg)
  }
  h
}

# Number of workers for a parallel step whose workers each need about `per_worker_gb`.
# Never more than requested and never more than the free memory allows (at least 1).
n_workers <- function(cfg, per_worker_gb, requested) {
  by_memory <- floor(memory_headroom_gb(cfg) / max(per_worker_gb, 0.05))
  max(1L, as.integer(min(requested, by_memory)))
}

# Evaluate `expr` in a short-lived forked worker process and return its value.
# Why: after a large model fit R does not hand the memory back to the operating system, so
# several fits in one session add up until the runtime runs out of memory and the kernel dies.
# A worker's memory is released completely when it exits; only the (small) result is copied back,
# and the input data are shared with the parent (fork), not serialised.
# Falls back to normal evaluation when forking is unavailable (Windows) or `enabled` is FALSE.
run_isolated <- function(expr, enabled = TRUE) {
  if (!enabled || .Platform$OS.type != "unix") return(expr)
  job <- parallel::mcparallel(try(expr, silent = TRUE))
  out <- parallel::mccollect(job)[[1L]]
  if (is.null(out)) {
    stop("a model worker process died - most likely the runtime ran out of memory. Use a runtime with ",
         "more memory, or set cfg$model$sample_persons to fit fewer people.", call. = FALSE)
  }
  if (inherits(out, "try-error")) stop(conditionMessage(attr(out, "condition")), call. = FALSE)
  out
}

# run_isolated() configured by cfg$resources.
isolated <- function(cfg, expr) run_isolated(expr, enabled = isTRUE(cfg$resources$isolate_fits))

# A value computed on first use and dropped on release(): lets the pipeline build a large table
# only when some model really has to run, and free it as soon as the stream is finished.
lazy <- function(build) {
  val <- NULL; done <- FALSE
  list(get     = function() { if (!done) { val <<- build(); done <<- TRUE }; val },
       built   = function() done,
       release = function() { val <<- NULL; done <<- FALSE; invisible(gc()) })
}

# ==============================================================================
# Checkpoints: one file per independently runnable analysis ("unit")
# ------------------------------------------------------------------------------
#   <dir>/<kind>__<stream>__<id>__<cfgtag>.rds         a finished unit
#   <dir>/<kind>__<stream>__<id>__<cfgtag>.failed.rds  the last failure of a unit (never reused)
# Every file stores  key (what the result depends on)  +  value  +  metadata (time, warnings,
# R and package versions). A file is REUSED only if its key equals the key of the run that asks:
#   cfg   the settings that change this kind of result (see ckpt_cfg_part)
#   data  fingerprints of the input tables (file md5 of L1, or content summary for in-memory data)
#   spec  the model specification (predictors, k, adjusted, ...)
#   sample  the explicit person sample, if any
# Files are written to a temporary name and renamed, so a crash never leaves a half-written
# file that looks valid. A failed unit is recorded with its error message and re-run next time;
# it never overwrites or removes the result of any other unit.
# ==============================================================================

# Short md5 of a plain R object (lists/vectors only; used for settings, not for data.tables).
ckpt_hash <- function(x) {
  f <- tempfile(); on.exit(unlink(f), add = TRUE)
  writeBin(serialize(x, NULL, version = 2L), f)
  unname(tools::md5sum(f))
}

# The settings that change a result of this kind. Resource settings are deliberately absent.
# The person-sample options are not here either: they are part of the unit's own key, so results
# made with another sample size are re-run, while Cell 17's exploratory sample does not hide the
# full-data results of the other models.
ckpt_cfg_part <- function(cfg, kind) {
  core <- cfg[c("predictors", "own_predictors", "activity_predictors", "seed", "cohorts", "overlap")]
  model_part <- cfg$model[c("engine", "nAGQ", "min_events", "min_persons_with_events")]
  if (kind %in% c("gapB", "gapC")) return(c(core, list(gap_models = cfg$gap_models, boot = cfg$boot,
                                                      predict = cfg$predict, kind = kind)))
  extra <- list(onset = cfg$onset, model = model_part, primary_term = cfg$primary_term, kind = kind)
  if (kind == "sensitivity") extra$boot <- cfg$boot
  c(core, extra)
}

ckpt_key <- function(cfg, kind, data, spec = NULL, sample = NULL) {
  list(cfg = ckpt_cfg_part(cfg, kind), data = data, spec = spec, sample = sample)
}

ckpt_path <- function(dir, unit, key, failed = FALSE) {
  safe <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
  file.path(dir, sprintf("%s__%s__%s__%s%s.rds", unit$kind, safe(unit$stream), safe(unit$id),
                         substr(ckpt_hash(key$cfg), 1L, 8L), if (failed) ".failed" else ""))
}

pkg_versions <- function() {
  pk <- c("data.table", "lme4", "glmmTMB", "ordinal", "survival", "sandwich")
  vapply(pk, function(p) if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p)) else NA_character_, "")
}

# Which parts of two keys differ ("" = compatible).
ckpt_mismatch <- function(stored, wanted) {
  parts <- c("cfg", "data", "spec", "sample")
  paste(parts[!vapply(parts, function(p) identical(stored[[p]], wanted[[p]]), logical(1))], collapse = ", ")
}

# Atomic write: a temporary file in the same folder, then rename.
ckpt_write <- function(path, obj) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp-", Sys.getpid())
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(obj, tmp)
  if (!file.rename(tmp, path)) stop("could not move the finished checkpoint into place: ", path, call. = FALSE)
  invisible(path)
}

# Read a checkpoint and say whether it may be reused.
#   status: "missing" | "corrupt" | "stale" (reason = which key parts differ) | "valid"
ckpt_read <- function(path, key) {
  if (!file.exists(path)) return(list(status = "missing"))
  obj <- tryCatch(readRDS(path), error = function(e) NULL)
  if (!is.list(obj) || is.null(obj$key) || is.null(obj$unit)) {
    log_info("checkpoint %s is unreadable (interrupted write?) - it will be recomputed", basename(path))
    return(list(status = "corrupt"))
  }
  why <- ckpt_mismatch(obj$key, key)
  if (nzchar(why)) return(list(status = "stale", reason = why, obj = obj))
  list(status = "valid", obj = obj)
}

# Did a unit's value report a failed model? (status "failed" in a model-info table or row)
value_failed <- function(value) {
  st <- if (is.list(value) && !is.data.frame(value) && !is.null(value$info)) value$info$status
        else if (is.data.frame(value) && "status" %in% names(value)) value$status
        else character()
  length(st) > 0L && any(grepl("^failed", st))
}
value_error <- function(value) {
  tab <- if (is.list(value) && !is.data.frame(value) && !is.null(value$info)) value$info else value
  st <- tab$status; nt <- if ("note" %in% names(tab)) tab$note else rep("", length(st))
  paste(unique(nt[grepl("^failed", st)] [nzchar(nt[grepl("^failed", st)])]), collapse = " | ")
}

# Resume-or-compute one unit.
#   unit  list(kind, stream, id)           key  from ckpt_key()
#   compute()  returns the unit's value (a list or data.table); only called when no valid checkpoint exists
# Returns list(status = "resumed" | "computed" | "failed", value, error, path).
# Errors inside compute() are caught HERE (one place), written to the .failed file and logged; with
# cfg$resources$stop_on_error = TRUE they are re-thrown after being recorded.
run_unit <- function(dir, unit, key, compute, cfg = NULL) {
  path <- ckpt_path(dir, unit, key); fpath <- ckpt_path(dir, unit, key, failed = TRUE)
  label <- paste(unit$kind, unit$stream, unit$id, sep = " / ")
  r <- ckpt_read(path, key)
  if (r$status == "valid") {
    log_info("[skip] %s (checkpoint from %s)", label, format(r$obj$created, "%Y-%m-%d %H:%M"))
    return(list(status = "resumed", value = r$obj$value, error = NULL, path = path))
  }
  if (r$status == "stale") log_info("[stale: %s changed] %s", r$reason, label)
  log_info("[run] %s  (session %.1f GB, free %.1f GB)", label, rss_gb(), available_memory_gb())

  t0 <- Sys.time(); warns <- character()
  out <- tryCatch(withCallingHandlers(compute(), warning = function(w) {
                    warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") }),
                  error = function(e) e)
  fail <- function(msg) {
    log_info("[FAILED] %s: %s", label, msg)
    ckpt_write(fpath, list(unit = unit, key = key, error = msg, created = Sys.time(), warnings = unique(warns)))
    list(status = "failed", value = NULL, error = msg, path = fpath)
  }
  if (inherits(out, "error")) {
    res <- fail(conditionMessage(out))
    if (isTRUE(cfg$resources$stop_on_error)) stop(out)
    return(res)
  }
  if (value_failed(out)) return(fail(value_error(out)))

  ckpt_write(path, list(unit = unit, key = key, value = out, created = Sys.time(),
                        elapsed_sec = as.numeric(difftime(Sys.time(), t0, units = "secs")),
                        warnings = unique(warns), r_version = R.version.string, packages = pkg_versions()))
  unlink(fpath)
  list(status = "computed", value = out, error = NULL, path = path)
}

# All reusable results of one kind. `cfg_part` is ckpt_cfg_part(cfg, kind); `data_keys` a named list
# (by stream) of the CURRENT input fingerprints. Files made with other settings or other data are
# ignored. Returns list(values, units, failures).
collect_units <- function(dir, kind, cfg, data_keys = NULL) {
  cfg_now <- ckpt_cfg_part(cfg, kind)
  files <- list.files(dir, pattern = paste0("^", kind, "__.*\\.rds$"), full.names = TRUE)
  rows <- lapply(files, function(f) {
    o <- tryCatch(readRDS(f), error = function(e) NULL)
    if (!is.list(o) || is.null(o$key) || !identical(o$key$cfg, cfg_now)) return(NULL)
    if (!is.null(data_keys) && !identical(o$key$data, data_keys[[o$unit$stream]] %||% o$key$data)) return(NULL)
    list(obj = o, failed = grepl("\\.failed\\.rds$", f))
  })
  rows <- Filter(Negate(is.null), rows)
  ok <- Filter(function(r) !r$failed, rows); bad <- Filter(function(r) r$failed, rows)
  fail_tab <- rbindlist(lapply(bad, function(r) data.table(kind = r$obj$unit$kind, stream = r$obj$unit$stream,
                                                           id = r$obj$unit$id, error = r$obj$error,
                                                           time = format(r$obj$created))), fill = TRUE)
  list(values = lapply(ok, function(r) r$obj$value),
       units = rbindlist(lapply(ok, function(r) data.table(kind = r$obj$unit$kind, stream = r$obj$unit$stream,
                                                           id = r$obj$unit$id)), fill = TRUE),
       failures = fail_tab)
}

# One row per checkpoint file in `dir`: what is reusable now, what is stale, failed or unreadable.
# `data_keys` / `cfg` are optional; without them a file is only checked for readability.
checkpoint_report <- function(dir, cfg = NULL, data_keys = NULL) {
  files <- list.files(dir, pattern = "\\.rds$", full.names = TRUE)
  if (length(files) == 0L) return(data.table(file = character(), kind = character(), stream = character(),
                                             id = character(), state = character(), detail = character(),
                                             created = character(), elapsed_sec = numeric(), size_mb = numeric()))
  rbindlist(lapply(files, function(f) {
    o <- tryCatch(readRDS(f), error = function(e) NULL)
    base <- data.table(file = basename(f), size_mb = round(file.size(f) / 1e6, 2))
    if (!is.list(o) || is.null(o$unit)) return(cbind(base, kind = NA_character_, stream = NA_character_, id = NA_character_,
                                                      state = "corrupt", detail = "unreadable", created = NA_character_, elapsed_sec = NA_real_))
    state <- "valid"; detail <- ""
    if (grepl("\\.failed\\.rds$", f)) { state <- "failed"; detail <- o$error }
    else if (!is.null(cfg)) {
      cfg_ok  <- identical(o$key$cfg, ckpt_cfg_part(cfg, o$unit$kind))
      data_ok <- is.null(data_keys) || identical(o$key$data, data_keys[[o$unit$stream]] %||% o$key$data)
      if (!cfg_ok) { state <- "stale"; detail <- "settings changed" }
      else if (!data_ok) { state <- "stale"; detail <- "input data changed" }
    }
    if (length(o$warnings)) detail <- paste0(detail, if (nzchar(detail)) " | ", length(o$warnings), " warning(s): ", o$warnings[1L])
    cbind(base, kind = o$unit$kind, stream = o$unit$stream, id = o$unit$id, state = state, detail = detail,
          created = format(o$created), elapsed_sec = round(o$elapsed_sec %||% NA_real_, 1))
  }), fill = TRUE)
}

# Delete checkpoints that can never be reused again (stale, unreadable, failed records, leftover
# temporary files). dry_run = TRUE only lists them.
prune_checkpoints <- function(dir, cfg, data_keys = NULL, dry_run = TRUE) {
  rep <- checkpoint_report(dir, cfg, data_keys)
  drop <- rep[state %in% c("stale", "corrupt", "failed")]
  tmp <- list.files(dir, pattern = "\\.tmp-[0-9]+$", full.names = FALSE)
  if (!dry_run) unlink(file.path(dir, c(drop$file, tmp)))
  data.table(file = c(drop$file, tmp), removed = !dry_run)
}

# ---- fingerprints of input tables ------------------------------------------------------------------

# Content summary of a data.table: rows, people, date range and rounded column sums. Cheap (one
# pass) and changes whenever the data change; used for tables that are not read from a file.
fingerprint_table <- function(dt) {
  num <- names(dt)[vapply(dt, is.numeric, logical(1))]
  list(n_rows = nrow(dt), n_persons = if ("person_id" %in% names(dt)) uniqueN(dt$person_id) else NA_integer_,
       dates = if ("date" %in% names(dt)) as.character(range(dt$date)) else NULL,
       sums = vapply(num, function(cl) signif(sum(as.numeric(dt[[cl]]), na.rm = TRUE), 8L), numeric(1)))
}
