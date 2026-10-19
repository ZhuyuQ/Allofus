# ==============================================================================
# io.R
# ------------------------------------------------------------------------------
# All reading/writing of files lives here, so the analysis functions stay pure
# (data in -> data out).
#
#  * Raw BigQuery exports (L0)  : list_export_files(), read_export_files()
#  * Compact day tables (L1)    : write_day_table(), read_day_table()
#  * Cohort assignment / CCI    : read_final_pid(), read_cci()
#
# The reader and the file lister are injectable (`lister=`, `reader=`) so the
# same code is tested on local CSV files and run on gs:// paths in the workbench.
# ==============================================================================

# ---- L0: raw exports ---------------------------------------------------------

# Build the gs:// path (with wildcard when sharded) for one cohort x stream.
# Mirrors the layout used in the export notebooks:
#   <bucket>/bq_exports/<owner>/<date>/<prefix>_<stem>/<prefix>_<stem>[_*].csv
export_path <- function(cfg, cohort, stream) {
  prefix <- cfg$cohort_prefix[[cohort]]
  spec   <- if (stream == "meta") list(stem = cfg$meta_stem, sharded = TRUE)
            else if (stream == "device") list(stem = cfg$device_stem, sharded = TRUE)
            else cfg$streams[[stream]]
  folder <- paste0(prefix, "_", spec$stem)
  file   <- if (spec$sharded) paste0(folder, "_*.csv") else paste0(folder, ".csv")
  file.path(cfg$export$bucket, "bq_exports", cfg$export$owner,
            cfg$export$date, folder, file)
}

# `gsutil ls <pattern>` -> character vector of files. Stops if nothing matches.
gsutil_ls <- function(pattern) {
  out <- suppressWarnings(system2("gsutil", c("ls", shQuote(pattern)),
                                  stdout = TRUE, stderr = FALSE))
  status <- attr(out, "status")
  if (!is.null(status) || length(out) == 0L) {
    stop("no files found for: ", pattern, call. = FALSE)
  }
  out
}

# Read one file from GCS through `gsutil cat`.
gsutil_fread <- function(file, ...) {
  data.table::fread(cmd = paste("gsutil cat", shQuote(file)), ...)
}

# Read and stack all files matching `pattern`.
#  * every shard is read with the SAME column classes (the old reader re-inferred
#    types per shard, which can break bind_rows);
#  * `select` keeps only the needed columns (smaller memory).
read_export_files <- function(pattern, col_classes = NULL, select = NULL,
                              lister = gsutil_ls, reader = gsutil_fread) {
  files <- lister(pattern)
  log_info("reading %d file(s) for %s", length(files), basename(pattern))
  parts <- lapply(files, function(f) {
    reader(f, colClasses = col_classes, select = select,
           na.strings = c("", "NA", "NULL"))
  })
  data.table::rbindlist(parts, use.names = TRUE, fill = FALSE)
}

# ---- L1: compact day tables --------------------------------------------------

# Resolve "auto" to a concrete format.
resolve_format <- function(format) {
  if (format == "auto") {
    return(if (requireNamespace("arrow", quietly = TRUE)) "parquet" else "rds")
  }
  format
}

# Path of the day table for one stream: <dir>/day_<stream>.<ext>
day_table_path <- function(dir, stream, format) {
  file.path(dir, paste0("day_", stream, ".", if (format == "parquet") "parquet" else "rds"))
}

# Write a day table; the list of columns the ETL dropped is saved next to it (both formats).
write_day_table <- function(dt, dir, stream, format = "auto") {
  format <- resolve_format(format)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  path <- day_table_path(dir, stream, format)
  tmp <- paste0(path, ".tmp")                      # written under a temporary name, then renamed
  if (format == "parquet") arrow::write_parquet(as.data.frame(dt), tmp) else saveRDS(dt, tmp)
  file.rename(tmp, path)
  saveRDS(attr(dt, "dropped_columns"), file.path(dir, paste0("day_", stream, ".dropped.rds")))
  invisible(path)
}

# Read a day table back as a keyed data.table with Date dates. `cols` (optional) keeps only
# person_id, date and those columns - Parquet reads just them from disk; RDS has to read the whole
# file but the other columns are dropped at once.
read_day_table <- function(dir, stream, format = "auto", cols = NULL) {
  format <- resolve_format(format)
  path <- day_table_path(dir, stream, format)
  if (!file.exists(path)) stop("day table not found: ", path, call. = FALSE)
  keep <- if (is.null(cols)) NULL else unique(c("person_id", "date", cols))
  dt <- if (format == "parquet") {
    data.table::as.data.table(arrow::read_parquet(path, col_select = keep))
  } else {
    readRDS(path)
  }
  if (!is.null(keep) && format != "parquet") {
    missing_cols <- setdiff(keep, names(dt))
    if (length(missing_cols)) stop(sprintf("day table '%s' has no column(s): %s", stream, paste(missing_cols, collapse = ", ")), call. = FALSE)
    dt <- dt[, keep, with = FALSE]
  }
  dt[, date := as.Date(date)]
  data.table::setkey(dt, person_id, date)
  dt
}

read_dropped_columns <- function(dir, stream) {
  f <- file.path(dir, paste0("day_", stream, ".dropped.rds"))
  if (file.exists(f)) readRDS(f) else NULL
}

# ---- the day store: how the pipeline asks for day tables ------------------------------------------
#
# The pipeline never holds all streams at once. It asks the store for ONE stream (and only the
# columns it needs) and drops the table when the stream is finished.
#   open_day_store(cfg, dir)        reads from the L1 folder written by run_etl()
#   as_day_store(day_tables, persons)  wraps tables that are already in memory (demo, tests)
# store$get(stream, cols)   store$persons()   store$fingerprint(name)   store$streams
# fingerprint(): md5 of the L1 file (or a content summary for in-memory tables) - what a
# checkpoint is tied to, so changed input data can never be mixed with old results.

open_day_store <- function(cfg, dir) {
  format <- resolve_format(cfg$storage_format)
  streams <- names(cfg$streams)
  paths <- c(setNames(vapply(streams, function(s) day_table_path(dir, s, format), ""), streams),
             persons = file.path(dir, "persons.rds"))
  missing_files <- paths[!file.exists(paths)]
  if (length(missing_files)) stop("L1 file(s) not found - run the ETL (Cell 8) first: ",
                                  paste(basename(missing_files), collapse = ", "), call. = FALSE)
  fp <- list(); persons <- NULL
  structure(list(
    kind = "files", dir = dir, streams = streams, paths = paths,
    get = function(stream, cols = NULL) read_day_table(dir, stream, format, cols),
    persons = function() { if (is.null(persons)) persons <<- readRDS(paths[["persons"]]); persons },
    device_person = function() { f <- file.path(dir, "device_person.rds"); if (file.exists(f)) readRDS(f) else NULL },       # written by run_etl_device(); NULL if not built
    device_timeline = function() { f <- file.path(dir, "device_timeline.rds"); if (file.exists(f)) readRDS(f) else NULL },
    columns = function(stream) {                    # column names, from the ETL manifest when it has them
      mf <- file.path(dir, "manifest.tsv"); st <- stream
      cl <- if (file.exists(mf)) { m <- data.table::fread(mf); if ("columns" %in% names(m)) m$columns[m$stream == st] }
      if (length(cl) == 1L && nzchar(cl)) strsplit(cl, ",", fixed = TRUE)[[1L]] else names(read_day_table(dir, stream, format))
    },
    fingerprint = function(name) {
      if (is.null(fp[[name]])) fp[[name]] <<- list(md5 = unname(tools::md5sum(paths[[name]])), bytes = file.size(paths[[name]]))
      fp[[name]]
    }), class = "day_store")
}

as_day_store <- function(day_tables, persons = NULL, device = NULL) {     # device: list(person =, timeline =) from build_device_tables(), optional
  if (inherits(day_tables, "day_store")) return(day_tables)
  fp <- list()
  structure(list(
    kind = "memory", streams = names(day_tables),
    get = function(stream, cols = NULL) {
      d <- day_tables[[stream]]
      if (is.null(cols)) d else d[, unique(c("person_id", "date", cols)), with = FALSE]
    },
    persons = function() persons,
    device_person = function() device$person,
    device_timeline = function() device$timeline,
    columns = function(stream) names(day_tables[[stream]]),
    fingerprint = function(name) {
      if (is.null(fp[[name]])) fp[[name]] <<- fingerprint_table(if (name == "persons") persons else day_tables[[name]])
      fp[[name]]
    }), class = "day_store")
}

print.day_store <- function(x, ...) {
  cat("<day_store:", x$kind, ">  streams:", paste(x$streams, collapse = ", "), "\n")
  if (x$kind == "files") cat(sprintf("  %-9s %6.0f MB  %s\n", names(x$paths), file.size(x$paths) / 1e6, basename(x$paths)), sep = "")
  invisible(x)
}

# The input fingerprints a result of `stream` depends on: its own day table, the activity table
# (activity predictors enter every model), and the person table (covariates, cohorts).
data_key <- function(store, stream) {
  list(days = store$fingerprint(stream),
       act = if (stream != "activity") store$fingerprint("activity"),
       persons = store$fingerprint("persons"))
}

# ---- cohort assignment and CCI -----------------------------------------------

# final_pid.txt : tab-separated, columns person_id, cohort.
# One row per MEMBERSHIP: a person who belongs to two cohorts has two rows. (The exports are
# per cohort, so they use these memberships as they are; the ETL then resolves the overlap
# with resolve_cohort_overlap().)
read_final_pid <- function(path) {
  x <- data.table::fread(path, sep = "\t")
  assert_columns(x, c("person_id", "cohort"), "final_pid")
  unique(x[, .(person_id, cohort)])
}

# pre_final_cci.txt : tab-separated, columns person_id, cci_score (cohort optional).
read_cci <- function(path) {
  x <- data.table::fread(path, sep = "\t")
  assert_columns(x, c("person_id", "cci_score"), "CCI file")
  x[, .(cci_score = cci_score[1L]), by = person_id]
}
